#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(DESeq2)
  library(fs)
  library(ggplot2)
})

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L) y else x
}

project_root <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0L) {
    script_path <- normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/", mustWork = FALSE)
    return(dirname(dirname(script_path)))
  }
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

root_dir <- project_root()

unified_metadata_path <- file.path(root_dir, "test_data", "orf_unification_mouse_Mucosal_Immunity", "unified_orfs.metadata.tsv")
riboseqc_dir <- file.path(root_dir, "test_data", "riboseqc")
outdir <- file.path(root_dir, "test_results", "quant_analysis_CHX")
dir_create(outdir, recurse = TRUE)
padj_threshold <- 0.05
lfc_threshold <- 1

# sample annotations
sample_annotations <- data.table(
  sample_id = c("SRR7956050", "SRR7956051", "SRR7956052", "SRR7956053"),
  condition = c("CHX_NT", "CHX_LPS", "CHX_NT", "CHX_LPS"),
  replicate = c(1, 1, 2, 2)
)

parse_exon_blocks <- function(block_string) {
  if (is.na(block_string) || !nzchar(block_string)) {
    return(data.table(start = integer(), end = integer()))
  }
  parts <- strsplit(block_string, ",", fixed = TRUE)[[1]]
  coords <- tstrsplit(parts, "-", fixed = TRUE)
  data.table(
    start = as.integer(coords[[1]]),
    end = as.integer(coords[[2]])
  )
}

prepare_orf_catalog <- function(metadata_path) {
  cols <- c(
    "orf_id", "chrom", "strand", "start", "end", "length_aa", "exon_blocks",
    "gene_id", "transcript_id", "tools", "samples", "unique_psites", "pN"
  )
  orfs <- fread(metadata_path, sep = "\t", select = cols, showProgress = FALSE)
  orfs[, `:=`(
    start = as.integer(start),
    end = as.integer(end),
    length_aa = as.numeric(length_aa),
    unique_psites = as.numeric(unique_psites),
    pN = as.numeric(pN)
  )]
  orfs[, length_nt := vapply(strsplit(exon_blocks, ",", fixed = TRUE), function(parts) {
    sum(vapply(parts, function(x) {
      coords <- strsplit(x, "-", fixed = TRUE)[[1]]
      as.integer(coords[[2]]) - as.integer(coords[[1]]) + 1L
    }, integer(1)))
  }, integer(1))]
  orfs[length_aa >= 16 & pN >= 1]
}

make_block_index <- function(orf_catalog) {
  block_rows <- vector("list", nrow(orf_catalog))
  for (i in seq_len(nrow(orf_catalog))) {
    blocks <- parse_exon_blocks(orf_catalog$exon_blocks[[i]])
    if (nrow(blocks) == 0L) {
      next
    }
    block_rows[[i]] <- data.table(
      orf_id = orf_catalog$orf_id[[i]],
      chrom = orf_catalog$chrom[[i]],
      strand = orf_catalog$strand[[i]],
      block_start = blocks$start,
      block_end = blocks$end
    )
  }
  blocks_dt <- rbindlist(block_rows, use.names = TRUE, fill = TRUE)
  if (nrow(blocks_dt) == 0L) {
    return(blocks_dt)
  }
  blocks_dt[, `:=`(
    start1 = block_start,
    end1 = block_end,
    start0 = block_start - 1L,
    end0 = block_end
  )]
  setkey(blocks_dt, chrom, start0, end0)
  blocks_dt
}

read_uniq_bedgraph <- function(sample_id, strand, riboseqc_dir) {
  strand_suffix <- if (identical(strand, "+")) "plus" else "minus"
  path <- file.path(riboseqc_dir, sprintf("%s_P_sites_uniq_%s.bedgraph", sample_id, strand_suffix))
  if (!file.exists(path)) {
    stop(sprintf("Missing bedgraph file: %s", path), call. = FALSE)
  }
  dt <- fread(
    path,
    sep = "\t",
    header = FALSE,
    col.names = c("chrom", "bg_start0", "bg_end0", "value"),
    showProgress = FALSE
  )
  dt <- dt[!startsWith(chrom, "track") & !startsWith(chrom, "#")]
  dt[, `:=`(
    bg_start0 = as.integer(bg_start0),
    bg_end0 = as.integer(bg_end0),
    value = as.numeric(value)
  )]
  dt[value > 0]
}

count_sample_unique_psites <- function(sample_id, blocks_dt, riboseqc_dir) {
  counts <- vector("list", 2L)
  strands <- c("+", "-")
  for (j in seq_along(strands)) {
    strand_local <- strands[[j]]
    sample_blocks <- blocks_dt[strand == strand_local]
    if (nrow(sample_blocks) == 0L) {
      counts[[j]] <- data.table(orf_id = character(), unique_psites = numeric())
      next
    }
    bg <- read_uniq_bedgraph(sample_id, strand_local, riboseqc_dir)
    if (nrow(bg) == 0L) {
      counts[[j]] <- unique(sample_blocks[, .(orf_id)])[,
        .(orf_id, unique_psites = 0)]
      next
    }
    setkey(bg, chrom, bg_start0, bg_end0)
    ov <- foverlaps(bg, sample_blocks, by.x = c("chrom", "bg_start0", "bg_end0"), by.y = c("chrom", "start0", "end0"), nomatch = 0L)
    if (nrow(ov) == 0L) {
      counts[[j]] <- unique(sample_blocks[, .(orf_id)])[,
        .(orf_id, unique_psites = 0)]
      next
    }
    ov[, overlap_width := pmin(bg_end0, end0) - pmax(bg_start0, start0)]
    ov <- ov[overlap_width > 0]
    counts[[j]] <- ov[, .(unique_psites = sum(value * overlap_width)), by = orf_id]
  }

  sample_counts <- rbindlist(counts, use.names = TRUE, fill = TRUE)[, .(unique_psites = sum(unique_psites)), by = orf_id]
  sample_counts[]
}

build_count_and_metric_tables <- function(orf_catalog, sample_annotations, riboseqc_dir) {
  blocks_dt <- make_block_index(orf_catalog)
  per_sample <- vector("list", nrow(sample_annotations))

  for (i in seq_len(nrow(sample_annotations))) {
    sample_id <- sample_annotations$sample_id[[i]]
    message(sprintf("Counting sample-level unique P-sites for %s ...", sample_id))
    sample_counts <- count_sample_unique_psites(sample_id, blocks_dt, riboseqc_dir)
    sample_counts <- merge(
      orf_catalog[, .(orf_id, length_nt)],
      sample_counts,
      by = "orf_id",
      all.x = TRUE
    )
    sample_counts[is.na(unique_psites), unique_psites := 0]
    sample_counts[, `:=`(
      sample_id = sample_id,
      pN = unique_psites / length_nt
    )]
    per_sample[[i]] <- sample_counts[, .(orf_id, sample_id, unique_psites, pN)]
  }

  metrics_long <- rbindlist(per_sample, use.names = TRUE)
  counts_wide <- dcast(metrics_long, orf_id ~ sample_id, value.var = "unique_psites", fill = 0)
  pn_wide <- dcast(metrics_long, orf_id ~ sample_id, value.var = "pN", fill = 0)

  list(
    metrics_long = metrics_long,
    counts_wide = counts_wide,
    pn_wide = pn_wide
  )
}

run_deseq2 <- function(counts_wide, sample_annotations) {
  stopifnot(length(unique(counts_wide$orf_id)) == nrow(counts_wide))
  count_df <- as.data.frame(counts_wide)
  rownames(count_df) <- count_df$orf_id
  count_mat <- as.matrix(count_df[, sample_annotations$sample_id, drop = FALSE])
  storage.mode(count_mat) <- "integer"

  col_data <- as.data.frame(sample_annotations)
  rownames(col_data) <- col_data$sample_id
  col_data$condition <- factor(col_data$condition, levels = c("CHX_NT", "CHX_LPS"))

  dds <- DESeqDataSetFromMatrix(
    countData = count_mat,
    colData = col_data,
    design = ~ condition
  )
  dds <- dds[rowSums(counts(dds)) > 0, ]
  dds <- DESeq(dds, quiet = TRUE)
  res <- results(dds, contrast = c("condition", "CHX_LPS", "CHX_NT"))
  res_df <- as.data.frame(res)
  res_df$orf_id <- rownames(res_df)
  res_dt <- as.data.table(res_df)
  list(dds = dds, results = res_dt)
}

write_outputs <- function(orf_catalog, metrics_long, counts_wide, pn_wide, deseq_results, outdir) {
  analysis_orfs <- unique(metrics_long[, .(
    orf_id,
    max_sample_unique_psites = max(unique_psites, na.rm = TRUE),
    max_sample_pN = max(pN, na.rm = TRUE),
    mean_sample_pN = mean(pN, na.rm = TRUE)
  ), by = .(orf_id)])
  analysis_orfs <- analysis_orfs[, !duplicated(names(analysis_orfs)), with = FALSE]

  filtered_orfs <- merge(
    orf_catalog,
    analysis_orfs,
    by = "orf_id",
    all.x = FALSE
  )[max_sample_pN >= 1]

  filtered_ids <- filtered_orfs$orf_id
  filtered_counts <- counts_wide[orf_id %in% filtered_ids]
  filtered_metrics <- metrics_long[orf_id %in% filtered_ids]

  deseq_results <- as.data.table(deseq_results)
  deseq_results <- deseq_results[, !duplicated(names(deseq_results)), with = FALSE]

  deseq_dt <- merge(
    filtered_orfs,
    deseq_results[orf_id %in% filtered_ids],
    by = "orf_id",
    all.x = TRUE
  )
  deseq_dt[, significant := !is.na(padj) & padj < padj_threshold]
  deseq_dt[, abs_log2FoldChange := abs(log2FoldChange)]
  deseq_dt[, direction := fifelse(
    significant & log2FoldChange >= lfc_threshold, "Up in CHX_LPS",
    fifelse(significant & log2FoldChange <= -lfc_threshold, "Up in CHX_NT", "Not significant")
  )]
  setorder(deseq_dt, padj, -abs_log2FoldChange, -max_sample_pN)

  up_dt <- deseq_dt[direction == "Up in CHX_LPS"]
  down_dt <- deseq_dt[direction == "Up in CHX_NT"]

  fwrite(metrics_long, file.path(outdir, "orf_sample_metrics.long.tsv"), sep = "\t")
  fwrite(counts_wide, file.path(outdir, "orf_unique_psites_matrix.tsv"), sep = "\t")
  fwrite(pn_wide, file.path(outdir, "orf_pN_matrix.tsv"), sep = "\t")
  fwrite(filtered_metrics, file.path(outdir, "orf_sample_metrics.filtered.tsv"), sep = "\t")
  fwrite(filtered_counts, file.path(outdir, "orf_unique_psites_matrix.filtered.tsv"), sep = "\t")
  fwrite(deseq_dt, file.path(outdir, "deseq2_CHX_LPS_vs_CHX_NT.orf_results.tsv"), sep = "\t")
  fwrite(deseq_dt[significant == TRUE], file.path(outdir, "deseq2_CHX_LPS_vs_CHX_NT.significant_orfs.tsv"), sep = "\t")
  fwrite(up_dt, file.path(outdir, "deseq2_CHX_LPS_vs_CHX_NT.up_in_CHX_LPS.tsv"), sep = "\t")
  fwrite(down_dt, file.path(outdir, "deseq2_CHX_LPS_vs_CHX_NT.up_in_CHX_NT.tsv"), sep = "\t")

  summary_lines <- c(
    sprintf("Input ORFs (length_aa >= 16): %s", nrow(orf_catalog)),
    sprintf("ORFs passing sample-level pN >= 1: %s", length(filtered_ids)),
    sprintf("Samples: %s", paste(sample_annotations$sample_id, collapse = ", ")),
    sprintf("Contrast: %s", "CHX_LPS vs CHX_NT"),
    sprintf("Significant ORFs (padj < %.2f): %s", padj_threshold, nrow(deseq_dt[significant == TRUE])),
    sprintf("Up in CHX_LPS (padj < %.2f and log2FC >= %.1f): %s", padj_threshold, lfc_threshold, nrow(up_dt)),
    sprintf("Up in CHX_NT (padj < %.2f and log2FC <= -%.1f): %s", padj_threshold, lfc_threshold, nrow(down_dt))
  )
  writeLines(summary_lines, file.path(outdir, "analysis_summary.txt"))

  invisible(deseq_dt)
}

plot_volcano <- function(deseq_dt, outdir) {
  plot_dt <- copy(deseq_dt)
  plot_dt[, neg_log10_padj := -log10(pmax(padj, 1e-300))]
  plot_dt[!is.finite(neg_log10_padj), neg_log10_padj := NA_real_]
  top_labels <- plot_dt[direction != "Not significant" & !is.na(padj)][order(padj, -abs_log2FoldChange)][1:min(.N, 12)]

  p <- ggplot(plot_dt, aes(x = log2FoldChange, y = neg_log10_padj, color = direction)) +
    geom_point(alpha = 0.65, size = 1.2, na.rm = TRUE) +
    geom_vline(xintercept = c(-lfc_threshold, lfc_threshold), linetype = "dashed", linewidth = 0.3, color = "grey50") +
    geom_hline(yintercept = -log10(padj_threshold), linetype = "dashed", linewidth = 0.3, color = "grey50") +
    geom_text(
      data = top_labels,
      aes(x = log2FoldChange, y = neg_log10_padj, label = orf_id),
      inherit.aes = FALSE,
      nudge_y = 0.3,
      size = 2.5,
      check_overlap = TRUE
    ) +
    scale_color_manual(
      values = c(
        "Up in CHX_LPS" = "#c23b22",
        "Up in CHX_NT" = "#2166ac",
        "Not significant" = "grey70"
      )
    ) +
    labs(
      title = "CHX_LPS vs CHX_NT ORF Differential Translation",
      x = "log2 fold change",
      y = expression(-log[10](adjusted~p))
    ) +
    theme_minimal(base_size = 11) +
    theme(
      legend.title = element_blank(),
      plot.title = element_text(face = "bold")
    )

  ggsave(
    filename = file.path(outdir, "deseq2_CHX_LPS_vs_CHX_NT.volcano.png"),
    plot = p,
    width = 8.5,
    height = 6.5,
    dpi = 300,
    bg = "white"
  )
}

plot_significant_heatmap <- function(filtered_counts, sample_annotations, deseq_dt, outdir) {
  if (!requireNamespace("ComplexHeatmap", quietly = TRUE) || !requireNamespace("circlize", quietly = TRUE)) {
    warning("ComplexHeatmap/circlize not available; skipping significant ORF heatmap.")
    return(invisible(NULL))
  }

  heat_ids <- deseq_dt[direction != "Not significant" & !is.na(padj)][order(direction, padj, -abs_log2FoldChange), orf_id]
  if (length(heat_ids) == 0L) {
    return(invisible(NULL))
  }

  heat_dt <- filtered_counts[match(heat_ids, filtered_counts$orf_id)]
  heat_mat <- as.matrix(as.data.frame(heat_dt[, sample_annotations$sample_id, with = FALSE]))
  rownames(heat_mat) <- heat_dt$orf_id
  heat_mat <- log2(heat_mat + 1)
  heat_mat <- t(scale(t(heat_mat)))
  heat_mat[is.na(heat_mat)] <- 0

  row_groups <- deseq_dt[match(heat_ids, orf_id), direction]
  column_anno <- data.frame(condition = sample_annotations$condition)
  rownames(column_anno) <- sample_annotations$sample_id

  top_anno <- ComplexHeatmap::HeatmapAnnotation(
    df = column_anno,
    col = list(condition = c("CHX_NT" = "#2166ac", "CHX_LPS" = "#c23b22")),
    annotation_name_gp = grid::gpar(fontsize = 10, fontface = "bold")
  )

  ht <- ComplexHeatmap::Heatmap(
    heat_mat,
    name = "row z-score",
    col = circlize::colorRamp2(c(-2, 0, 2), c("#2166ac", "#f7f7f7", "#c23b22")),
    top_annotation = top_anno,
    cluster_columns = TRUE,
    cluster_rows = FALSE,
    show_row_names = FALSE,
    show_column_names = TRUE,
    column_names_gp = grid::gpar(fontsize = 10),
    row_split = factor(row_groups, levels = c("Up in CHX_LPS", "Up in CHX_NT")),
    row_title_gp = grid::gpar(fontsize = 10, fontface = "bold"),
    use_raster = TRUE,
    raster_quality = 2,
    column_title = sprintf("All Significant Differential ORFs (n = %s)", length(heat_ids)),
    column_title_gp = grid::gpar(fontsize = 12, fontface = "bold")
  )

  png(
    filename = file.path(outdir, "deseq2_CHX_LPS_vs_CHX_NT.all_significant_orfs_heatmap.png"),
    width = 2600,
    height = 3200,
    res = 300
  )
  ComplexHeatmap::draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
  dev.off()
}

message("Preparing ORF catalog ...")
orf_catalog <- prepare_orf_catalog(unified_metadata_path)
message(sprintf("Retained %s ORFs after metadata filtering (length_aa >= 16 and aggregated pN >= 1).", nrow(orf_catalog)))

message("Building sample-level count and pN matrices from RiboseQC unique P-site bedgraphs ...")
tables <- build_count_and_metric_tables(
  orf_catalog = orf_catalog,
  sample_annotations = sample_annotations,
  riboseqc_dir = riboseqc_dir
)

analysis_filter <- tables$metrics_long[, .(max_sample_pN = max(pN, na.rm = TRUE)), by = orf_id]
message(sprintf("Running DESeq2 on %s ORFs after sample-level recount ...", nrow(analysis_filter)))
deseq_input <- tables$counts_wide[orf_id %in% analysis_filter$orf_id]
deseq_fit <- run_deseq2(deseq_input, sample_annotations)

message("Writing outputs ...")
final_results <- write_outputs(
  orf_catalog = orf_catalog,
  metrics_long = tables$metrics_long,
  counts_wide = tables$counts_wide,
  pn_wide = tables$pn_wide,
  deseq_results = deseq_fit$results,
  outdir = outdir
)

plot_volcano(final_results, outdir)
plot_significant_heatmap(
  filtered_counts = tables$counts_wide[orf_id %in% final_results$orf_id],
  sample_annotations = sample_annotations,
  deseq_dt = final_results,
  outdir = outdir
)

message(sprintf("Analysis complete. Results written to %s", outdir))
message(sprintf("Final significant ORF count (padj < %.2f): %s", padj_threshold, nrow(final_results[significant == TRUE])))
