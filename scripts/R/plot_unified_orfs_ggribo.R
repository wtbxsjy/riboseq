#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

usage <- function() {
  cat(
    "Usage:\n",
    "  Rscript scripts/R/plot_unified_orfs_ggribo.R \\\n",
    "    --metadata test_data/orf_unification_mouse_Mucosal_Immunity/unified_orfs.metadata.tsv \\\n",
    "    --gtf test_data/orf_unification_mouse_Mucosal_Immunity/unified_orfs.gtf \\\n",
    "    --riboseqc-dir test_data/riboseqc \\\n",
    "    --outdir results/ggribo_plots \\\n",
    "    [--orf-ids ORF_1,ORF_2] [--orf-id-file ids.txt] [--top-n 10] \\\n",
    "    [--samples SRR7956050,SRR7956051] [--signal unique|all] \\\n",
    "    [--backend auto|ggribo|manual] [--min-unique-psites 10] \\\n",
    "    [--max-samples-per-orf 4] [--format png|pdf|both]\n",
    sep = ""
  )
}

parse_args <- function(args) {
  if (length(args) == 0L || any(args %in% c("-h", "--help"))) {
    usage()
    quit(status = 0L)
  }

  out <- list(
    metadata = NULL,
    gtf = NULL,
    riboseqc_dir = NULL,
    outdir = NULL,
    orf_ids = NULL,
    orf_id_file = NULL,
    top_n = 10L,
    samples = NULL,
    signal = "unique",
    min_unique_psites = 10L,
    max_samples_per_orf = 4L,
    format = "both",
    backend = "auto",
    width = 12,
    height_base = 2.8,
    dpi = 150
  )

  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) {
      stop(sprintf("Unexpected argument: %s", key), call. = FALSE)
    }
    if (i == length(args)) {
      stop(sprintf("Missing value for argument: %s", key), call. = FALSE)
    }
    value <- args[[i + 1L]]
    switch(
      key,
      "--metadata" = out$metadata <- value,
      "--gtf" = out$gtf <- value,
      "--riboseqc-dir" = out$riboseqc_dir <- value,
      "--outdir" = out$outdir <- value,
      "--orf-ids" = out$orf_ids <- value,
      "--orf-id-file" = out$orf_id_file <- value,
      "--top-n" = out$top_n <- as.integer(value),
      "--samples" = out$samples <- value,
      "--signal" = out$signal <- value,
      "--min-unique-psites" = out$min_unique_psites <- as.integer(value),
      "--max-samples-per-orf" = out$max_samples_per_orf <- as.integer(value),
      "--format" = out$format <- value,
      "--backend" = out$backend <- value,
      "--width" = out$width <- as.numeric(value),
      "--height-base" = out$height_base <- as.numeric(value),
      "--dpi" = out$dpi <- as.integer(value),
      stop(sprintf("Unknown argument: %s", key), call. = FALSE)
    )
    i <- i + 2L
  }

  required <- c("metadata", "riboseqc_dir", "outdir")
  missing <- required[vapply(required, function(x) is.null(out[[x]]) || !nzchar(out[[x]]), logical(1))]
  if (length(missing) > 0L) {
    stop(sprintf("Missing required arguments: %s", paste(missing, collapse = ", ")), call. = FALSE)
  }

  if (!out$signal %in% c("unique", "all")) {
    stop("--signal must be 'unique' or 'all'", call. = FALSE)
  }
  if (!out$format %in% c("png", "pdf", "both")) {
    stop("--format must be 'png', 'pdf', or 'both'", call. = FALSE)
  }
  if (!out$backend %in% c("auto", "ggribo", "manual")) {
    stop("--backend must be 'auto', 'ggribo', or 'manual'", call. = FALSE)
  }

  out
}

split_csv <- function(x) {
  if (is.null(x) || !nzchar(x) || identical(x, "NA")) {
    return(character())
  }
  trimws(unlist(strsplit(x, ",", fixed = TRUE)))
}

parse_score_map <- function(x) {
  if (is.na(x) || !nzchar(x) || identical(x, "NA")) {
    return(character())
  }
  parts <- trimws(unlist(strsplit(x, ",", fixed = TRUE)))
  parts[nzchar(parts)]
}

clean_attrs <- function(x) {
  x <- sub('^"', "", x)
  x <- sub('"$', "", x)
  gsub('""', '"', x, fixed = TRUE)
}

extract_attr <- function(attrs, key) {
  m <- regexec(sprintf('%s "([^"]+)"', key), attrs)
  hit <- regmatches(attrs, m)[[1]]
  if (length(hit) < 2L) {
    return(NA_character_)
  }
  hit[2]
}

parse_exon_blocks <- function(blocks_string, strand) {
  blocks <- rbindlist(lapply(strsplit(blocks_string, ",", fixed = TRUE)[[1]], function(chunk) {
    bounds <- as.integer(strsplit(chunk, "-", fixed = TRUE)[[1]])
    data.table(start = bounds[1], end = bounds[2])
  }))
  setorder(blocks, start, end)
  blocks[, block_index := .I]
  blocks[, width := end - start + 1L]

  tx_blocks <- copy(blocks)
  if (strand == "-") {
    tx_blocks <- tx_blocks[.N:1]
  }
  tx_blocks[, tx_index := .I]
  tx_blocks[, tx_start := cumsum(shift(width, fill = 0L)) + 1L]
  tx_blocks[, tx_end := tx_start + width - 1L]
  tx_blocks[]
}

map_positions_to_tx <- function(pos, tx_blocks, strand) {
  local <- integer(length(pos))
  local[] <- NA_integer_
  for (i in seq_len(nrow(tx_blocks))) {
    block <- tx_blocks[i]
    idx <- pos >= block$start & pos <= block$end
    if (!any(idx)) {
      next
    }
    if (strand == "+") {
      local[idx] <- block$tx_start + (pos[idx] - block$start)
    } else {
      local[idx] <- block$tx_start + (block$end - pos[idx])
    }
  }
  local
}

detect_available_samples <- function(riboseqc_dir, signal) {
  suffix <- if (signal == "unique") "_ggribo.tsv" else "_P_sites_plus.bedgraph"
  files <- list.files(riboseqc_dir, pattern = if (signal == "unique") "_ggribo\\.tsv$" else "_P_sites_plus\\.bedgraph$", full.names = FALSE)
  sub(suffix, "", files, fixed = TRUE)
}

read_signal_track <- function(sample_id, riboseqc_dir, signal, cache_env) {
  key <- paste(sample_id, signal, sep = "::")
  if (exists(key, envir = cache_env, inherits = FALSE)) {
    return(get(key, envir = cache_env, inherits = FALSE))
  }

  dt <- if (signal == "unique") {
    tsv_path <- file.path(riboseqc_dir, sprintf("%s_ggribo.tsv", sample_id))
    if (file.exists(tsv_path)) {
      x <- fread(tsv_path, header = FALSE, col.names = c("count", "chrom", "pos", "strand"))
      x[, pos := as.integer(pos)]
      x[, count := as.numeric(count)]
      x[]
    } else {
      plus_path <- file.path(riboseqc_dir, sprintf("%s_P_sites_uniq_plus.bedgraph", sample_id))
      minus_path <- file.path(riboseqc_dir, sprintf("%s_P_sites_uniq_minus.bedgraph", sample_id))
      read_bedgraph_pair(plus_path, minus_path)
    }
  } else {
    plus_path <- file.path(riboseqc_dir, sprintf("%s_P_sites_plus.bedgraph", sample_id))
    minus_path <- file.path(riboseqc_dir, sprintf("%s_P_sites_minus.bedgraph", sample_id))
    read_bedgraph_pair(plus_path, minus_path)
  }

  assign(key, dt, envir = cache_env)
  dt
}

read_bedgraph_pair <- function(plus_path, minus_path) {
  read_one <- function(path, strand) {
    if (!file.exists(path)) {
      return(data.table(count = numeric(), chrom = character(), pos = integer(), strand = character()))
    }
    x <- fread(path, header = FALSE, col.names = c("chrom", "start0", "end1", "count"))
    x[, `:=`(pos = as.integer(end1), strand = strand, count = as.numeric(count))]
    x[, .(count, chrom, pos, strand)]
  }
  rbindlist(list(read_one(plus_path, "+"), read_one(minus_path, "-")))
}

intersect_samples <- function(row_samples, available_samples, requested_samples = NULL) {
  x <- intersect(split_csv(row_samples), available_samples)
  if (!is.null(requested_samples) && length(requested_samples) > 0L) {
    x <- intersect(x, requested_samples)
  }
  unique(x)
}

make_orf_summary_label <- function(row, selected_samples) {
  lines <- c(
    sprintf("%s | %s:%s-%s (%s)", row$orf_id, row$chrom, format(row$start, big.mark = ","), format(row$end, big.mark = ","), row$strand),
    sprintf("tools: %s", row$tools),
    sprintf("samples shown: %s", paste(selected_samples, collapse = ", ")),
    sprintf("unique_psites=%s | pN=%s", row$unique_psites, row$pN)
  )

  tool_scores <- parse_score_map(row$tool_scores)
  if (length(tool_scores) > 0L) {
    lines <- c(lines, sprintf("tool_scores: %s", paste(tool_scores, collapse = "; ")))
  }
  tool_pvalues <- parse_score_map(row$tool_pvalues)
  if (length(tool_pvalues) > 0L) {
    lines <- c(lines, sprintf("tool_pvalues: %s", paste(tool_pvalues, collapse = "; ")))
  }
  paste(lines, collapse = "\n")
}

build_sample_plot <- function(sample_dt, sample_id, tx_blocks, total_orf_len, signal) {
  if (nrow(sample_dt) == 0L) {
    base <- data.table(local_pos = seq_len(total_orf_len), count = 0, frame = factor("frame0", levels = c("frame0", "frame1", "frame2")))
  } else {
    setorder(sample_dt, local_pos)
    sample_dt[, frame := factor(paste0("frame", (local_pos - 1L) %% 3L), levels = c("frame0", "frame1", "frame2"))]
    base <- sample_dt
  }

  exon_boundaries <- tx_blocks$tx_end[-nrow(tx_blocks)]
  sample_total <- if (nrow(sample_dt) > 0L) sum(sample_dt$count, na.rm = TRUE) else 0

  ggplot(base, aes(x = local_pos, y = count, fill = frame)) +
    geom_col(width = 1) +
    geom_vline(xintercept = exon_boundaries + 0.5, linetype = "dashed", color = "grey70", linewidth = 0.3) +
    scale_fill_manual(values = c(frame0 = "#1b9e77", frame1 = "#d95f02", frame2 = "#7570b3"), drop = FALSE) +
    labs(
      title = sample_id,
      subtitle = sprintf("track_sum=%s | signal=%s", format(round(sample_total, 2), trim = TRUE), signal),
      x = "ORF transcript coordinate",
      y = "P-site count",
      fill = "Frame"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 11),
      plot.subtitle = element_text(size = 9),
      legend.position = "top",
      panel.grid.minor = element_blank()
    )
}

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L || all(is.na(x))) y else x
}

make_annotation_plot <- function(tx_blocks, label_text, total_orf_len) {
  exon_dt <- copy(tx_blocks)
  exon_dt[, y := 1]
  boundary_dt <- data.table(x = tx_blocks$tx_end[-nrow(tx_blocks)])

  ggplot() +
    geom_rect(
      data = exon_dt,
      aes(xmin = tx_start, xmax = tx_end, ymin = 0.15, ymax = 0.45),
      fill = "#4c78a8",
      color = "#2f4b6c"
    ) +
    geom_vline(data = boundary_dt, aes(xintercept = x + 0.5), linetype = "dashed", color = "grey60", linewidth = 0.3) +
    annotate("text", x = 1, y = 2.05, label = label_text, hjust = 0, vjust = 1, size = 3.2, family = "mono") +
    scale_x_continuous(limits = c(1, total_orf_len), expand = c(0, 0)) +
    coord_cartesian(ylim = c(0, 2.15), clip = "off") +
    labs(x = NULL, y = NULL) +
    theme_void(base_size = 11) +
    theme(plot.margin = margin(5.5, 5.5, 10, 5.5))
}

sanitize_filename <- function(x) {
  gsub("[^A-Za-z0-9._-]", "_", x)
}

load_gtf_cache <- function(gtf_path, cache_env) {
  key <- paste0("gtf::", normalizePath(gtf_path, winslash = "/", mustWork = FALSE))
  if (exists(key, envir = cache_env, inherits = FALSE)) {
    return(get(key, envir = cache_env, inherits = FALSE))
  }
  gtf <- fread(gtf_path, sep = "\t", header = FALSE, quote = "")
  setnames(gtf, paste0("V", 1:9))
  gtf[, V9 := clean_attrs(V9)]
  assign(key, gtf, envir = cache_env)
  gtf
}

augment_gtf_for_txdb <- function(gtf_path, row, outdir, cache_env) {
  cache_key <- paste0("augmented_gtf::", row$orf_id)
  if (exists(cache_key, envir = cache_env, inherits = FALSE)) {
    return(get(cache_key, envir = cache_env, inherits = FALSE))
  }

  gtf <- load_gtf_cache(gtf_path, cache_env)
  subset_gtf <- gtf[grepl(sprintf('orf_id "%s"', row$orf_id), V9, fixed = TRUE)]
  if (nrow(subset_gtf) == 0L) {
    stop(sprintf("ORF %s not found in %s", row$orf_id, gtf_path), call. = FALSE)
  }

  attrs <- subset_gtf$V9[1]
  gene_id <- extract_attr(attrs, "gene_id")
  tx_id <- extract_attr(attrs, "transcript_id")
  chrom <- subset_gtf$V1[1]
  source <- subset_gtf$V2[1]
  strand <- subset_gtf$V7[1]
  start <- min(subset_gtf$V4)
  end <- max(subset_gtf$V5)

  gene_row <- data.table(
    V1 = chrom, V2 = source, V3 = "gene", V4 = start, V5 = end,
    V6 = ".", V7 = strand, V8 = ".", V9 = sprintf('gene_id "%s";', gene_id)
  )
  tx_row <- data.table(
    V1 = chrom, V2 = source, V3 = "transcript", V4 = start, V5 = end,
    V6 = ".", V7 = strand, V8 = ".", V9 = sprintf('gene_id "%s"; transcript_id "%s";', gene_id, tx_id)
  )

  augmented <- rbindlist(list(gene_row, tx_row, subset_gtf), use.names = TRUE)
  out_path <- file.path(outdir, sprintf("%s.ggRibo_input.gtf", sanitize_filename(row$orf_id)))
  fwrite(augmented, out_path, sep = "\t", col.names = FALSE, quote = FALSE)
  out <- list(path = out_path, gene_id = gene_id, tx_id = tx_id)
  assign(cache_key, out, envir = cache_env)
  out
}

build_range_info_ggribo <- function(annotation, orf_id, cache_env) {
  cache_key <- paste0("range_info::", orf_id)
  if (exists(cache_key, envir = cache_env, inherits = FALSE)) {
    return(get(cache_key, envir = cache_env, inherits = FALSE))
  }
  Range_info <- get("Range_info", envir = asNamespace("ggRibo"))
  txdb <- suppressWarnings(txdbmaker::makeTxDbFromGFF(file = annotation, format = "gtf"))
  exonsByTx <- GenomicFeatures::exonsBy(txdb, by = "tx", use.names = TRUE)
  txByGene <- GenomicFeatures::transcriptsBy(txdb, by = "gene")
  cdsByTx <- GenomicFeatures::cdsBy(txdb, by = "tx", use.names = TRUE)
  fiveUTR <- GenomicFeatures::fiveUTRsByTranscript(txdb, use.names = TRUE)
  threeUTR <- GenomicFeatures::threeUTRsByTranscript(txdb, use.names = TRUE)
  tx_to_gene <- AnnotationDbi::select(
    txdb,
    keys = AnnotationDbi::keys(txdb, keytype = "TXNAME"),
    columns = c("TXNAME", "GENEID"),
    keytype = "TXNAME"
  )
  colnames(tx_to_gene) <- c("tx_id", "gene_id")
  tx_to_gene <- tx_to_gene[order(tx_to_gene$tx_id), ]

  out <- Range_info$new(
    exonsByTx = exonsByTx,
    txByGene = txByGene,
    cdsByTx = cdsByTx,
    fiveUTR = fiveUTR,
    threeUTR = threeUTR,
    tx_to_gene = tx_to_gene
  )
  assign(cache_key, out, envir = cache_env)
  out
}

build_ggribo_inputs <- function(samples, riboseqc_dir, signal, cache_env) {
  if (signal == "unique") {
    lapply(samples, function(sample_id) {
      df <- copy(read_signal_track(sample_id, riboseqc_dir, signal, cache_env))
      setnames(df, c("count", "chrom", "pos", "strand"), c("count", "chr", "position", "strand"), skip_absent = TRUE)
      list(type = "tabular", data = as.data.frame(df[, .(count, chr, position, strand)]))
    })
  } else {
    lapply(samples, function(sample_id) {
      list(
        type = "bedgraph",
        plus = file.path(riboseqc_dir, sprintf("%s_P_sites_plus.bedgraph", sample_id)),
        minus = file.path(riboseqc_dir, sprintf("%s_P_sites_minus.bedgraph", sample_id))
      )
    })
  }
}

build_orf_plot_ggribo <- function(row, gtf_path, riboseqc_dir, available_samples, requested_samples, signal, max_samples_per_orf, outdir, cache_env) {
  selected_samples <- intersect_samples(row$samples, available_samples, requested_samples)
  if (length(selected_samples) == 0L) {
    return(NULL)
  }
  selected_samples <- head(selected_samples, max_samples_per_orf)

  annotation_info <- augment_gtf_for_txdb(gtf_path, row, outdir, cache_env)
  Riboseq_inputs <- build_ggribo_inputs(selected_samples, riboseqc_dir, signal, cache_env)
  GRangeInfo <- build_range_info_ggribo(annotation_info$path, row$orf_id, cache_env)

  ggRibo::ggRibo_tx(
    gene_id = annotation_info$gene_id,
    tx_id = annotation_info$tx_id,
    NAME = sprintf("%s | tools=%s | unique_psites=%s | pN=%s", row$orf_id, row$tools, row$unique_psites, row$pN),
    RNAseq = NULL,
    Riboseq = Riboseq_inputs,
    SampleNames = selected_samples,
    GRangeInfo = GRangeInfo,
    Y_scale = "each",
    plot_ORF_ranges = FALSE,
    show_seq = FALSE,
    title_font_size = 10,
    sample_label_font_size = 4
  )
}

build_orf_plot <- function(row, riboseqc_dir, available_samples, requested_samples, signal, max_samples_per_orf, cache_env) {
  selected_samples <- intersect_samples(row$samples, available_samples, requested_samples)
  if (length(selected_samples) == 0L) {
    return(NULL)
  }
  selected_samples <- head(selected_samples, max_samples_per_orf)

  tx_blocks <- parse_exon_blocks(row$exon_blocks, row$strand)
  total_orf_len <- sum(tx_blocks$width)

  sample_plots <- lapply(selected_samples, function(sample_id) {
    track <- copy(read_signal_track(sample_id, riboseqc_dir, signal, cache_env))
    track <- track[chrom == row$chrom & strand == row$strand]
    if (nrow(track) > 0L) {
      keep <- Reduce(`|`, lapply(seq_len(nrow(tx_blocks)), function(i) track$pos >= tx_blocks$start[i] & track$pos <= tx_blocks$end[i]))
      track <- track[keep]
    }
    if (nrow(track) > 0L) {
      track[, local_pos := map_positions_to_tx(pos, tx_blocks, row$strand)]
      track <- track[!is.na(local_pos), .(count = sum(count)), by = .(local_pos)]
    } else {
      track <- data.table(local_pos = integer(), count = numeric())
    }
    build_sample_plot(track, sample_id, tx_blocks, total_orf_len, signal)
  })

  info_label <- make_orf_summary_label(row, selected_samples)
  annotation_plot <- make_annotation_plot(tx_blocks, info_label, total_orf_len)

  wrap_plots(sample_plots, ncol = 1, guides = "collect") / annotation_plot +
    plot_layout(heights = c(rep(1, length(sample_plots)), 0.9))
}

select_orfs <- function(metadata_path, available_samples, requested_orf_ids, requested_samples, top_n, min_unique_psites) {
  dt <- fread(
    metadata_path,
    sep = "\t",
    select = c("orf_id", "chrom", "strand", "start", "end", "exon_blocks", "gene_id", "transcript_id", "tools", "samples", "tool_scores", "tool_pvalues", "unique_psites", "pN")
  )

  if (!is.null(requested_orf_ids) && length(requested_orf_ids) > 0L) {
    dt <- dt[orf_id %in% requested_orf_ids]
  } else {
    dt <- dt[unique_psites >= min_unique_psites]
  }

  dt[, matched_samples := lapply(samples, intersect_samples, available_samples = available_samples, requested_samples = requested_samples)]
  dt[, matched_sample_count := lengths(matched_samples)]
  dt <- dt[matched_sample_count > 0L]

  if (nrow(dt) == 0L) {
    return(dt)
  }

  if (is.null(requested_orf_ids) || length(requested_orf_ids) == 0L) {
    setorder(dt, -unique_psites, -matched_sample_count, orf_id)
    dt <- dt[seq_len(min(top_n, .N))]
  } else {
    dt[, requested_order := match(orf_id, requested_orf_ids)]
    setorder(dt, requested_order, orf_id)
    dt[, requested_order := NULL]
  }

  dt[]
}

main <- function() {
  opt <- parse_args(commandArgs(trailingOnly = TRUE))
  dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

  available_samples <- detect_available_samples(opt$riboseqc_dir, opt$signal)
  if (length(available_samples) == 0L) {
    stop("No compatible RiboseQC signal files found.", call. = FALSE)
  }

  requested_samples <- split_csv(opt$samples)
  requested_orf_ids <- split_csv(opt$orf_ids)
  if (!is.null(opt$orf_id_file) && nzchar(opt$orf_id_file)) {
    requested_orf_ids <- unique(c(requested_orf_ids, trimws(readLines(opt$orf_id_file, warn = FALSE))))
    requested_orf_ids <- requested_orf_ids[nzchar(requested_orf_ids)]
  }
  if (length(requested_orf_ids) == 0L) {
    requested_orf_ids <- NULL
  }
  if (length(requested_samples) == 0L) {
    requested_samples <- NULL
  }

  ggribo_available <- requireNamespace("ggRibo", quietly = TRUE) &&
    requireNamespace("txdbmaker", quietly = TRUE) &&
    requireNamespace("GenomicFeatures", quietly = TRUE) &&
    requireNamespace("AnnotationDbi", quietly = TRUE)
  backend <- opt$backend
  if (backend == "auto") {
    backend <- if (ggribo_available && !is.null(opt$gtf) && nzchar(opt$gtf)) "ggribo" else "manual"
  }
  if (backend == "ggribo" && !ggribo_available) {
    stop("Requested --backend ggribo but ggRibo/txdbmaker/GenomicFeatures/AnnotationDbi are not available.", call. = FALSE)
  }
  if (backend == "ggribo" && (is.null(opt$gtf) || !nzchar(opt$gtf))) {
    stop("Requested --backend ggribo but --gtf was not provided.", call. = FALSE)
  }

  message(sprintf("Detected %d available RiboseQC samples: %s", length(available_samples), paste(sort(available_samples), collapse = ", ")))
  message(sprintf("Using plotting backend: %s", backend))

  selected <- select_orfs(
    metadata_path = opt$metadata,
    available_samples = available_samples,
    requested_orf_ids = requested_orf_ids,
    requested_samples = requested_samples,
    top_n = opt$top_n,
    min_unique_psites = opt$min_unique_psites
  )

  if (nrow(selected) == 0L) {
    stop("No ORFs matched the requested filters and available RiboseQC samples.", call. = FALSE)
  }

  fwrite(
    selected[, .(orf_id, chrom, strand, start, end, tools, samples, unique_psites, pN, matched_sample_count)],
    file.path(opt$outdir, "selected_orfs.tsv"),
    sep = "\t"
  )

  cache_env <- new.env(parent = emptyenv())
  pdf_path <- file.path(opt$outdir, "all_selected_orfs.pdf")
  if (opt$format %in% c("pdf", "both")) {
    pdf(pdf_path, width = opt$width, height = opt$height_base * max(2, opt$max_samples_per_orf + 1))
    on.exit(dev.off(), add = TRUE)
  }

  for (i in seq_len(nrow(selected))) {
    row <- selected[i]
    message(sprintf("[%d/%d] plotting %s", i, nrow(selected), row$orf_id))
    p <- if (backend == "ggribo") {
      build_orf_plot_ggribo(
        row = row,
        gtf_path = opt$gtf,
        riboseqc_dir = opt$riboseqc_dir,
        available_samples = available_samples,
        requested_samples = requested_samples,
        signal = opt$signal,
        max_samples_per_orf = opt$max_samples_per_orf,
        outdir = opt$outdir,
        cache_env = cache_env
      )
    } else {
      build_orf_plot(
        row = row,
        riboseqc_dir = opt$riboseqc_dir,
        available_samples = available_samples,
        requested_samples = requested_samples,
        signal = opt$signal,
        max_samples_per_orf = opt$max_samples_per_orf,
        cache_env = cache_env
      )
    }
    if (is.null(p)) {
      next
    }

    base_name <- sanitize_filename(row$orf_id)
    if (opt$format %in% c("png", "both")) {
      png_path <- file.path(opt$outdir, sprintf("%s.png", base_name))
      ggsave(
        filename = png_path,
        plot = p,
        width = opt$width,
        height = opt$height_base * (length(intersect_samples(row$samples, available_samples, requested_samples)) + 0.9),
        dpi = opt$dpi,
        limitsize = FALSE
      )
    }
    if (opt$format %in% c("pdf", "both")) {
      print(p)
    }
  }

  message(sprintf("Finished. Outputs written to %s", opt$outdir))
}

main()
