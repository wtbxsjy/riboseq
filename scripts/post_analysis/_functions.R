# ─── Shared functions for ORF post-analysis ──────────────────────────
# Source this from Quarto or standalone scripts:
#   source(here::here("scripts/post_analysis/_functions.R"))
#
# Depends on: dplyr, ggplot2, tidyr, yaml, data.table, scales
# Optional: ggRibo (for ggribo section)

# ── Internal: read TSV with best available backend ────────────────────
.read_tsv <- function(path, ...) {
  if (requireNamespace("data.table", quietly = TRUE)) {
    as_tibble(data.table::fread(path, sep = "\t", header = TRUE,
                                stringsAsFactors = FALSE, ...))
  } else {
    as_tibble(utils::read.delim(path, sep = "\t", header = TRUE,
                                stringsAsFactors = FALSE, ...))
  }
}

# ── 0. YAML config loader ─────────────────────────────────────────────

load_config <- function(config_path) {
  cfg <- yaml::read_yaml(config_path)

  # Expand {result_dir} template
  rd <- cfg$input$result_dir
  for (key in setdiff(names(cfg$input), "result_dir")) {
    cfg$input[[key]] <- glue::glue(cfg$input[[key]], result_dir = rd)
  }

  # Ensure output directory
  dir.create(cfg$input$output_dir, showWarnings = FALSE, recursive = TRUE)

  cfg
}


# ── 1. Data loading + ID mapping ──────────────────────────────────────

#' Load all pipeline outputs and join by (chrom, start, end, strand)
#'
#' Returns a list with elements: metadata, confidence, expression,
#' purity, gencode, merged (the unified cross-walk).
load_pipeline_data <- function(cfg) {
  message("Loading metadata ...")
  meta <- .read_tsv(cfg$input$metadata) |>
    select(orf_id, chrom, strand, start, end, length_aa,
           tools, samples, total_reads, total_psites, pN)

  # Map (chrom, strand) to character for consistent joining
  meta <- meta |>
    mutate(across(c(chrom, strand), as.character),
           coord_key = paste(chrom, start, end, strand, sep = ":"))

  message(sprintf("  %s ORFs in metadata", scales::comma(nrow(meta))))

  # ── Confidence (OCS) ──
  message("Loading confidence ...")
  conf <- .read_tsv(cfg$input$confidence) |>
    select(orf_id, tier, ocs, s_translation, s_agreement,
           s_coverage, s_periodicity, s_readlevel,
           detecting_tools, n_detecting)
  message(sprintf("  %s ORFs with OCS", scales::comma(nrow(conf))))

  # ── Expression ──
  message("Loading expression ...")
  expr <- .read_tsv(cfg$input$expression)
  # Compute per-ORF aggregates across samples
  read_cols <- grep("_reads$", names(expr), value = TRUE)
  pN_cols   <- grep("_pN$", names(expr), value = TRUE)

  expr_agg <- expr |>
    mutate(
      total_reads_expr = rowSums(pick(all_of(read_cols)), na.rm = TRUE),
      max_pN     = apply(pick(all_of(pN_cols)), 1, max, na.rm = TRUE),
      n_samples  = apply(pick(all_of(read_cols)), 1, function(x) sum(x > 0, na.rm = TRUE)),
      mean_ratio = rowMeans(pick(all_of(pN_cols)), na.rm = TRUE)
    ) |>
    select(orf_id, total_reads_expr, max_pN, n_samples, mean_ratio)
  message(sprintf("  %s ORFs with expression", scales::comma(nrow(expr_agg))))

  # ── GENCODE classification ──
  message("Loading GENCODE classification ...")
  gencode_path <- cfg$input$gencode
  if (grepl("\\.gz$", gencode_path)) {
    gencode <- .read_tsv(gencode_path)
  } else {
    gencode <- .read_tsv(gencode_path)
  }
  gencode <- gencode |>
    rename_with(~ case_when(
      . == "chrm"   ~ "chrom",
      . == "starts" ~ "start",
      . == "ends"   ~ "end",
      TRUE ~ .
    )) |>
    mutate(
      across(c(chrom, strand), as.character),
      coord_key = paste(chrom, start, end, strand, sep = ":")
    ) |>
    select(coord_key, orf_biotype, gene_biotype, trans, gene, gene_name,
           phaseI_biotype, pep, orf_length)
  message(sprintf("  %s ORFs with GENCODE classification", scales::comma(nrow(gencode))))

  # ── P-site purity ──
  purity_file <- file.path(cfg$input$output_dir, "psite_purity.tsv")
  if (file.exists(purity_file)) {
    message("Loading P-site purity ...")
    purity <- .read_tsv(purity_file) |>
      rename(orf_id = orf_id_raw) |>
      mutate(coord_key = paste(chrom, start, end, strand, sep = ":"))
    message(sprintf("  %s ORFs with P-site purity", scales::comma(nrow(purity))))
  } else {
    message("  WARNING: psite_purity.tsv not found — purity filters will be NA")
    purity <- NULL
  }

  # ── Merge via orf_id and coord_key ──
  message("Merging data sources ...")
  merged <- meta |>
    mutate(coord_key = paste(chrom, start, end, strand, sep = ":")) |>
    left_join(conf, by = "orf_id") |>
    left_join(expr_agg, by = "orf_id")

  # Join GENCODE by coordinate (cross-ID format mapping)
  merged <- merged |>
    left_join(gencode |> select(coord_key, orf_biotype, gene_biotype,
                                 phaseI_biotype, pep, orf_length),
              by = "coord_key")

  message(sprintf("  Merged: %s ORFs", scales::comma(nrow(merged))))
  message(sprintf("  With GENCODE match: %s (%.1f%%)",
    scales::comma(sum(!is.na(merged$orf_biotype))),
    100 * mean(!is.na(merged$orf_biotype))))

  # Return as named list
  list(
    metadata    = meta,
    confidence  = conf,
    expression  = expr,
    gencode     = gencode,
    purity      = purity,
    merged      = merged
  )
}


# ── 2. Expression-based filtering (replaces P-site purity when bedgraph
#       scanning is too expensive) ─────────────────────────────────────

#' Apply per-sample expression-based quality filters
#'
#' Uses unified_orfs_expression_summary.tsv ({sample}_reads and {sample}_pN
#' columns) as proxies for P-site purity. No bedgraph rescanning needed.
#'
#' Mapping:
#'   p_site_GSE  → {sample}_reads   (total reads overlapping ORF)
#'   p_site_pct  → {sample}_pN      (P-site periodicity, [0-1] or >1)
#'   p_site_pos  → dropped          (no direct proxy)
apply_expression_filters <- function(data, cfg, expr_raw) {
  filt <- cfg$filtering

  read_cols <- grep("_reads$", names(expr_raw), value = TRUE)
  pN_cols   <- grep("_pN$", names(expr_raw), value = TRUE)

  # Per-sample pass/fail from expression columns (using reframe)
  n_samples_total <- expr_raw |>
    reframe(n = rowSums(pick(all_of(read_cols)) > 0, na.rm = TRUE)) |>
    pull(n)

  n_pass_reads <- expr_raw |>
    reframe(across(all_of(read_cols), ~ . > filt$p_site_gse_min)) |>
    as.matrix() |> rowSums(na.rm = TRUE)

  n_pass_pN <- expr_raw |>
    reframe(across(all_of(pN_cols), ~ . > filt$p_site_pct_min)) |>
    as.matrix() |> rowSums(na.rm = TRUE)

  # Triple: reads > min AND pN > min (align by shared sample names)
  sample_names <- intersect(sub("_reads$", "", read_cols),
                            sub("_pN$", "", pN_cols))
  triple_pass <- rep(0L, nrow(expr_raw))
  for (s in sample_names) {
    rcol <- paste0(s, "_reads")
    pcol <- paste0(s, "_pN")
    if (rcol %in% names(expr_raw) && pcol %in% names(expr_raw)) {
      triple_pass <- triple_pass +
        (expr_raw[[rcol]] > filt$p_site_gse_min &
         expr_raw[[pcol]] > filt$p_site_pct_min)
    }
  }
  n_pass_triple <- triple_pass

  data |>
    mutate(
      n_pass_gse  = n_pass_reads,
      n_pass_pct  = n_pass_pN,
      n_pass_pos  = NA_integer_,   # no proxy for position offset
      n_pass_all  = n_pass_triple,
      n_samples_purity = n_samples_total,
      pass_gse_any   = coalesce(n_pass_gse, 0) >= 1,
      pass_gse_2     = coalesce(n_pass_gse, 0) >= 2,
      pass_gse_3     = coalesce(n_pass_gse, 0) >= 3,
      pass_pct_any   = coalesce(n_pass_pct, 0) >= 1,
      pass_pct_2     = coalesce(n_pass_pct, 0) >= 2,
      pass_pct_3     = coalesce(n_pass_pct, 0) >= 3,
      pass_pos_any   = NA,
      pass_pos_2     = NA,
      pass_pos_3     = NA,
      pass_triple_any = coalesce(n_pass_all, 0) >= 1,
      pass_triple_2   = coalesce(n_pass_all, 0) >= 2,
      pass_triple_3   = coalesce(n_pass_all, 0) >= 3
    )
}


# ── 3. Cross-analysis helpers ──────────────────────────────────────────

#' Build a comparison matrix: rows = filter strategies, cols = metrics
build_filter_matrix <- function(df) {
  filter_cols <- names(df)[grep("^pass_", names(df))]
  filter_labels <- c(
    pass_gse_any   = "GSE>9 (any)",    pass_gse_2   = "GSE>9 (≥2)",
    pass_gse_3     = "GSE>9 (≥3)",
    pass_pct_any   = "pct>50% (any)",  pass_pct_2   = "pct>50% (≥2)",
    pass_pct_3     = "pct>50% (≥3)",
    pass_pos_any   = "pos>2 (any)",    pass_pos_2   = "pos>2 (≥2)",
    pass_pos_3     = "pos>2 (≥3)",
    pass_triple_any = "Triple (any)",  pass_triple_2 = "Triple (≥2)",
    pass_triple_3   = "Triple (≥3)"
  )

  map_dfr(setNames(filter_cols, filter_cols), function(fc) {
    sub <- df |> filter(.data[[fc]])
    tibble(
      filter_id  = fc,
      filter_label = filter_labels[fc],
      n_retained = nrow(sub),
      pct_total  = 100 * nrow(sub) / nrow(df),
      n_OCS_High = sum(sub$tier == "High", na.rm = TRUE),
      n_OCS_Medium = sum(sub$tier == "Medium", na.rm = TRUE),
      n_OCS_Low  = sum(sub$tier == "Low", na.rm = TRUE),
      n_nonCDS   = sum(!is.na(sub$orf_biotype) & sub$orf_biotype != "CDS",
                       na.rm = TRUE),
      n_CDS      = sum(sub$orf_biotype == "CDS", na.rm = TRUE),
      n_multi_tool = sum(sub$n_detecting >= 2, na.rm = TRUE),
      median_pN  = median(sub$pN, na.rm = TRUE)
    )
  })
}


# ── 4. ggRibo helpers (reusable across projects) ──────────────────────

#' Minimal Range_info R6 class for ggRibo compatibility
Range_info <- R6::R6Class("Range_info",
  public = list(
    exonsByTx  = NULL, txByGene = NULL, cdsByTx = NULL,
    fiveUTR = NULL, threeUTR = NULL, tx_to_gene = NULL,
    initialize = function(exonsByTx, txByGene, cdsByTx,
                          fiveUTR, threeUTR, tx_to_gene) {
      self$exonsByTx  <- exonsByTx
      self$txByGene   <- txByGene
      self$cdsByTx    <- cdsByTx
      self$fiveUTR    <- fiveUTR
      self$threeUTR   <- threeUTR
      self$tx_to_gene <- tx_to_gene
    }
  )
)

#' Import GTF annotation for ggRibo
gtf_import_custom <- function(annotation, format = "gtf",
                               dataSource = "", organism = "") {
  txdb <- suppressWarnings(
    txdbmaker::makeTxDbFromGFF(file = annotation, format = format,
                               dataSource = dataSource, organism = organism)
  )
  exonsByTx  <- GenomicFeatures::exonsBy(txdb, by = "tx", use.names = TRUE)
  txByGene   <- GenomicFeatures::transcriptsBy(txdb, by = "gene")
  cdsByTx    <- GenomicFeatures::cdsBy(txdb, by = "tx", use.names = TRUE)
  fiveUTR    <- GenomicFeatures::fiveUTRsByTranscript(txdb, use.names = TRUE)
  threeUTR   <- GenomicFeatures::threeUTRsByTranscript(txdb, use.names = TRUE)
  tx_to_gene <- AnnotationDbi::select(
    txdb, keys = AnnotationDbi::keys(txdb, keytype = "TXNAME"),
    columns = c("TXNAME", "GENEID"), keytype = "TXNAME"
  )
  colnames(tx_to_gene) <- c("tx_id", "gene_id")
  tx_to_gene <- tx_to_gene[order(tx_to_gene$tx_id), ]

  Txome_Range <- Range_info$new(
    exonsByTx = exonsByTx, txByGene = txByGene, cdsByTx = cdsByTx,
    fiveUTR = fiveUTR, threeUTR = threeUTR, tx_to_gene = tx_to_gene
  )
  assign("Txome_Range", Txome_Range, envir = .GlobalEnv)
}

#' Create minimal single-ORF GTF for ggRibo
create_orf_gtf <- function(orf, out_path) {
  gene_id <- sprintf("AMP_%s_%s_%d_%d",
                     orf$org, orf$chrom, orf$start, orf$end)
  tx_id   <- paste0(gene_id, "_T001")
  exon_id <- paste0(gene_id, "_E001")
  prot_id <- paste0(gene_id, "_P001")

  lines <- c(
    sprintf('%s\tAMP\tgene\t%d\t%d\t.\t%s\t.\tgene_id "%s"; gene_source "AMP"; gene_biotype "protein_coding";',
            orf$chrom, orf$start, orf$end, orf$strand, gene_id),
    sprintf('%s\tAMP\ttranscript\t%d\t%d\t.\t%s\t.\tgene_id "%s"; transcript_id "%s"; gene_source "AMP"; gene_biotype "protein_coding"; transcript_source "AMP"; transcript_biotype "protein_coding";',
            orf$chrom, orf$start, orf$end, orf$strand, gene_id, tx_id),
    sprintf('%s\tAMP\tCDS\t%d\t%d\t.\t%s\t0\tgene_id "%s"; transcript_id "%s"; exon_number "1"; gene_source "AMP"; gene_biotype "protein_coding"; transcript_source "AMP"; transcript_biotype "protein_coding"; protein_id "%s";',
            orf$chrom, orf$start, orf$end, orf$strand, gene_id, tx_id, prot_id),
    sprintf('%s\tAMP\texon\t%d\t%d\t.\t%s\t.\tgene_id "%s"; transcript_id "%s"; exon_number "1"; gene_source "AMP"; gene_biotype "protein_coding"; transcript_source "AMP"; transcript_biotype "protein_coding"; exon_id "%s";',
            orf$chrom, orf$start, orf$end, orf$strand, gene_id, tx_id, exon_id)
  )
  writeLines(lines, out_path)
  list(gene_id = gene_id, tx_id = tx_id)
}

#' Find P-site bedgraph files for a sample
find_bedgraph_files <- function(cfg, sample_name) {
  ribo_dir <- cfg$input$riboseqc_dir
  plus_file  <- file.path(ribo_dir, paste0(sample_name, "_P_sites_plus.bedgraph"))
  minus_file <- file.path(ribo_dir, paste0(sample_name, "_P_sites_minus.bedgraph"))
  if (file.exists(plus_file) && file.exists(minus_file)) {
    return(list(plus = plus_file, minus = minus_file))
  }
  NULL
}

#' Find top N samples by read count for an ORF from expression data
find_top_samples <- function(expr_data, orf_id, n = 3) {
  matched <- expr_data[expr_data$orf_id == orf_id, ]
  if (nrow(matched) == 0) return(NULL)
  read_cols <- grep("_reads$", names(matched), value = TRUE)
  reads <- sort(as.numeric(matched[1, read_cols]), decreasing = TRUE)
  names(reads) <- read_cols
  reads <- reads[reads > 0]
  if (length(reads) > n) reads <- reads[1:n]
  result <- list()
  for (i in seq_along(reads)) {
    result[[i]] <- list(
      sample = sub("_reads$", "", names(reads)[i]),
      reads  = reads[i]
    )
  }
  result
}
