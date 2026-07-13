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

  # Apply defaults for filtering if not set
  if (is.null(cfg$filtering$prelim_reads_min))
    cfg$filtering$prelim_reads_min <- 9
  if (is.null(cfg$filtering$prelim_pN_min))
    cfg$filtering$prelim_pN_min <- 0.5
  if (is.null(cfg$filtering$prelim_cross_sample))
    cfg$filtering$prelim_cross_sample <- 1
  if (is.null(cfg$filtering$p_site_gse_min))
    cfg$filtering$p_site_gse_min <- 9
  if (is.null(cfg$filtering$p_site_pct_min))
    cfg$filtering$p_site_pct_min <- 0.5
  if (is.null(cfg$filtering$p_site_pos_min))
    cfg$filtering$p_site_pos_min <- 2
  if (is.null(cfg$filtering$final_cross_sample))
    cfg$filtering$final_cross_sample <- 2

  cfg
}


# ── 1. Data loading + ID mapping ──────────────────────────────────────

#' Load all pipeline outputs and join by (chrom, start, end, strand)
#'
#' Returns a list with elements: metadata, confidence, expression,
#' purity, gencode, merged (the unified cross-walk).
#'
#' Unclassified ORFs (no GENCODE match) are assigned biotype "intergenic"
#' — they are inherently non-CDS since they don't overlap any annotated
#' transcript. The original GENCODE biotype (including NA) is preserved
#' as `gencode_biotype`.
load_pipeline_data <- function(cfg) {
  message("Loading metadata ...")
  meta <- .read_tsv(cfg$input$metadata) |>
    select(orf_id, chrom, strand, start, end, length_aa,
           tools, samples, total_reads, total_psites, pN,
           any_of("is_cds_overlap"), any_of("overlapping_genes"))

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
  gencode <- .read_tsv(gencode_path)
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

  # ── P-site purity (real bedgraph-based, not proxies) ──
  purity_file <- file.path(cfg$input$output_dir, "psite_purity.tsv")
  if (file.exists(purity_file)) {
    message("Loading P-site purity (exact bedgraph computation) ...")
    purity <- .read_tsv(purity_file) |>
      mutate(coord_key = paste(chrom, start, end, strand, sep = ":"))
    message(sprintf("  %s ORFs with P-site purity", scales::comma(nrow(purity))))
  } else {
    message("  NOTE: psite_purity.tsv not found — run compute_psite_purity.py first")
    purity <- NULL
  }

  # ── Merge via orf_id and coord_key ──
  message("Merging data sources ...")
  merged <- meta |>
    mutate(coord_key = paste(chrom, start, end, strand, sep = ":")) |>
    left_join(conf, by = "orf_id") |>
    left_join(expr_agg, by = "orf_id")

  # Join GENCODE by coordinate (cross-ID format mapping)
  # Preserve original GENCODE biotype as gencode_biotype
  merged <- merged |>
    left_join(gencode |> select(coord_key, orf_biotype, gene_biotype,
                                 phaseI_biotype, pep, orf_length),
              by = "coord_key") |>
    rename(gencode_biotype = orf_biotype)

  # Assign final biotype: keep GENCODE classification where available,
  # unclassified ORFs → "intergenic" (inherently non-CDS)
  merged <- merged |>
    mutate(
      orf_biotype = if_else(
        !is.na(gencode_biotype), gencode_biotype, "intergenic"
      ),
      is_classified = !is.na(gencode_biotype)
    )

  # Ensure orf_length columns don't collide
  if ("orf_length.x" %in% names(merged)) {
    merged <- merged |>
      mutate(orf_length = coalesce(orf_length.y, orf_length.x)) |>
      select(-orf_length.x, -orf_length.y)
  }

  n_classified <- sum(merged$is_classified)
  n_intergenic <- sum(merged$orf_biotype == "intergenic")
  n_cds        <- sum(merged$orf_biotype == "CDS", na.rm = TRUE)
  n_noncds     <- sum(merged$orf_biotype != "CDS", na.rm = TRUE)

  message(sprintf("  Merged: %s ORFs", scales::comma(nrow(merged))))
  message(sprintf("  Classified (GENCODE): %s (%.1f%%)",
    scales::comma(n_classified), 100 * n_classified / nrow(merged)))
  message(sprintf("  Intergenic (unclassified): %s (%.1f%%)",
    scales::comma(n_intergenic), 100 * n_intergenic / nrow(merged)))
  message(sprintf("  CDS: %s  |  Non-CDS: %s (incl. intergenic)",
    scales::comma(n_cds), scales::comma(n_noncds)))

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


# ── 2. Stage 1: Preliminary expression-based filtering (cheap) ─────────
#       Uses expression-summary columns ({sample}_reads and {sample}_pN).
#       These are proxies — NOT real P-site purity. The purpose of this
#       stage is to reduce the ORF set before running the expensive
#       bedgraph-based P-site computation.
#
#       Mapping:
#         {sample}_reads  → signal presence (validated r=0.99 vs p_site_GSE)
#         {sample}_pN     → frame periodicity quality (NOT P-site fraction)
#                            pN measures 3-nt periodicity, r=-0.20 vs p_site_pct
#         position offset → dropped (no proxy; computed in Stage 2)

#' Apply preliminary expression-based quality filters
#'
#' Samples come from heterogeneous datasets, so per-sample pass/fail
#' is the default (≥1 sample). The goal is to eliminate ORFs with
#' insufficient signal before running the expensive bedgraph scan.
apply_prelim_filters <- function(data, cfg, expr_raw) {
  filt <- cfg$filtering

  read_cols <- grep("_reads$", names(expr_raw), value = TRUE)
  pN_cols   <- grep("_pN$", names(expr_raw), value = TRUE)

  # Per-sample pass/fail from expression columns
  n_samples_total <- expr_raw |>
    reframe(n = rowSums(pick(all_of(read_cols)) > 0, na.rm = TRUE)) |>
    pull(n)

  n_pass_reads <- expr_raw |>
    reframe(across(all_of(read_cols), ~ . > filt$prelim_reads_min)) |>
    as.matrix() |> rowSums(na.rm = TRUE)

  n_pass_pN <- expr_raw |>
    reframe(across(all_of(pN_cols), ~ . > filt$prelim_pN_min)) |>
    as.matrix() |> rowSums(na.rm = TRUE)

  # Combined: reads > min AND pN > min (align by shared sample names)
  sample_names <- intersect(sub("_reads$", "", read_cols),
                            sub("_pN$", "", pN_cols))
  combined_pass <- rep(0L, nrow(expr_raw))
  for (s in sample_names) {
    rcol <- paste0(s, "_reads")
    pcol <- paste0(s, "_pN")
    if (rcol %in% names(expr_raw) && pcol %in% names(expr_raw)) {
      combined_pass <- combined_pass +
        (expr_raw[[rcol]] > filt$prelim_reads_min &
         expr_raw[[pcol]] > filt$prelim_pN_min)
    }
  }

  cs <- filt$prelim_cross_sample  # default 1 (any sample)

  data |>
    mutate(
      n_pass_reads  = n_pass_reads,
      n_pass_pN     = n_pass_pN,
      n_pass_combined = combined_pass,
      n_samples_expr = n_samples_total,
      prelim_pass_reads   = coalesce(n_pass_reads, 0) >= cs,
      prelim_pass_pN      = coalesce(n_pass_pN, 0) >= cs,
      prelim_pass_combined = coalesce(combined_pass, 0) >= cs,
      prelim_pass_reads_any = coalesce(n_pass_reads, 0) >= 1,
      prelim_pass_reads_2   = coalesce(n_pass_reads, 0) >= 2,
      prelim_pass_reads_3   = coalesce(n_pass_reads, 0) >= 3,
      prelim_pass_combined_any = coalesce(combined_pass, 0) >= 1,
      prelim_pass_combined_2   = coalesce(combined_pass, 0) >= 2,
      prelim_pass_combined_3   = coalesce(combined_pass, 0) >= 3
    )
}


# ── 3. Stage 2: Real P-site purity filtering (bedgraph-based, exact) ──
#       Uses psite_purity.tsv generated by compute_psite_purity.py.
#       This script intersects ORF coordinates with actual RiboseQC
#       bedgraph files — no proxies, no approximations.

#' Apply real P-site purity filters from bedgraph computation
#'
#' Reads psite_purity.tsv (produced by compute_psite_purity.py) and
#' applies per-sample thresholds for p_site_GSE, p_site_pct, p_site_pos.
#'
#' The P-site purity data is computed by backtracking through RiboseQC
#' _P_sites_{plus,minus}.bedgraph and _coverage_{plus,minus}.bedgraph
#' files — exact per-ORF per-sample values, not proxies.
apply_psite_filters <- function(data, cfg, purity_tsv_path) {
  if (!file.exists(purity_tsv_path)) {
    warning("P-site purity file not found: ", purity_tsv_path,
            " — run compute_psite_purity.py first")
    return(data)
  }

  message("Loading P-site purity data ...")
  purity <- .read_tsv(purity_tsv_path)
  message(sprintf("  %s ORFs with P-site purity", scales::comma(nrow(purity))))

  filt <- cfg$filtering

  # Discover sample names from column suffixes
  psite_cols <- grep("_p_site_GSE$", names(purity), value = TRUE)
  pct_cols   <- grep("_p_site_pct$", names(purity), value = TRUE)
  pos_cols   <- grep("_p_site_pos$", names(purity), value = TRUE)

  # Use the orf_id column for joining — psite_purity.tsv uses the same orf_id
  # Per-sample pass/fail
  n_pass_gse <- purity |>
    reframe(across(all_of(psite_cols), ~ . > filt$p_site_gse_min)) |>
    as.matrix() |> rowSums(na.rm = TRUE)

  n_pass_pct <- purity |>
    reframe(across(all_of(pct_cols), ~ . > filt$p_site_pct_min)) |>
    as.matrix() |> rowSums(na.rm = TRUE)

  n_pass_pos <- purity |>
    reframe(across(all_of(pos_cols), ~ . > filt$p_site_pos_min)) |>
    as.matrix() |> rowSums(na.rm = TRUE)

  # Per-ORF aggregates
  purity_agg <- purity |>
    mutate(
      n_pass_gse = n_pass_gse,
      n_pass_pct = n_pass_pct,
      n_pass_pos = n_pass_pos,
      # Total P-site reads across all samples
      total_psites_real = rowSums(pick(all_of(psite_cols)), na.rm = TRUE),
      # Mean P-site fraction across samples with signal
      mean_p_site_pct = rowMeans(
        pick(all_of(pct_cols)) |> mutate(across(everything(), ~ if_else(. > 0, ., NA_real_))),
        na.rm = TRUE
      )
    )

  # Triple: all three criteria simultaneously
  sample_names <- intersect(
    intersect(sub("_p_site_GSE$", "", psite_cols),
              sub("_p_site_pct$", "", pct_cols)),
    sub("_p_site_pos$", "", pos_cols)
  )
  triple_pass <- rep(0L, nrow(purity))
  for (s in sample_names) {
    gse_col <- paste0(s, "_p_site_GSE")
    pct_col <- paste0(s, "_p_site_pct")
    pos_col <- paste0(s, "_p_site_pos")
    if (all(c(gse_col, pct_col, pos_col) %in% names(purity))) {
      triple_pass <- triple_pass +
        (purity[[gse_col]] > filt$p_site_gse_min &
         purity[[pct_col]] > filt$p_site_pct_min &
         purity[[pos_col]] > filt$p_site_pos_min)
    }
  }
  purity_agg$n_pass_triple <- triple_pass

  cs <- filt$final_cross_sample  # default 2

  purity_agg <- purity_agg |>
    mutate(
      psite_pass_gse   = coalesce(n_pass_gse, 0) >= cs,
      psite_pass_pct   = coalesce(n_pass_pct, 0) >= cs,
      psite_pass_pos   = coalesce(n_pass_pos, 0) >= cs,
      psite_pass_triple = coalesce(n_pass_triple, 0) >= cs,
      psite_pass_gse_any   = coalesce(n_pass_gse, 0) >= 1,
      psite_pass_gse_2     = coalesce(n_pass_gse, 0) >= 2,
      psite_pass_gse_3     = coalesce(n_pass_gse, 0) >= 3,
      psite_pass_triple_any = coalesce(n_pass_triple, 0) >= 1,
      psite_pass_triple_2   = coalesce(n_pass_triple, 0) >= 2,
      psite_pass_triple_3   = coalesce(n_pass_triple, 0) >= 3
    ) |>
    select(orf_id, n_pass_gse, n_pass_pct, n_pass_pos, n_pass_triple,
           total_psites_real, mean_p_site_pct,
           starts_with("psite_pass_"))

  # Merge back into data
  data |>
    left_join(purity_agg, by = "orf_id")
}


# ── 4. Cross-analysis helpers ──────────────────────────────────────────

#' Build a comparison matrix: rows = filter strategies, cols = metrics
#'
#' Handles both prelim_pass_* and psite_pass_* columns.
build_filter_matrix <- function(df, filter_pattern = "^(prelim_pass_|psite_pass_)") {
  filter_cols <- names(df)[grep(filter_pattern, names(df))]
  if (length(filter_cols) == 0) {
    filter_cols <- names(df)[grep("^pass_", names(df))]  # legacy fallback
  }

  filter_labels <- c(
    # Stage 1: Preliminary
    prelim_pass_reads    = "Prelim: reads (≥1)",
    prelim_pass_reads_any = "Prelim: reads (any)",
    prelim_pass_reads_2  = "Prelim: reads (≥2)",
    prelim_pass_reads_3  = "Prelim: reads (≥3)",
    prelim_pass_pN       = "Prelim: pN (≥1)",
    prelim_pass_combined = "Prelim: reads+pN (≥1)",
    prelim_pass_combined_any = "Prelim: reads+pN (any)",
    prelim_pass_combined_2   = "Prelim: reads+pN (≥2)",
    prelim_pass_combined_3   = "Prelim: reads+pN (≥3)",
    # Stage 2: Real P-site
    psite_pass_gse       = "P-site: GSE (≥2)",
    psite_pass_gse_any   = "P-site: GSE (any)",
    psite_pass_gse_2     = "P-site: GSE (≥2)",
    psite_pass_gse_3     = "P-site: GSE (≥3)",
    psite_pass_pct       = "P-site: pct (≥2)",
    psite_pass_pos       = "P-site: pos (≥2)",
    psite_pass_triple    = "P-site: Triple (≥2)",
    psite_pass_triple_any = "P-site: Triple (any)",
    psite_pass_triple_2   = "P-site: Triple (≥2)",
    psite_pass_triple_3   = "P-site: Triple (≥3)",
    # Legacy (backward-compatible)
    pass_gse_any   = "GSE>9 (any)",    pass_gse_2   = "GSE>9 (≥2)",
    pass_gse_3     = "GSE>9 (≥3)",
    pass_pct_any   = "pct>50% (any)",  pass_pct_2   = "pct>50% (≥2)",
    pass_pct_3     = "pct>50% (≥3)",
    pass_pos_any   = "pos>2 (any)",    pass_pos_2   = "pos>2 (≥2)",
    pass_pos_3     = "pos>2 (≥3)",
    pass_triple_any = "Triple (any)",  pass_triple_2 = "Triple (≥2)",
    pass_triple_3   = "Triple (≥3)"
  )

  results <- lapply(setNames(filter_cols, filter_cols), function(fc) {
    sub <- df |> filter(.data[[fc]])
    data.frame(
      filter_id  = fc,
      filter_label = ifelse(fc %in% names(filter_labels),
                            filter_labels[fc], fc),
      n_retained = nrow(sub),
      pct_total  = 100 * nrow(sub) / nrow(df),
      n_OCS_High = sum(sub$tier == "High", na.rm = TRUE),
      n_OCS_Medium = sum(sub$tier == "Medium", na.rm = TRUE),
      n_OCS_Low  = sum(sub$tier == "Low", na.rm = TRUE),
      n_nonCDS   = sum(sub$orf_biotype != "CDS", na.rm = TRUE),
      n_CDS      = sum(sub$orf_biotype == "CDS", na.rm = TRUE),
      n_intergenic = sum(sub$orf_biotype == "intergenic", na.rm = TRUE),
      n_multi_tool = sum(sub$n_detecting >= 2, na.rm = TRUE),
      median_reads = median(sub$total_reads_expr, na.rm = TRUE),
      median_pN  = median(sub$pN, na.rm = TRUE)
    )
  })
  do.call(rbind, results) |> as.data.frame(stringsAsFactors = FALSE)
}


# ── 5. ggRibo helpers (reusable across projects) ──────────────────────

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

#' Find top N samples by a chosen metric for an ORF
#'
#' @param expr_data Expression summary table (has {sample}_reads and {sample}_pN)
#' @param psite_data P-site purity table (has {sample}_p_site_GSE, {sample}_p_site_pct, {sample}_p_site_pos)
#' @param orf_id ORF identifier
#' @param n Number of top samples to return
#' @param sort_by Metric to sort by: "reads", "p_site_GSE", "p_site_pct", "pN"
#' @return List of lists, each with sample, reads, and metric-specific value
find_top_samples <- function(expr_data, orf_id, n = 10,
                             psite_data = NULL,
                             sort_by = c("reads", "p_site_GSE", "p_site_pct", "pN")) {
  sort_by <- match.arg(sort_by)

  matched <- expr_data[expr_data$orf_id == orf_id, ]
  if (nrow(matched) == 0) return(NULL)

  read_cols <- grep("_reads$", names(matched), value = TRUE)

  # Determine sort values based on metric
  if (sort_by == "reads") {
    sort_vals <- as.numeric(matched[1, read_cols])
    names(sort_vals) <- sub("_reads$", "", read_cols)
  } else if (sort_by == "pN") {
    pN_cols <- grep("_pN$", names(matched), value = TRUE)
    sort_vals <- as.numeric(matched[1, pN_cols])
    names(sort_vals) <- sub("_pN$", "", pN_cols)
  } else if (sort_by %in% c("p_site_GSE", "p_site_pct") && !is.null(psite_data)) {
    # Look up from real P-site data
    psm <- psite_data[psite_data$orf_id == orf_id, ]
    if (nrow(psm) == 0) {
      # Fall back to reads
      sort_vals <- as.numeric(matched[1, read_cols])
      names(sort_vals) <- sub("_reads$", "", read_cols)
    } else {
      metric_cols <- grep(paste0("_", sort_by, "$"), names(psm), value = TRUE)
      sort_vals <- as.numeric(psm[1, metric_cols])
      names(sort_vals) <- sub(paste0("_", sort_by, "$"), "", metric_cols)
    }
  } else {
    # Default: reads
    sort_vals <- as.numeric(matched[1, read_cols])
    names(sort_vals) <- sub("_reads$", "", read_cols)
  }

  # Get corresponding reads for each sample
  reads_named <- as.numeric(matched[1, read_cols])
  names(reads_named) <- sub("_reads$", "", read_cols)

  # Sort descending, keep top N with positive sort values
  sort_vals <- sort(sort_vals, decreasing = TRUE)
  sort_vals <- sort_vals[sort_vals > 0]
  if (length(sort_vals) > n) sort_vals <- sort_vals[1:n]

  result <- list()
  for (i in seq_along(sort_vals)) {
    sn <- names(sort_vals)[i]
    result[[i]] <- list(
      sample      = sn,
      reads       = unname(reads_named[sn]),
      sort_value  = unname(sort_vals[i]),
      sort_metric = sort_by
    )
  }
  result
}
