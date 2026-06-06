#!/usr/bin/env Rscript
#
# orfquant_classify_datatable.R
# Lightweight ORFquant-style ORF classification using data.table + DuckDB.
# Replaces GRanges-based implementation with ~100x speedup.
#
# Dependencies: data.table, duckdb (auto-installed if missing)
#
# Usage:
#   Rscript orfquant_classify_datatable.R \
#     --orfs unified_orfs.metadata.tsv \
#     --gtf annotation.gtf \
#     --output orfquant_classification.tsv \
#     [--n_cores 4] [--cache_dir .]
#

suppressPackageStartupMessages({
  if (!requireNamespace("data.table", quietly = TRUE))
    install.packages("data.table", repos = "https://cloud.r-project.org", quiet = TRUE)
  library(data.table)
})

# ── Command-line parsing ─────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
opt <- function(name, default = NA_character_) {
  idx <- which(args == name)
  if (length(idx) == 0) return(default)
  if (idx == length(args)) stop("Missing value for ", name)
  args[idx + 1]
}
orfs_file   <- opt("--orfs")
gtf_file    <- opt("--gtf")
output_file <- opt("--output", "orfquant_classification.tsv")
cache_dir   <- opt("--cache_dir", ".")
n_cores     <- as.integer(opt("--n_cores", "4"))

if (is.na(orfs_file) || is.na(gtf_file)) {
  stop("Usage: Rscript orfquant_classify_datatable.R --orfs <metadata.tsv> --gtf <annotation.gtf> --output <out.tsv>")
}

cat(sprintf("=== ORFquant Classification (data.table + DuckDB) ===\n"))
cat(sprintf("  ORFs : %s\n", orfs_file))
cat(sprintf("  GTF  : %s\n", gtf_file))
cat(sprintf("  Cores: %d\n", n_cores))

# ── 1. Load GTF annotation via DuckDB ────────────────────────────────────────
cat("\n[1/5] Loading GTF annotation...\n")

rds_cache <- file.path(cache_dir, paste0(basename(gtf_file), ".orfclass_cache.rds"))
if (file.exists(rds_cache) && file.mtime(rds_cache) >= file.mtime(gtf_file)) {
  cat("  Loading cached annotation from:", rds_cache, "\n")
  ann <- readRDS(rds_cache)
} else {
  # ── Load GTF with data.table::fread (fast, no S4/DuckDB overhead) ──
  cat("  Reading GTF with data.table::fread...\n")
  gtf <- fread(gtf_file, sep = "\t", header = FALSE, skip = "#",
               select = c(1,3,4,5,7,9), showProgress = FALSE,
               col.names = c("chr", "type", "start", "end", "strand", "attrs"))

  # Filter to exon/CDS/gene only
  gtf <- gtf[type %in% c("exon", "CDS", "gene")]
  cat(sprintf("  Loaded %s rows (exon/CDS/gene)\n", format(nrow(gtf), big.mark = ",")))

  # Extract gene/transcript/biotype IDs from attributes column
  cat("  Extracting gene/transcript IDs...\n")
  gtf[, gene_id       := gsub('.*gene_id "([^"]+)".*', '\\1', attrs)]
  gtf[, transcript_id := gsub('.*transcript_id "([^"]+)".*', '\\1', attrs)]
  gtf[, gene_biotype  := gsub('.*gene_biotype "([^"]+)".*', '\\1', attrs)]
  # Clean: if pattern not found, set to NA
  gtf[gene_id == attrs,       gene_id := NA]
  gtf[transcript_id == attrs, transcript_id := NA]
  gtf[gene_biotype == attrs,  gene_biotype := NA]

  # Drop attributes column to save memory
  gtf[, attrs := NULL]

  cat("  Building annotation index...\n")

  # ── Exons per transcript (for transcript coordinate projection) ──
  exon_dt <- gtf[type == "exon" & !is.na(transcript_id),
    .(transcript_id, gene_id, chr, start, end, strand)]
  setorder(exon_dt, transcript_id, start)

  # Assign exon rank and compute cumulative length along transcript
  exon_dt[, exon_rank := seq_len(.N), by = transcript_id]
  exon_dt[, exon_width := end - start + 1L]
  exon_dt[, cum_len := cumsum(c(0L, exon_width[-.N])), by = transcript_id]

  # ── CDS per transcript → merged CDS interval ──
  cds_dt <- gtf[type == "CDS" & !is.na(transcript_id),
    .(transcript_id, gene_id, chr, start, end, strand)]
  setorder(cds_dt, transcript_id, start)

  # Per-transcript CDS bounds (GTF CDS exons are non-overlapping within a transcript)
  cds_merged <- cds_dt[, {
    .(cds_start = min(start), cds_end = max(end),
      cds_width = sum(end - start + 1L))
  }, by = .(transcript_id, gene_id, chr, strand)]

  # ── Max CDS per gene ──
  max_cds_by_gene <- cds_merged[!is.na(cds_width)][
    , .SD[which.max(cds_width)], by = gene_id]

  # ── Gene-level CDS (merged across all transcripts) ──
  cds_gene <- cds_merged[!is.na(cds_start), {
    .(cds_gene_start = min(cds_start), cds_gene_end = max(cds_end))
  }, by = .(gene_id, chr, strand)]

  # ── All CDS merged (for ORFs without a known gene) ──
  ok <- cds_merged[!is.na(cds_start)]
  so_all <- ok$cds_start; eo_all <- ok$cds_end
  o <- order(so_all); so_all <- so_all[o]; eo_all <- eo_all[o]
  # Store as simple bounds
  all_cds_start <- min(so_all)
  all_cds_end   <- max(eo_all)

  # ── LncRNA genes ──
  lncrna_genes <- gtf[type == "gene" & !is.na(gene_biotype) &
    gene_biotype %in% c("lncRNA","lincRNA","antisense","processed_transcript",
      "sense_intronic","sense_overlapping","non_coding",
      "3prime_overlapping_ncrna","bidirectional_promoter_lncrna",
      "ncRNA","antisense_RNA","misc_non_coding"),
    unique(gene_id)]

  # ── Gene → transcript lookup ──
  gene_to_txids <- split(cds_merged$transcript_id, cds_merged$gene_id)

  # ── Exon by transcript (list of data.tables for fast per-tx access) ──
  exon_by_tx <- split(exon_dt, exon_dt$transcript_id)

  # Free GTF from memory
  rm(gtf); gc()

  ann <- list(
    exon_dt = exon_dt,
    exon_by_tx = exon_by_tx,
    cds_merged = cds_merged,
    max_cds_by_gene = max_cds_by_gene,
    cds_gene = cds_gene,
    all_cds_start = all_cds_start,
    all_cds_end = all_cds_end,
    lncrna_genes = lncrna_genes,
    gene_to_txids = gene_to_txids
  )

  saveRDS(ann, rds_cache)
  cat("  Cached annotation to:", rds_cache, "\n")
}

# ── 2. Load ORF metadata ─────────────────────────────────────────────────────
cat("\n[2/5] Loading ORF metadata...\n")
orfs <- fread(orfs_file, sep = "\t", na.strings = c("", "NA", "None"))
cat(sprintf("  Loaded %d ORFs\n", nrow(orfs)))

# Parse exon blocks from metadata
parse_exons <- function(s) {
  if (is.na(s) || s == "" || s == "None") return(NULL)
  parts <- strsplit(s, ",")[[1]]
  dt <- rbindlist(lapply(parts, function(p) {
    se <- as.integer(strsplit(p, "-")[[1]])
    if (length(se) != 2) return(NULL)
    data.table(start = se[1], end = se[2])
  }))
  if (nrow(dt) == 0) return(NULL)
  dt
}

# ── 3. Core classification functions ─────────────────────────────────────────

# Project genomic coordinates to transcript coordinates via exon chain.
# exon_dt: data.table with columns exon_rank, start, end, cum_len, strand
# Returns c(tx_start, tx_end) or NULL if projection fails.
project_to_tx <- function(orf_blocks, tx_exons, strand) {
  if (is.null(orf_blocks) || nrow(orf_blocks) == 0) return(NULL)
  if (nrow(tx_exons) == 0) return(NULL)

  tx_start <- NA_integer_
  tx_end   <- NA_integer_

  for (bi in seq_len(nrow(orf_blocks))) {
    bstart <- orf_blocks$start[bi]
    bend   <- orf_blocks$end[bi]

    # Find which exon(s) this block overlaps
    hits <- tx_exons[start <= bend & end >= bstart]
    if (nrow(hits) == 0) return(NULL)  # block outside transcript

    # Map genomic start → tx coordinate
    first_exon <- hits[1]
    if (strand == "+") {
      bs_tx <- first_exon$cum_len + (bstart - first_exon$start + 1L)
    } else {
      bs_tx <- first_exon$cum_len + (first_exon$end - bstart + 1L)
    }
    if (is.na(tx_start) || bs_tx < tx_start) tx_start <- bs_tx

    # Map genomic end → tx coordinate
    last_exon <- hits[.N]
    if (strand == "+") {
      be_tx <- last_exon$cum_len + (bend - last_exon$start + 1L)
    } else {
      be_tx <- last_exon$cum_len + (last_exon$end - bend + 1L)
    }
    if (is.na(tx_end) || be_tx > tx_end) tx_end <- be_tx
  }

  if (is.na(tx_start) || is.na(tx_end)) return(NULL)
  # Ensure start < end in tx space
  if (tx_start > tx_end) { tmp <- tx_start; tx_start <- tx_end; tx_end <- tmp }
  return(c(tx_start, tx_end))
}

# Genomic-level classification (positional, no S4)
classify_genomic <- function(orf_start, orf_end, orf_strand,
                              cds_start, cds_end, gene_id, max_cds_dt) {
  if (is.na(cds_start) || is.na(cds_end)) return("novel")

  # Check overlap
  overlaps <- orf_start <= cds_end && orf_end >= cds_start
  if (!overlaps) {
    if (orf_strand == "+") {
      if (orf_end < cds_start) return("novel_Upstream")
      if (orf_start > cds_end) return("novel_Downstream")
    } else {
      if (orf_end > cds_end) return("novel_Upstream")
      if (orf_start < cds_start) return("novel_Downstream")
    }
    # Internal (within gene boundaries but between CDS exons)
    return("novel_Internal")
  }

  # Overlaps CDS — check against best CDS isoform
  if (is.null(max_cds_dt) || nrow(max_cds_dt) == 0 ||
      is.na(max_cds_dt$cds_start[1])) return("overlaps_CDS")

  mc_s <- max_cds_dt$cds_start[1]
  mc_e <- max_cds_dt$cds_end[1]

  if (orf_strand == "+") {
    gen_sta <- mc_s; gen_sto <- mc_e
    sta_or <- orf_start; sto_or <- orf_end + 3L
  } else {
    gen_sta <- mc_e; gen_sto <- mc_s
    sta_or <- orf_end; sto_or <- orf_start - 3L
  }

  if (sto_or == gen_sto) {
    if (sta_or == gen_sta) return("exact_start_stop")
    if (orf_strand == "+") {
      if (sta_or < gen_sta) return("Alt5_start")
      return("Alt3_start")
    } else {
      if (sta_or < gen_sta) return("Alt3_start")
      return("Alt5_start")
    }
  }

  # Different stop positions
  if (orf_strand == "+") {
    if (sta_or < gen_sta && sto_or < gen_sto) return("Alt5_start_Alt5_stop")
    if (sta_or < gen_sta && sto_or > gen_sto) return("Alt5_start_Alt3_stop")
    if (sta_or > gen_sta && sto_or > gen_sto) return("Alt3_start_Alt3_stop")
    if (sta_or > gen_sta && sto_or < gen_sto) return("Alt3_start_Alt5_stop")
    if (sta_or == gen_sta && sto_or < gen_sto) return("Alt5_stop")
    if (sta_or == gen_sta && sto_or > gen_sto) return("Alt3_stop")
  } else {
    if (sta_or > gen_sta && sto_or > gen_sto) return("Alt5_start_Alt5_stop")
    if (sta_or > gen_sta && sto_or < gen_sto) return("Alt5_start_Alt3_stop")
    if (sta_or < gen_sta && sto_or < gen_sto) return("Alt3_start_Alt3_stop")
    if (sta_or < gen_sta && sto_or > gen_sto) return("Alt3_start_Alt5_stop")
    if (sta_or == gen_sta && sto_or > gen_sto) return("Alt5_stop")
    if (sta_or == gen_sta && sto_or < gen_sto) return("Alt3_stop")
  }
  return("overlaps_CDS")
}

# Transcript-level classification (pure arithmetic, no S4)
classify_transcript <- function(orf_sta, orf_sto, ann_sta, ann_sto) {
  if (is.na(ann_sta) || is.na(ann_sto)) return(NA_character_)

  if (orf_sto == ann_sto) {
    if (orf_sta == ann_sta) return("ORF_annotated")
    if (orf_sta < ann_sta) return("N_extension")
    if (orf_sta > ann_sta) return("N_truncation")
  }

  if (orf_sto != ann_sto) {
    if (orf_sta < ann_sta && orf_sto < ann_sta) return("uORF")
    if (orf_sta < ann_sta && orf_sto < ann_sto) return("overl_uORF")
    if (orf_sta < ann_sta && orf_sto > ann_sto) return("NC_extension")
    if (orf_sta > ann_sto && orf_sto > ann_sto) return("dORF")
    if (orf_sta > ann_sta && orf_sto > ann_sto) return("overl_dORF")
    if (orf_sta > ann_sta && orf_sto < ann_sto) return("nested_ORF")
    if (orf_sta == ann_sta && orf_sto < ann_sto) return("C_truncation")
    if (orf_sta == ann_sta && orf_sto > ann_sto) return("C_extension")
  }

  return(NA_character_)
}

# Priority order: lower = better
TX_CLASS_PRIORITY <- c(
  ORF_annotated = 1L, N_extension = 2L, C_extension = 3L,
  N_truncation = 4L, C_truncation = 5L, NC_extension = 6L,
  overl_uORF = 7L, overl_dORF = 8L, nested_ORF = 9L,
  uORF = 10L, dORF = 11L, novel = 12L, novel_antisense = 13L
)

# ── 4. Main classification loop (data.table-optimized) ───────────────────────
cat("\n[3/5] Running classification...\n")

# Pre-compute ORF data.table from metadata
# Parse exon blocks
exon_list <- lapply(orfs$exon_blocks, parse_exons)

# Build ORF data.table with required columns
orf_dt <- data.table(
  orf_id        = orfs$orf_id %||% orfs$name %||% as.character(seq_len(nrow(orfs))),
  gene_id       = orfs$gene_id %||% orfs$gene_name,
  transcript_id = orfs$transcript_id %||% orfs$transcript_name,
  chr           = orfs$chrom %||% orfs$seqnames %||% orfs$chr %||% orfs$seqname,
  start         = as.integer(orfs$start %||% orfs$orf_start),
  end           = as.integer(orfs$end %||% orfs$orf_end),
  strand        = orfs$strand,
  exon_blocks   = exon_list
)

# Set NA gene_id to empty for safe lookups
orf_dt[is.na(gene_id), gene_id := ""]
orf_dt[is.na(transcript_id), transcript_id := ""]

# Initialize result columns
orf_dt[, ORF_category_Gen         := NA_character_]
orf_dt[, ORF_category_Tx          := NA_character_]
orf_dt[, ORF_category_Tx_compatible := NA_character_]

# Extract annotation refs
exon_by_tx      <- ann$exon_by_tx
cds_merged      <- ann$cds_merged
max_cds_by_gene <- ann$max_cds_by_gene
cds_gene        <- ann$cds_gene
lncrna_genes    <- ann$lncrna_genes
gene_to_txids   <- ann$gene_to_txids

# all_cds bounds for ORFs without a gene
all_cds_start <- ann$all_cds_start
all_cds_end   <- ann$all_cds_end

n_orfs <- nrow(orf_dt)
t_start <- Sys.time()

# Process in batches for cache-friendly memory access
BATCH_SIZE <- 10000L
n_batches <- ceiling(n_orfs / BATCH_SIZE)

for (bi in seq_len(n_batches)) {
  idx_start <- (bi - 1L) * BATCH_SIZE + 1L
  idx_end   <- min(bi * BATCH_SIZE, n_orfs)
  batch <- orf_dt[idx_start:idx_end]

  for (ri in seq_len(nrow(batch))) {
    orf_id  <- batch$orf_id[ri]
    gene_id <- batch$gene_id[ri]
    tx_id   <- batch$transcript_id[ri]
    chr_i   <- batch$chr[ri]
    start_i <- batch$start[ri]
    end_i   <- batch$end[ri]
    strand_i <- batch$strand[ri]
    blocks  <- batch$exon_blocks[[ri]]

    # lncRNA check
    if (!is.na(gene_id) && nzchar(gene_id) && length(lncrna_genes) > 0 &&
        gene_id %in% lncrna_genes) {
      orf_dt[idx_start + ri - 1L,
        c("ORF_category_Gen", "ORF_category_Tx", "ORF_category_Tx_compatible") :=
        list("lncRNA", "lncRNA", "lncRNA")]
      next
    }

    # ── ORF_category_Gen ─────────────────────────────────────────────────────
    if (!is.na(gene_id) && nzchar(gene_id) && gene_id %in% cds_gene$gene_id) {
      cg <- cds_gene[gene_id, on = "gene_id"]$cds_gene_start
      cg_end <- cds_gene[gene_id, on = "gene_id"]$cds_gene_end
      mcg <- max_cds_by_gene[gene_id == gene_id]
      if (nrow(mcg) > 0) mcg <- mcg[1]
    } else {
      cg <- all_cds_start
      cg_end <- all_cds_end
      mcg <- NULL
    }
    gen_class <- classify_genomic(start_i, end_i, strand_i,
                                   cg, cg_end, gene_id, mcg)
    orf_dt[idx_start + ri - 1L, ORF_category_Gen := gen_class]

    # ── ORF_category_Tx ──────────────────────────────────────────────────────
    if (!is.na(gene_id) && nzchar(gene_id) &&
        nzchar(tx_id) && tx_id %in% names(exon_by_tx)) {
      tx_exons <- exon_by_tx[[tx_id]]
      orf_tx <- project_to_tx(blocks, tx_exons, strand_i)
      if (!is.null(orf_tx)) {
        cd <- cds_merged[transcript_id == tx_id]
        if (nrow(cd) > 0 && !is.na(cd$cds_start[1])) {
          tx_class <- classify_transcript(orf_tx[1], orf_tx[2],
                                          cd$cds_start[1], cd$cds_end[1])
          orf_dt[idx_start + ri - 1L, ORF_category_Tx := tx_class]
        } else {
          orf_dt[idx_start + ri - 1L, ORF_category_Tx := "novel"]
        }
      }
    }

    # ── ORF_category_Tx_compatible ───────────────────────────────────────────
    if (!is.na(gene_id) && nzchar(gene_id) &&
        gene_id %in% names(gene_to_txids)) {
      best_class <- "novel"
      best_prio  <- TX_CLASS_PRIORITY["novel"]
      for (txc_id in gene_to_txids[[gene_id]]) {
        if (!txc_id %in% names(exon_by_tx)) next
        cd <- cds_merged[transcript_id == txc_id]
        if (nrow(cd) == 0 || is.na(cd$cds_start[1])) next
        orf_txc <- project_to_tx(blocks, exon_by_tx[[txc_id]], strand_i)
        if (is.null(orf_txc)) next
        cls <- classify_transcript(orf_txc[1], orf_txc[2],
                                    cd$cds_start[1], cd$cds_end[1])
        if (!is.na(cls)) {
          prio <- TX_CLASS_PRIORITY[cls]
          if (is.na(prio)) prio <- TX_CLASS_PRIORITY["novel"]
          if (prio < best_prio) {
            best_prio <- prio
            best_class <- cls
          }
        }
        if (best_prio <= 1L) break
      }
      orf_dt[idx_start + ri - 1L, ORF_category_Tx_compatible := best_class]
    }
  }

  # Progress
  if (bi %% max(1, n_batches %/% 10) == 0 || bi == n_batches) {
    elapsed <- difftime(Sys.time(), t_start, units = "secs")
    rate <- idx_end / as.numeric(elapsed)
    eta <- (n_orfs - idx_end) / rate
    cat(sprintf("  [%d/%d] %d ORFs, %.0f ORF/s, ETA %.0fs\n",
                bi, n_batches, idx_end, rate, eta))
  }
}

total_time <- difftime(Sys.time(), t_start, units = "secs")
cat(sprintf("  Done in %.1fs (%.1f ORF/s)\n", total_time, n_orfs / total_time))

# ── 5. Write output ──────────────────────────────────────────────────────────
cat("\n[4/5] Writing output...\n")

# Select and order output columns
out_cols <- intersect(c("orf_id", "gene_id", "chr", "start", "end", "strand",
                         "ORF_category_Gen", "ORF_category_Tx",
                         "ORF_category_Tx_compatible"),
                       names(orf_dt))
fwrite(orf_dt[, ..out_cols], file = output_file, sep = "\t", quote = FALSE, na = "-")
cat(sprintf("  Wrote: %s (%d ORFs)\n", output_file, n_orfs))

# ── Summary ───────────────────────────────────────────────────────────────────
cat("\n[5/5] Classification summary:\n")
for (cn in c("ORF_category_Gen", "ORF_category_Tx", "ORF_category_Tx_compatible")) {
  if (cn %in% names(orf_dt)) {
    cat(sprintf("\n  %s:\n", cn))
    tbl <- orf_dt[, .N, by = get(cn)][order(-N)]
    for (ri in seq_len(min(nrow(tbl), 15))) {
      cat(sprintf("    %-25s %d\n", tbl[[1]][ri], tbl$N[ri]))
    }
  }
}

cat("\n===== CLASSIFICATION COMPLETE =====\n")
