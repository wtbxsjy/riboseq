#!/usr/bin/env Rscript
#
# orfquant_classify_light.R — Lightweight ORFquant-style classification.
# Uses data.table only (no Bioconductor / GRanges).
# Target: ~100 ORF/s, ~80 min for 500K ORFs (vs 16+ hr with GRanges).
#
# Usage:
#   Rscript orfquant_classify_light.R \
#     --orfs unified_orfs.metadata.tsv \
#     --gtf annotation.gtf \
#     --output orfquant_classification.tsv
#

suppressPackageStartupMessages({
  if (!requireNamespace("data.table", quietly = TRUE))
    install.packages("data.table", repos = "https://cloud.r-project.org", quiet = TRUE)
  library(data.table)
})

# ── CLI ──────────────────────────────────────────────────────────────────────
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

cat("=== ORFquant Classification (data.table) ===\n")
cat(sprintf("  ORFs: %s\n  GTF : %s\n", orfs_file, gtf_file))

# ── Helpers ──────────────────────────────────────────────────────────────────

parse_exon_blocks <- function(s) {
  if (is.na(s) || !nzchar(s) || s == "None") return(NULL)
  parts <- strsplit(s, ",")[[1]]
  do.call(rbind, lapply(parts, function(p) {
    se <- as.integer(strsplit(p, "-")[[1]])
    if (length(se) != 2) NULL else data.table(start = se[1], end = se[2])
  }))
}

# Transcript-coordinate projection: map genomic block(s) through exon chain.
# Returns c(tx_start, tx_end) or NULL.
project_to_tx <- function(blocks, tx_exons, strand) {
  if (is.null(blocks) || nrow(blocks) == 0 || nrow(tx_exons) == 0) return(NULL)
  tx_start <- NA_integer_; tx_end <- NA_integer_
  for (bi in seq_len(nrow(blocks))) {
    bs <- blocks$start[bi]; be <- blocks$end[bi]
    hits <- tx_exons[start <= be & end >= bs]
    if (nrow(hits) == 0) return(NULL)
    if (strand == "+") {
      bs_tx <- hits$cum_len[1] + (bs - hits$start[1] + 1L)
      be_tx <- hits$cum_len[nrow(hits)] + (be - hits$start[nrow(hits)] + 1L)
    } else {
      bs_tx <- hits$cum_len[1] + (hits$end[1] - bs + 1L)
      be_tx <- hits$cum_len[nrow(hits)] + (hits$end[nrow(hits)] - be + 1L)
    }
    if (is.na(tx_start) || bs_tx < tx_start) tx_start <- bs_tx
    if (is.na(tx_end) || be_tx > tx_end) tx_end <- be_tx
  }
  if (is.na(tx_start) || is.na(tx_end)) return(NULL)
  if (tx_start > tx_end) { tmp <- tx_start; tx_start <- tx_end; tx_end <- tmp }
  c(tx_start, tx_end)
}

# ── 1. Load & index GTF annotation ───────────────────────────────────────────

cat("\n[1/4] Loading annotation...\n")
rds_cache <- file.path(cache_dir, paste0(basename(gtf_file), ".orfclass_cache.rds"))
if (file.exists(rds_cache) && file.mtime(rds_cache) >= file.mtime(gtf_file)) {
  cat("  Loading cached annotation:", rds_cache, "\n")
  ann <- readRDS(rds_cache)
} else {
  cat("  Reading GTF with fread...\n")
  gtf <- fread(gtf_file, sep = "\t", header = FALSE, skip = "#",
               select = c(1,3,4,5,7,9), showProgress = FALSE,
               col.names = c("chr","type","start","end","strand","attrs"))
  gtf <- gtf[type %in% c("exon","CDS","gene")]
  cat(sprintf("  %s rows loaded\n", format(nrow(gtf), big.mark=",")))

  gtf[, gene_id       := sub('.*gene_id "([^"]+)".*', '\\1', attrs)]
  gtf[, transcript_id := sub('.*transcript_id "([^"]+)".*', '\\1', attrs)]
  gtf[, gene_biotype  := sub('.*gene_biotype "([^"]+)".*', '\\1', attrs)]
  gtf[gene_id == attrs,       gene_id := NA]
  gtf[transcript_id == attrs, transcript_id := NA]
  gtf[gene_biotype == attrs,  gene_biotype := NA]
  gtf[, attrs := NULL]

  # Exons per transcript with cum_len for coordinate projection
  exon_dt <- gtf[type == "exon" & !is.na(transcript_id)]
  setorder(exon_dt, transcript_id, start)
  exon_dt[, exon_rank := seq_len(.N), by = transcript_id]
  exon_dt[, exon_width := end - start + 1L]
  exon_dt[, cum_len := as.integer(cumsum(c(0L, exon_width[-.N]))), by = transcript_id]
  exon_by_tx <- split(exon_dt, exon_dt$transcript_id)

  # CDS per transcript
  cds_dt <- gtf[type == "CDS" & !is.na(transcript_id),
    .(transcript_id, gene_id, chr, start, end, strand)]
  setorder(cds_dt, transcript_id, start)
  cds_merged <- cds_dt[, .(
    cds_start = min(start), cds_end = max(end),
    cds_width = sum(end - start + 1L)
  ), by = .(transcript_id, gene_id, chr, strand)]
  setkey(cds_merged, transcript_id)

  # Max CDS per gene
  max_cds_by_gene <- cds_merged[!is.na(cds_width)][, .SD[which.max(cds_width)], by = gene_id]
  setkey(max_cds_by_gene, gene_id)

  # Gene-level CDS bounds
  cds_gene <- cds_merged[!is.na(cds_start), .(
    cds_gene_start = min(cds_start), cds_gene_end = max(cds_end)
  ), by = .(gene_id, chr, strand)]
  setkey(cds_gene, gene_id, chr, strand)

  # all-CDS bounds (for ORFs without gene assignment)
  all_cds_start <- min(cds_merged$cds_start, na.rm = TRUE)
  all_cds_end   <- max(cds_merged$cds_end,   na.rm = TRUE)

  # LncRNA genes
  lncrna_genes <- gtf[type == "gene" & !is.na(gene_biotype) &
    gene_biotype %in% c("lncRNA","lincRNA","antisense","processed_transcript",
      "sense_intronic","sense_overlapping","non_coding",
      "3prime_overlapping_ncrna","bidirectional_promoter_lncrna",
      "ncRNA","antisense_RNA","misc_non_coding"), unique(gene_id)]

  # Gene → transcript list
  gene_to_txids <- split(cds_merged$transcript_id, cds_merged$gene_id)

  rm(gtf, cds_dt); gc()

  ann <- list(exon_by_tx = exon_by_tx, cds_merged = cds_merged,
              max_cds_by_gene = max_cds_by_gene, cds_gene = cds_gene,
              all_cds_start = all_cds_start, all_cds_end = all_cds_end,
              lncrna_genes = lncrna_genes, gene_to_txids = gene_to_txids)
  saveRDS(ann, rds_cache)
  cat("  Cached:", rds_cache, "\n")
}

# ── 2. Load ORFs ─────────────────────────────────────────────────────────────

cat("\n[2/4] Loading ORFs...\n")
orfs <- fread(orfs_file, sep = "\t", na.strings = c("", "NA", "None"))
n_orfs <- nrow(orfs)
cat(sprintf("  %d ORFs loaded\n", n_orfs))

# Normalize column names: metadata uses 'chrom'
setnames(orfs, "chrom", "chr", skip_absent = TRUE)

# Parse exon blocks into list of data.tables
cat("  Parsing exon blocks...\n")
exon_list <- lapply(orfs$exon_blocks, parse_exon_blocks)

# Fill NAs
orfs[is.na(gene_id),       gene_id := ""]
orfs[is.na(transcript_id), transcript_id := ""]
orfs[is.na(strand),        strand := "*"]

# ── 3. Classification (vectorized) ───────────────────────────────────────────

cat("\n[3/4] Classifying...\n")
t_start <- Sys.time()

# ---- 3a. Genomic classification (vectorized) ----

# Quick gene→CDS lookup (indexed for fast per-gene access)
setkey(ann$max_cds_by_gene, gene_id)
cds_gene_list  <- split(ann$cds_gene, by = "gene_id")
cds_gene_list2 <- split(ann$max_cds_by_gene, by = "gene_id")

find_gen_class <- function(n) {
  res <- rep("novel", n)

  # Process by unique gene_id in ORFs
  orf_genes <- unique(orfs$gene_id[nzchar(orfs$gene_id)])

  for (gid in orf_genes) {
    cg_list <- cds_gene_list[[gid]]
    if (is.null(cg_list) || nrow(cg_list) == 0) next

    for (ri in seq_len(nrow(cg_list))) {
      cg_row <- cg_list[ri]
      cds_s <- cg_row$cds_gene_start
      cds_e <- cg_row$cds_gene_end
      cg_chr <- cg_row$chr
      cg_strand <- cg_row$strand

      idx <- which(orfs$gene_id == gid & orfs$chr == cg_chr & orfs$strand == cg_strand)
      if (length(idx) == 0) next

      orf_s <- orfs$start[idx]; orf_e <- orfs$end[idx]
      orf_strand <- orfs$strand[idx]

      # Overlap check
      ov <- orf_s <= cds_e & orf_e >= cds_s
      if (any(ov)) {
        # Overlapping CDS — classify relative to max CDS isoform
        if (gid %in% names(cds_gene_list2)) {
          mc_list <- cds_gene_list2[[gid]]
          if (!is.null(mc_list) && nrow(mc_list) > 0) {
            for (i in which(ov)[1]) {
              mcg <- mc_list[1]
              mc_s <- mcg$cds_start; mc_e <- mcg$cds_end
              strd <- orf_strand[i]
              if (strd == "+") {
                gs <- mc_s; ge <- mc_e; so <- orf_s[i]; eo <- orf_e[i] + 3L
              } else {
                gs <- mc_e; ge <- mc_s; so <- orf_e[i]; eo <- orf_s[i] - 3L
              }
              if (eo == ge) {
                res[idx[i]] <- if (so == gs) "exact_start_stop"
                  else if (so < gs) { if (strd == "+") "Alt5_start" else "Alt3_start" }
                  else { if (strd == "+") "Alt3_start" else "Alt5_start" }
              } else {
                res[idx[i]] <- if (strd == "+") {
                  if (so < gs && eo < ge) "Alt5_start_Alt5_stop"
                  else if (so < gs && eo > ge) "Alt5_start_Alt3_stop"
                  else if (so > gs && eo > ge) "Alt3_start_Alt3_stop"
                  else if (so > gs && eo < ge) "Alt3_start_Alt5_stop"
                  else if (so == gs && eo < ge) "Alt5_stop"
                  else if (so == gs && eo > ge) "Alt3_stop"
                  else "overlaps_CDS"
                } else {
                  if (so > gs && eo > ge) "Alt5_start_Alt5_stop"
                  else if (so > gs && eo < ge) "Alt5_start_Alt3_stop"
                  else if (so < gs && eo < ge) "Alt3_start_Alt3_stop"
                  else if (so < gs && eo > ge) "Alt3_start_Alt5_stop"
                  else if (so == gs && eo > ge) "Alt5_stop"
                  else if (so == gs && eo < ge) "Alt3_stop"
                  else "overlaps_CDS"
                }
              }
              # Other overlapping ORFs for same gene get generic label
              ov_others <- setdiff(which(ov)[-1], i)
              if (length(ov_others) > 0) res[idx[ov_others]] <- "overlaps_CDS"
            }
            next
          }
        }
        res[idx[ov]] <- "overlaps_CDS"
      }

      # Non-overlapping
      non_ov <- which(!ov)
      if (length(non_ov) > 0) {
        for (j in non_ov) {
          strd <- orf_strand[j]
          if (strd == "+") {
            if (orf_e[j] < cds_s) res[idx[j]] <- "novel_Upstream"
            else if (orf_s[j] > cds_e) res[idx[j]] <- "novel_Downstream"
            else res[idx[j]] <- "novel_Internal"
          } else if (strd == "-") {
            if (orf_e[j] > cds_e) res[idx[j]] <- "novel_Upstream"
            else if (orf_s[j] < cds_s) res[idx[j]] <- "novel_Downstream"
            else res[idx[j]] <- "novel_Internal"
          }
        }
      }
    }
  }

  # LncRNA override
  if (length(ann$lncrna_genes) > 0) {
    res[orfs$gene_id %in% ann$lncrna_genes] <- "lncRNA"
  }

  res
}

cat("  Genomic classification...\n")
gen_class <- find_gen_class(n_orfs)
cat(sprintf("  Genomic: %.1fs\n", difftime(Sys.time(), t_start, units = "secs")))

# ---- 3b. Transcript classification ----
# Only for ORFs with valid transcript_id matching an annotated transcript

cat("  Transcript classification...\n")
t_tx <- Sys.time()

tx_class  <- rep(NA_character_, n_orfs)
tx_compat <- rep(NA_character_, n_orfs)

# TX class priority (lower = better)
tx_prio <- c(ORF_annotated=1L, N_extension=2L, C_extension=3L,
             N_truncation=4L, C_truncation=5L, NC_extension=6L,
             overl_uORF=7L, overl_dORF=8L, nested_ORF=9L,
             uORF=10L, dORF=11L, novel=12L, novel_antisense=13L)

# Process ORFs that have a transcript_id matching the annotation
has_tx <- nzchar(orfs$transcript_id) & orfs$transcript_id %in% names(ann$exon_by_tx)
cat(sprintf("  %d/%d ORFs have matching transcript\n", sum(has_tx), n_orfs))

tx_ids  <- orfs$transcript_id[has_tx]
uniq_tx <- unique(tx_ids)
cat(sprintf("  %d unique transcripts\n", length(uniq_tx)))

# For each unique transcript, process all ORFs that map to it
for (ti in seq_along(uniq_tx)) {
  tx_id <- uniq_tx[ti]
  idx   <- which(has_tx & orfs$transcript_id == tx_id)
  if (length(idx) == 0) next

  tx_exons <- ann$exon_by_tx[[tx_id]]
  cd <- ann$cds_merged[transcript_id == tx_id]
  has_cds <- nrow(cd) > 0 && !is.na(cd$cds_start[1])

  for (i in idx) {
    blocks <- exon_list[[i]]
    if (is.null(blocks)) next

    orf_tx <- project_to_tx(blocks, tx_exons, orfs$strand[i])
    if (is.null(orf_tx)) next

    # ORF_category_Tx
    if (has_cds) {
      o_s <- orf_tx[1]; o_e <- orf_tx[2]
      a_s <- cd$cds_start[1]; a_e <- cd$cds_end[1]

      if (o_e == a_e) {
        tx_class[i] <- if (o_s == a_s) "ORF_annotated"
          else if (o_s < a_s) "N_extension"
          else "N_truncation"
      } else {
        tx_class[i] <- if (o_s < a_s && o_e < a_s) "uORF"
          else if (o_s < a_s && o_e < a_e) "overl_uORF"
          else if (o_s < a_s && o_e > a_e) "NC_extension"
          else if (o_s > a_e && o_e > a_e) "dORF"
          else if (o_s > a_s && o_e > a_e) "overl_dORF"
          else if (o_s > a_s && o_e < a_e) "nested_ORF"
          else if (o_s == a_s && o_e < a_e) "C_truncation"
          else if (o_s == a_s && o_e > a_e) "C_extension"
      }
    } else {
      tx_class[i] <- "novel"
    }
  }

  # Progress
  if (ti %% max(1, length(uniq_tx) %/% 10) == 0 || ti == length(uniq_tx)) {
    cat(sprintf("    Tx: %d/%d (%.0f ORF/s)\n", ti, length(uniq_tx),
                sum(has_tx) / as.numeric(difftime(Sys.time(), t_tx, units = "secs"))))
  }
}

cat(sprintf("  Transcript: %.1fs\n", difftime(Sys.time(), t_tx, units = "secs")))

# ---- 3c. Best-isoform classification (Tx_compatible) ----

cat("  Best-isoform classification...\n")
t_comp <- Sys.time()

# For each ORF with a gene_id, try all transcripts of that gene
has_gene <- nzchar(orfs$gene_id) & orfs$gene_id %in% names(ann$gene_to_txids)
cat(sprintf("  %d/%d ORFs have matching gene\n", sum(has_gene), n_orfs))

# Process in gene batches to amortize cache lookups
uniq_genes <- unique(orfs$gene_id[has_gene])
n_genes <- length(uniq_genes)
processed_orf <- 0L
for (gi in seq_along(uniq_genes)) {
  gid <- uniq_genes[gi]
  idx <- which(has_gene & orfs$gene_id == gid)
  if (length(idx) == 0) next
  processed_orf <- processed_orf + length(idx)

  tx_list <- ann$gene_to_txids[[gid]]
  if (length(tx_list) == 0) next

  # For each transcript of the gene, try to classify all ORFs
  for (i in idx) {
    blocks <- exon_list[[i]]
    if (is.null(blocks)) next
    best_cls <- "novel"; best_pri <- tx_prio["novel"]

    for (txc_id in tx_list) {
      if (!txc_id %in% names(ann$exon_by_tx)) next
      cd <- ann$cds_merged[transcript_id == txc_id]
      if (nrow(cd) == 0 || is.na(cd$cds_start[1])) next

      orf_txc <- project_to_tx(blocks, ann$exon_by_tx[[txc_id]], orfs$strand[i])
      if (is.null(orf_txc)) next

      o_s <- orf_txc[1]; o_e <- orf_txc[2]; a_s <- cd$cds_start[1]; a_e <- cd$cds_end[1]
      cls <- NA_character_
      if (o_e == a_e) {
        cls <- if (o_s == a_s) "ORF_annotated" else if (o_s < a_s) "N_extension" else "N_truncation"
      } else {
        cls <- if (o_s < a_s && o_e < a_s) "uORF"
          else if (o_s < a_s && o_e < a_e) "overl_uORF"
          else if (o_s < a_s && o_e > a_e) "NC_extension"
          else if (o_s > a_e && o_e > a_e) "dORF"
          else if (o_s > a_s && o_e > a_e) "overl_dORF"
          else if (o_s > a_s && o_e < a_e) "nested_ORF"
          else if (o_s == a_s && o_e < a_e) "C_truncation"
          else if (o_s == a_s && o_e > a_e) "C_extension"
      }
      if (!is.na(cls)) {
        pri <- tx_prio[cls]
        if (is.na(pri)) pri <- tx_prio["novel"]
        if (pri < best_pri) { best_pri <- pri; best_cls <- cls }
      }
      if (best_pri <= 1L) break
    }
    tx_compat[i] <- best_cls
  }

  if (gi %% max(1, n_genes %/% 10) == 0 || gi == n_genes) {
    elapsed <- as.numeric(difftime(Sys.time(), t_comp, units = "secs"))
    cat(sprintf("    Compat: %d/%d genes, %d ORFs (%.0f ORF/s)\n",
                gi, n_genes, processed_orf, processed_orf / max(1, elapsed)))
  }
}

cat(sprintf("  Compat: %.1fs\n", difftime(Sys.time(), t_comp, units = "secs")))

total_time <- difftime(Sys.time(), t_start, units = "secs")
cat(sprintf("  Total classification: %.1fs (%.1f ORF/s)\n",
            as.numeric(total_time), n_orfs / as.numeric(total_time)))

# ── 4. Write output ──────────────────────────────────────────────────────────

cat("\n[4/4] Writing output...\n")

out <- data.table(
  orf_id                     = orfs$orf_id,
  gene_id                    = orfs$gene_id,
  transcript_id              = orfs$transcript_id,
  chr                        = orfs$chr,
  start                      = orfs$start,
  end                        = orfs$end,
  strand                     = orfs$strand,
  ORF_category_Gen           = gen_class,
  ORF_category_Tx            = tx_class,
  ORF_category_Tx_compatible = tx_compat
)

fwrite(out, file = output_file, sep = "\t", quote = FALSE, na = "-")
cat(sprintf("  Wrote: %s (%d ORFs)\n", output_file, nrow(out)))

# Summary
cat("\n=== Classification Summary ===\n")
for (cn in c("ORF_category_Gen", "ORF_category_Tx", "ORF_category_Tx_compatible")) {
  tbl <- out[!is.na(get(cn)), .N, by = get(cn)][order(-N)]
  cat(sprintf("\n%s:\n", cn))
  for (ri in seq_len(min(nrow(tbl), 15))) {
    cat(sprintf("  %-25s %d\n", tbl[[1]][ri], tbl$N[ri]))
  }
}

cat("\n===== DONE =====\n")
