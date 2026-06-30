#!/usr/bin/env Rscript
#############################################################################
# ORF Quantification from RiboseQC P-site Bedgraphs
#
# Purpose: Quantify unified GENCODE-annotated ORFs using RiboseQC P-site data
# Input:
#   - GENCODE ORF GTF (from gencode-riboseqORFs)
#   - RiboseQC P-site bedgraph files (plus/minus strands)
# Output:
#   - Count matrix (ORFs × samples)
#   - TPM normalized matrix
#   - Summary statistics
#
# Usage:
#   Rscript quantify_orfs_from_psites.R \
#     --gtf results/Mouse_Unified.orfs.gtf \
#     --bedgraph-dir riboseqc_results/ \
#     --sample-pattern "(.+)_P_sites_unique_(plus|minus)\\.bedgraph$" \
#     --outdir quantification_results
#############################################################################

suppressPackageStartupMessages({
  library(optparse)
  library(GenomicRanges)
  library(rtracklayer)
  library(tidyverse)
  library(data.table)
})

# ========== Command Line Arguments ==========
option_list <- list(
  make_option(c("--gtf"), type = "character",
              help = "GENCODE ORF GTF file (from gencode-riboseqORFs)"),
  make_option(c("--bedgraph-dir"), type = "character",
              help = "Directory containing RiboseQC bedgraph files"),
  make_option(c("--sample-pattern"), type = "character",
              default = "(.+)_P_sites_unique_(plus|minus)\\.bedgraph$",
              help = "Regex pattern to extract sample names [default: %default]"),
  make_option(c("--outdir"), type = "character", default = "quantification_results",
              help = "Output directory [default: %default]"),
  make_option(c("--min-count"), type = "integer", default = 5,
              help = "Minimum total P-site count to keep an ORF [default: %default]"),
  make_option(c("--threads"), type = "integer", default = 4,
              help = "Number of threads for parallel processing [default: %default]")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Validate inputs
if (is.null(opt$gtf) || is.null(opt$`bedgraph-dir`)) {
  print_help(opt_parser)
  stop("Both --gtf and --bedgraph-dir are required", call. = FALSE)
}

if (!file.exists(opt$gtf)) {
  stop("GTF file not found: ", opt$gtf, call. = FALSE)
}

if (!dir.exists(opt$`bedgraph-dir`)) {
  stop("Bedgraph directory not found: ", opt$`bedgraph-dir`, call. = FALSE)
}

# Create output directory
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

cat("========================================\n")
cat("ORF Quantification from P-site Bedgraphs\n")
cat("========================================\n\n")

# ========== Step 1: Load GENCODE ORF GTF ==========
cat("[Step 1/5] Loading GENCODE ORF GTF...\n")

# Import GTF as GRanges
orf_gr <- import(opt$gtf)

# For gencode-riboseqORFs output, the ORF features are typically "CDS" or "exon"
# We want the ORF-level annotation (use "transcript" if available, otherwise collapse CDS)
if ("transcript" %in% mcols(orf_gr)$type) {
  orf_gr <- orf_gr[mcols(orf_gr)$type == "transcript"]
} else {
  # Collapse by transcript_id to get ORF-level ranges
  orf_gr <- reduce(split(orf_gr, mcols(orf_gr)$transcript_id))
  orf_gr <- unlist(orf_gr)
  names(mcols(orf_gr)) <- NULL
}

# Extract ORF metadata
orf_metadata <- mcols(orf_gr) %>%
  as_tibble() %>%
  mutate(
    orf_id = transcript_id,  # gencode-riboseqORFs uses transcript_id for ORF ID
    seqname = as.character(seqnames(orf_gr)),
    start = start(orf_gr),
    end = end(orf_gr),
    strand = as.character(strand(orf_gr)),
    width = width(orf_gr)
  )

cat(sprintf("  ✓ Loaded %d ORFs\n", length(orf_gr)))
cat(sprintf("  ✓ ORF biotypes: %s\n",
            paste(unique(orf_metadata$gene_type), collapse = ", ")))

# ========== Step 2: Discover and Load P-site Bedgraphs ==========
cat("\n[Step 2/5] Discovering P-site bedgraph files...\n")

# Find all bedgraph files
bedgraph_files <- list.files(
  opt$`bedgraph-dir`,
  pattern = opt$`sample-pattern`,
  full.names = TRUE
)

if (length(bedgraph_files) == 0) {
  stop("No bedgraph files found matching pattern: ", opt$`sample-pattern`, call. = FALSE)
}

# Parse filenames to extract sample names and strands
bedgraph_info <- tibble(path = bedgraph_files) %>%
  mutate(
    filename = basename(path),
    sample = str_match(filename, opt$`sample-pattern`)[, 2],
    strand_label = str_match(filename, opt$`sample-pattern`)[, 3]
  ) %>%
  filter(!is.na(sample)) %>%
  mutate(
    strand = case_when(
      str_detect(strand_label, "plus") ~ "+",
      str_detect(strand_label, "minus") ~ "-",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(strand))

# Check for paired plus/minus files
sample_counts <- bedgraph_info %>%
  group_by(sample) %>%
  summarise(n_strands = n(), .groups = "drop")

if (any(sample_counts$n_strands != 2)) {
  warning("Some samples missing plus or minus bedgraph files")
  cat("  Sample strand counts:\n")
  print(sample_counts)
}

samples <- unique(bedgraph_info$sample)
cat(sprintf("  ✓ Found %d samples: %s\n",
            length(samples),
            paste(head(samples, 5), collapse = ", ")))

if (length(samples) > 5) {
  cat(sprintf("    ... and %d more\n", length(samples) - 5))
}

# ========== Step 3: Quantify P-sites per ORF (Strand-Specific) ==========
cat("\n[Step 3/5] Quantifying P-sites per ORF (strand-specific)...\n")

# Function to load bedgraph and convert to GRanges
load_bedgraph_as_gr <- function(path, strand_char) {
  # Read bedgraph (chr, start, end, count)
  # Bedgraph is 0-based, half-open [start, end)
  bg <- fread(path, col.names = c("seqnames", "start", "end", "score"))

  # Convert to GRanges (GenomicRanges uses 1-based coordinates internally)
  # For bedgraph: add 1 to start to convert 0-based to 1-based
  gr <- GRanges(
    seqnames = bg$seqnames,
    ranges = IRanges(start = bg$start + 1, end = bg$end),
    strand = strand_char,
    score = bg$score
  )

  return(gr)
}

# Function to count P-sites overlapping ORFs (same strand only)
count_psites_in_orfs <- function(orf_gr, psite_gr) {
  # Find overlaps (strand-specific)
  overlaps <- findOverlaps(psite_gr, orf_gr, type = "within", ignore.strand = FALSE)

  # Sum P-site scores per ORF
  counts <- tapply(
    mcols(psite_gr)$score[queryHits(overlaps)],
    subjectHits(overlaps),
    sum
  )

  # Create full count vector (0 for ORFs with no P-sites)
  full_counts <- rep(0, length(orf_gr))
  full_counts[as.integer(names(counts))] <- counts

  return(full_counts)
}

# Initialize count matrix
count_matrix <- matrix(
  0,
  nrow = length(orf_gr),
  ncol = length(samples),
  dimnames = list(orf_metadata$orf_id, samples)
)

# Process each sample
pb <- txtProgressBar(min = 0, max = length(samples), style = 3)
for (i in seq_along(samples)) {
  sample_name <- samples[i]

  # Get plus and minus bedgraph paths
  plus_file <- bedgraph_info %>%
    filter(sample == sample_name, strand == "+") %>%
    pull(path)

  minus_file <- bedgraph_info %>%
    filter(sample == sample_name, strand == "-") %>%
    pull(path)

  # Load bedgraphs
  if (length(plus_file) == 1) {
    plus_gr <- load_bedgraph_as_gr(plus_file, "+")
    plus_counts <- count_psites_in_orfs(orf_gr, plus_gr)
  } else {
    plus_counts <- rep(0, length(orf_gr))
  }

  if (length(minus_file) == 1) {
    minus_gr <- load_bedgraph_as_gr(minus_file, "-")
    minus_counts <- count_psites_in_orfs(orf_gr, minus_gr)
  } else {
    minus_counts <- rep(0, length(orf_gr))
  }

  # Combine strand-specific counts (each ORF has only one strand)
  count_matrix[, i] <- plus_counts + minus_counts

  setTxtProgressBar(pb, i)
}
close(pb)

cat("\n  ✓ Quantification complete\n")

# ========== Step 4: Filter and Normalize ==========
cat("\n[Step 4/5] Filtering and normalizing...\n")

# Filter ORFs with low total counts
total_counts <- rowSums(count_matrix)
keep_orfs <- total_counts >= opt$`min-count`

cat(sprintf("  ✓ Keeping %d / %d ORFs (≥%d total P-sites)\n",
            sum(keep_orfs), length(keep_orfs), opt$`min-count`))

count_matrix_filtered <- count_matrix[keep_orfs, ]
orf_metadata_filtered <- orf_metadata[keep_orfs, ]

# Calculate TPM (Transcripts Per Million)
# TPM = (counts / ORF_length) / sum(counts / ORF_length) * 1e6
orf_lengths <- orf_metadata_filtered$width / 1000  # Convert to kb

tpm_matrix <- sweep(count_matrix_filtered, 1, orf_lengths, FUN = "/")
tpm_matrix <- sweep(tpm_matrix, 2, colSums(tpm_matrix), FUN = "/") * 1e6

cat("  ✓ TPM normalization complete\n")

# ========== Step 5: Save Results ==========
cat("\n[Step 5/5] Saving results...\n")

# Save raw counts
count_df <- count_matrix_filtered %>%
  as_tibble(rownames = "orf_id") %>%
  left_join(
    orf_metadata_filtered %>% select(orf_id, seqname, start, end, strand, gene_name, gene_type),
    by = "orf_id"
  ) %>%
  select(orf_id, seqname, start, end, strand, gene_name, gene_type, everything())

write_tsv(count_df, file.path(opt$outdir, "orf_counts_raw.tsv"))
cat(sprintf("  ✓ Raw counts: %s\n", file.path(opt$outdir, "orf_counts_raw.tsv")))

# Save TPM-normalized counts
tpm_df <- tpm_matrix %>%
  as_tibble(rownames = "orf_id") %>%
  left_join(
    orf_metadata_filtered %>% select(orf_id, seqname, start, end, strand, gene_name, gene_type),
    by = "orf_id"
  ) %>%
  select(orf_id, seqname, start, end, strand, gene_name, gene_type, everything())

write_tsv(tpm_df, file.path(opt$outdir, "orf_counts_tpm.tsv"))
cat(sprintf("  ✓ TPM counts: %s\n", file.path(opt$outdir, "orf_counts_tpm.tsv")))

# Save count matrix only (for DESeq2)
count_matrix_only <- count_matrix_filtered
write.csv(count_matrix_only, file.path(opt$outdir, "orf_counts_matrix.csv"))
cat(sprintf("  ✓ Count matrix: %s\n", file.path(opt$outdir, "orf_counts_matrix.csv")))

# Save summary statistics
summary_stats <- tibble(
  sample = colnames(count_matrix_filtered),
  total_psites = colSums(count_matrix_filtered),
  n_orfs_detected = colSums(count_matrix_filtered > 0),
  median_counts = apply(count_matrix_filtered, 2, function(x) median(x[x > 0])),
  mean_counts = colMeans(count_matrix_filtered)
)

write_tsv(summary_stats, file.path(opt$outdir, "sample_summary_stats.tsv"))
cat(sprintf("  ✓ Summary stats: %s\n", file.path(opt$outdir, "sample_summary_stats.tsv")))

# Save ORF-level statistics
orf_stats <- tibble(
  orf_id = orf_metadata_filtered$orf_id,
  seqname = orf_metadata_filtered$seqname,
  start = orf_metadata_filtered$start,
  end = orf_metadata_filtered$end,
  strand = orf_metadata_filtered$strand,
  width = orf_metadata_filtered$width,
  gene_name = orf_metadata_filtered$gene_name,
  gene_type = orf_metadata_filtered$gene_type,
  total_psites = rowSums(count_matrix_filtered),
  n_samples_detected = rowSums(count_matrix_filtered > 0),
  mean_tpm = rowMeans(tpm_matrix),
  max_tpm = apply(tpm_matrix, 1, max)
)

write_tsv(orf_stats, file.path(opt$outdir, "orf_summary_stats.tsv"))
cat(sprintf("  ✓ ORF stats: %s\n", file.path(opt$outdir, "orf_summary_stats.tsv")))

# ========== Session Info ==========
cat("\n========================================\n")
cat("✅ Quantification Complete!\n")
cat("========================================\n\n")
cat("Output files:\n")
cat(sprintf("  - Raw counts:       %s\n", "orf_counts_raw.tsv"))
cat(sprintf("  - TPM counts:       %s\n", "orf_counts_tpm.tsv"))
cat(sprintf("  - Count matrix:     %s\n", "orf_counts_matrix.csv"))
cat(sprintf("  - Sample stats:     %s\n", "sample_summary_stats.tsv"))
cat(sprintf("  - ORF stats:        %s\n", "orf_summary_stats.tsv"))
cat("\nNext steps:\n")
cat("  1. Differential expression: Use orf_counts_matrix.csv with DESeq2\n")
cat("  2. Visualization: Load orf_counts_tpm.tsv for heatmaps\n")
cat("  3. ORF filtering: Use orf_summary_stats.tsv to filter by detection\n")

# Save session info
sink(file.path(opt$outdir, "session_info.txt"))
sessionInfo()
sink()
