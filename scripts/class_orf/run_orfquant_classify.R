#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(rtracklayer))

# Locate the core function script
script_dir <- dirname(sub("--file=", "", commandArgs(trailingOnly=FALSE)[grep("--file=", commandArgs(trailingOnly=FALSE))]))
if (length(script_dir) == 0) script_dir <- "."
source(file.path(script_dir, "orfquant_orf_classify.R"))

option_list <- list(
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="Input ORFs (GTF/GFF/BED file)", metavar="file"),
  make_option(c("-a", "--annotation"), type="character", default=NULL,
              help="Reference annotation (GTF/GFF file)", metavar="file"),
  make_option(c("-o", "--output"), type="character", default=NULL,
              help="Output classification TSV file", metavar="file")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input) || is.null(opt$annotation) || is.null(opt$output)) {
  print_help(opt_parser)
  stop("Missing required arguments", call.=FALSE)
}

message("Loading ORFs from: ", opt$input)
# Import ORFs
# If GTF, we import it.
# The unified tool outputs GTF where features are 'CDS' or 'exon'.
# orfquant_classify_orfs expects GRanges or GRangesList or file path.
# If we pass file path, rtracklayer::import handles it.

message("Loading Annotation from: ", opt$annotation)
# orfquant_classify_orfs handles file path for annotation too.

message("Running classification...")
# Call the function
# Note: unified GTF has 'orf_id' in attributes, which will become a column in mcols(gr).
# The function uses 'orf_id_col' parameter if provided, or row names.
# Let's see how `orfquant_classify_orfs` handles it.
# It seems `orfquant_classify_orfs` (based on previous `view`) takes `orfs` and `annotation`.

# If we import input first, we can control how it is structured
orfs_gr <- rtracklayer::import(opt$input)
# The unified GTF uses 'CDS' and 'exon'.
# We should filter for CDS usually for ORFs, but let's see.
# The unified GTF has both.
orfs_gr <- orfs_gr[orfs_gr$type == "CDS"]

# We need to split into GRangesList by orf_id
orfs_grl <- split(orfs_gr, orfs_gr$orf_id)

results <- orfquant_classify_orfs(
  orfs = orfs_grl,
  annotation = opt$annotation
)

message("Writing results to: ", opt$output)
write.table(results, file = opt$output, sep = "\t", quote = FALSE, row.names = FALSE)
message("Done.")
