suppressPackageStartupMessages({
  library(GenomicRanges)
  library(IRanges)
  library(S4Vectors)
})

cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- cmd_args[grepl("^--file=", cmd_args)]
script_path <- if (length(file_arg) > 0) sub("^--file=", "", file_arg[1]) else ""
script_dir <- if (nzchar(script_path)) dirname(normalizePath(script_path)) else getwd()

source(file.path(dirname(script_dir), "orfquant_orf_classify.R"))

parse_coords <- function(coord_str) {
  if (is.na(coord_str) || coord_str == "") return(list())
  parts <- unlist(strsplit(coord_str, ","))
  lapply(parts, function(p) {
    se <- unlist(strsplit(p, "-"))
    if (length(se) != 2) return(NULL)
    c(as.integer(se[1]), as.integer(se[2]))
  })
}

cds_path <- file.path(script_dir, "test_cds.txt")
cds_df <- read.table(cds_path, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
colnames(cds_df) <- c("gene_id", "dummy", "coords", "strand")

gene_to_seq <- c(geneA = "chr1", geneB = "chr2", geneC = "chr3")

cds_genes <- GRangesList()
cds_txs <- GRangesList()

for (i in seq_len(nrow(cds_df))) {
  gene_id <- cds_df$gene_id[i]
  strand <- cds_df$strand[i]
  coords <- parse_coords(cds_df$coords[i])
  coords <- coords[!vapply(coords, is.null, logical(1))]
  starts <- vapply(coords, function(v) v[1], numeric(1))
  ends <- vapply(coords, function(v) v[2], numeric(1))
  seqname <- gene_to_seq[[gene_id]]
  gr <- GRanges(seqnames = seqname, ranges = IRanges(start = starts, end = ends), strand = strand)
  gr$gene_id <- gene_id
  gr$transcript_id <- gene_id
  cds_genes[[gene_id]] <- reduce(gr)
  cds_txs[[gene_id]] <- gr
}

annotation <- list(cds_genes = cds_genes, cds_txs = cds_txs)

orfs_df <- read.table(file.path(script_dir, "test_orfs.tsv"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)

res <- orfquant_classify_orfs_py(
  orfs = orfs_df,
  annotation = annotation,
  gene_id_col = "gene_id",
  transcript_id_col = "gene_id",
  seqname_col = "seqnames",
  strand_col = "strand",
  exons_col = "exons",
  orf_id_col = "orf_id"
)

write.csv(res, file = file.path(script_dir, "r_results.csv"), row.names = FALSE)
