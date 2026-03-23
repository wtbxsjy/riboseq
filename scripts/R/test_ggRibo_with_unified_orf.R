#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggRibo)
  library(txdbmaker)
  library(GenomicFeatures)
  library(AnnotationDbi)
  library(data.table)
  library(ggplot2)
})

usage <- function() {
  cat(
    "Usage:\n",
    "  Rscript scripts/R/test_ggRibo_with_unified_orf.R \\\n",
    "    --metadata test_data/orf_unification_mouse_Mucosal_Immunity/unified_orfs.metadata.tsv \\\n",
    "    --gtf test_data/orf_unification_mouse_Mucosal_Immunity/unified_orfs.gtf \\\n",
    "    --riboseqc-dir test_data/riboseqc \\\n",
    "    --orf-id ORF_8_ENSMUSG00000033793.13 \\\n",
    "    --outdir test_results/ggRibo_pkg_demo\n",
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
    orf_id = NULL,
    outdir = NULL
  )

  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    value <- args[[i + 1L]]
    switch(
      key,
      "--metadata" = out$metadata <- value,
      "--gtf" = out$gtf <- value,
      "--riboseqc-dir" = out$riboseqc_dir <- value,
      "--orf-id" = out$orf_id <- value,
      "--outdir" = out$outdir <- value,
      stop(sprintf("Unknown argument: %s", key), call. = FALSE)
    )
    i <- i + 2L
  }

  missing <- names(out)[vapply(out, is.null, logical(1))]
  if (length(missing) > 0L) {
    stop(sprintf("Missing required arguments: %s", paste(missing, collapse = ", ")), call. = FALSE)
  }

  out
}

clean_attrs <- function(x) {
  x <- sub('^"', "", x)
  x <- sub('"$', "", x)
  gsub('""', '"', x, fixed = TRUE)
}

extract_attr <- function(attrs, key) {
  match <- regexec(sprintf('%s "([^"]+)"', key), attrs)
  hit <- regmatches(attrs, match)[[1]]
  if (length(hit) < 2L) {
    return(NA_character_)
  }
  hit[2]
}

augment_gtf_for_txdb <- function(gtf_path, orf_id, outdir) {
  gtf <- fread(gtf_path, sep = "\t", header = FALSE, quote = "")
  setnames(gtf, paste0("V", 1:9))
  gtf[, V9 := clean_attrs(V9)]
  subset_gtf <- gtf[grepl(sprintf('orf_id "%s"', orf_id), V9, fixed = TRUE)]
  if (nrow(subset_gtf) == 0L) {
    stop(sprintf("ORF %s not found in %s", orf_id, gtf_path), call. = FALSE)
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
  augmented_path <- file.path(outdir, sprintf("%s.ggRibo_input.gtf", orf_id))
  fwrite(augmented, augmented_path, sep = "\t", col.names = FALSE, quote = FALSE)

  list(path = augmented_path, gene_id = gene_id, tx_id = tx_id)
}

build_range_info <- function(annotation) {
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

  Range_info$new(
    exonsByTx = exonsByTx,
    txByGene = txByGene,
    cdsByTx = cdsByTx,
    fiveUTR = fiveUTR,
    threeUTR = threeUTR,
    tx_to_gene = tx_to_gene
  )
}

main <- function() {
  opt <- parse_args(commandArgs(trailingOnly = TRUE))
  dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

  meta <- fread(
    opt$metadata,
    sep = "\t",
    select = c("orf_id", "gene_id", "transcript_id", "samples", "tools", "unique_psites", "pN")
  )
  row <- meta[orf_id == opt$orf_id][1]
  if (nrow(row) == 0L) {
    stop(sprintf("ORF %s not found in metadata", opt$orf_id), call. = FALSE)
  }

  annotation_info <- augment_gtf_for_txdb(opt$gtf, opt$orf_id, opt$outdir)
  samples <- strsplit(row$samples, ",", fixed = TRUE)[[1]]
  ribo_files <- file.path(opt$riboseqc_dir, paste0(samples, "_ggribo.tsv"))
  missing_files <- ribo_files[!file.exists(ribo_files)]
  if (length(missing_files) > 0L) {
    stop(sprintf("Missing ggRibo tabular files: %s", paste(missing_files, collapse = ", ")), call. = FALSE)
  }

  inputs_full <- create_seq_input(
    ribo_files = as.list(ribo_files),
    sample_names = samples,
    include_rna = FALSE
  )
  GRangeInfo <- build_range_info(annotation_info$path)

  p <- ggRibo_tx(
    gene_id = annotation_info$gene_id,
    tx_id = annotation_info$tx_id,
    NAME = sprintf("%s | tools=%s | unique_psites=%s | pN=%s", row$orf_id, row$tools, row$unique_psites, row$pN),
    RNAseq = NULL,
    Riboseq = inputs_full$Riboseq,
    SampleNames = samples,
    GRangeInfo = GRangeInfo,
    Y_scale = "each",
    plot_ORF_ranges = FALSE,
    show_seq = FALSE,
    title_font_size = 10,
    sample_label_font_size = 4
  )

  out_png <- file.path(opt$outdir, sprintf("%s.ggRibo_tx.png", opt$orf_id))
  ggsave(out_png, p, width = 13, height = max(6, 2 + 2 * length(samples)), dpi = 150, limitsize = FALSE)

  message(sprintf("Wrote %s", out_png))
  message(sprintf("Augmented GTF: %s", annotation_info$path))
}

main()
