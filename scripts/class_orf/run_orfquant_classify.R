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
              help="Output classification TSV file", metavar="file"),
  make_option(c("-m", "--metadata"), type="character", default=NULL,
              help="Unified ORF metadata TSV (for FA/per-sample outputs)", metavar="file"),
  make_option(c("-p", "--output_prefix"), type="character", default="orfquant_results",
              help="Prefix for extra output files (BED/GTF/FA/logs/out) [default: orfquant_results]")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input) || is.null(opt$annotation) || is.null(opt$output)) {
  print_help(opt_parser)
  stop("Missing required arguments", call.=FALSE)
}

message("Loading ORFs from: ", opt$input)
message("Loading Annotation from: ", opt$annotation)
message("Running classification...")

orfs_gr <- rtracklayer::import(opt$input)
orfs_gr <- orfs_gr[orfs_gr$type == "CDS"]
orfs_grl <- split(orfs_gr, orfs_gr$orf_id)

results <- orfquant_classify_orfs(
  orfs = orfs_grl,
  annotation = opt$annotation
)

message("Writing results to: ", opt$output)
write.table(results, file = opt$output, sep = "\t", quote = FALSE, row.names = FALSE)

# ─────────────────────────────────────────────────────────────
# Extra outputs: BED12, GTF, nucleotide FA, protein FA, logs, per-sample out
# ─────────────────────────────────────────────────────────────
prefix <- opt$output_prefix

# -- BED12 from the input GRangesList -------------------------
message("Writing BED12: ", prefix, ".orfs.bed")
bed_lines <- character(0)
for (nm in names(orfs_grl)) {
  gr <- orfs_grl[[nm]]
  if (length(gr) == 0) next
  chrom   <- as.character(GenomicRanges::seqnames(gr)[1])
  strand  <- as.character(GenomicRanges::strand(gr)[1])
  starts  <- IRanges::start(gr)
  ends    <- IRanges::end(gr)
  chrom_start <- min(starts) - 1L   # 0-based
  chrom_end   <- max(ends)           # 0-based exclusive = 1-based inclusive
  sizes  <- ends - starts + 1L
  bstarts <- starts - 1L - chrom_start
  gene_id <- if (!is.null(gr$gene_id)) gr$gene_id[1] else "NA"
  bed_lines <- c(bed_lines,
    paste(chrom, chrom_start, chrom_end, nm, 0, strand,
          chrom_start, chrom_end, "0,0,0",
          length(gr),
          paste(sizes,  collapse=","),
          paste(bstarts, collapse=","),
          sep="\t"))
}
writeLines(bed_lines, con=paste0(prefix, ".orfs.bed"))

# -- GTF from the input GRangesList ---------------------------
message("Writing GTF: ", prefix, ".orfs.gtf")
ann_gr <- unlist(orfs_grl, use.names=FALSE)
# Ensure orf_id is preserved as a metadata column
if (is.null(ann_gr$orf_id)) {
  ann_gr$orf_id <- rep(names(orfs_grl), lengths(orfs_grl))
}
# Annotate with classification results
orf_class_map <- setNames(results$ORF_type_py, results$orf_id)
ann_gr$orf_category <- orf_class_map[ann_gr$orf_id]
rtracklayer::export(ann_gr, con=paste0(prefix, ".orfs.gtf"), format="GTF")

# -- Nucleotide and protein FASTA (from metadata) -------------
if (!is.null(opt$metadata) && file.exists(opt$metadata)) {
  message("Writing FASTA from metadata: ", opt$metadata)
  meta <- tryCatch(
    read.table(opt$metadata, sep="\t", header=TRUE, stringsAsFactors=FALSE,
               quote="", comment.char=""),
    error = function(e) { message("Warning: could not read metadata: ", e$message); NULL }
  )
  if (!is.null(meta) && "sequence" %in% colnames(meta) && "orf_id" %in% colnames(meta)) {
    codon_table <- c(
      TTT="F",TTC="F",TTA="L",TTG="L",CTT="L",CTC="L",CTA="L",CTG="L",
      ATT="I",ATC="I",ATA="I",ATG="M",GTT="V",GTC="V",GTA="V",GTG="V",
      TCT="S",TCC="S",TCA="S",TCG="S",CCT="P",CCC="P",CCA="P",CCG="P",
      ACT="T",ACC="T",ACA="T",ACG="T",GCT="A",GCC="A",GCA="A",GCG="A",
      TAT="Y",TAC="Y",TAA="*",TAG="*",CAT="H",CAC="H",CAA="Q",CAG="Q",
      AAT="N",AAC="N",AAA="K",AAG="K",GAT="D",GAC="D",GAA="E",GAG="E",
      TGT="C",TGC="C",TGA="*",TGG="W",CGT="R",CGC="R",CGA="R",CGG="R",
      AGT="S",AGC="S",AGA="R",AGG="R",GGT="G",GGC="G",GGA="G",GGG="G"
    )
    translate_nt <- function(nt) {
      nt <- toupper(chartr("U", "T", nt))
      n  <- nchar(nt)
      codons <- substring(nt, seq(1, n-2, 3), seq(3, n, 3))
      paste(ifelse(codons %in% names(codon_table), codon_table[codons], "X"), collapse="")
    }
    fa_nt  <- file(paste0(prefix, ".orfs.fa"),     "w")
    fa_pep <- file(paste0(prefix, ".orfs.pep.fa"), "w")
    for (i in seq_len(nrow(meta))) {
      oid <- meta$orf_id[i]
      seq <- meta$sequence[i]
      if (is.na(seq) || seq == "") next
      chrom_col <- if ("chrom" %in% colnames(meta)) meta$chrom[i] else ""
      writeLines(c(paste0(">", oid, " ", chrom_col), seq), fa_nt)
      writeLines(c(paste0(">", oid, " ", chrom_col), translate_nt(seq)), fa_pep)
    }
    close(fa_nt); close(fa_pep)
  }

  # -- Per-sample stats .orfs.out --------------------------------
  if (!is.null(meta) && "samples" %in% colnames(meta)) {
    message("Writing per-sample stats: ", prefix, ".orfs.out")
    all_samples <- sort(unique(unlist(strsplit(
      meta$samples[!is.na(meta$samples) & meta$samples != ""], ","))))
    all_samples <- trimws(all_samples)
    all_samples <- all_samples[all_samples != ""]

    # Merge classification into metadata
    class_cols <- results[, c("orf_id", "ORF_type_py",
                               if ("ORF_category_Tx" %in% colnames(results)) "ORF_category_Tx" else NULL,
                               if ("ORF_category_Tx_compatible" %in% colnames(results)) "ORF_category_Tx_compatible" else NULL),
                          drop=FALSE]
    merged <- merge(meta, class_cols, by="orf_id", all.x=TRUE)

    # Build per-sample 0/1 columns
    sample_mat <- do.call(cbind, lapply(all_samples, function(s) {
      as.integer(sapply(strsplit(as.character(merged$samples), ","), function(v) s %in% trimws(v)))
    }))
    colnames(sample_mat) <- all_samples

    n_samples_vec <- rowSums(sample_mat)

    base_cols <- c("orf_id", "chrom", "start", "end", "strand",
                   "gene_id", "transcript_id", "length_aa", "tools")
    class_col_names <- c("ORF_type_py",
                         if ("ORF_category_Tx" %in% colnames(merged)) "ORF_category_Tx" else NULL,
                         if ("ORF_category_Tx_compatible" %in% colnames(merged)) "ORF_category_Tx_compatible" else NULL)
    avail_cols <- base_cols[base_cols %in% colnames(merged)]
    avail_class <- class_col_names[class_col_names %in% colnames(merged)]

    out_df <- cbind(
      merged[, avail_cols, drop=FALSE],
      n_samples = n_samples_vec,
      merged[, avail_class, drop=FALSE],
      sample_mat
    )
    write.table(out_df, file=paste0(prefix, ".orfs.out"),
                sep="\t", quote=FALSE, row.names=FALSE)
  }
}

# -- Logs -------------------------------------------------------
message("Writing logs: ", prefix, ".logs")
log_lines <- c(
  paste0("#input gtf: ",    opt$input),
  paste0("#annotation gtf: ", opt$annotation),
  paste0("#metadata: ",     if (!is.null(opt$metadata)) opt$metadata else "NA"),
  paste0("#total_orfs: ",   length(orfs_grl)),
  paste0("#output_prefix: ", prefix)
)
writeLines(log_lines, con=paste0(prefix, ".logs"))

message("Done.")

