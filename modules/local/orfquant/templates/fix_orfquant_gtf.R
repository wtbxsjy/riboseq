# fix_orfquant_gtf.R — Rebuild ORFquant GTF with proper attributes from final_results
#
# ORFquant's run_ORFquant() produces ORFs_gen without metadata columns when
# ORFs_tx is empty (common even with mirai backend).  rtracklayer::export()
# then writes GTF attribute column as ".".  This script rebuilds the GTF
# using base R I/O, eliminating the rtracklayer dependency (which fails to
# load in the ORFquant mirai container due to host-R library-path leakage of
# an R-4.6.1 curl.so into R 4.4.3).

args <- commandArgs(trailingOnly = TRUE)
prefix <- args[1]

cat("[GTF fix] Loading ORFquant results...\n")
load(paste0(prefix, "_final_ORFquant_results"))

if (!exists("ORFquant_results") || length(ORFquant_results[["ORFs_gen"]]) == 0) {
    cat("[GTF fix] No ORFs found, skipping GTF rewrite\n")
    quit(save = "no", status = 0)
}

g <- ORFquant_results[["ORFs_gen"]]

# --- Parse protein FASTA for per-ORF metadata --------------------------------
fa_file <- paste0(prefix, "_Protein_sequences.fasta")
if (!file.exists(fa_file)) {
    cat("[GTF fix] WARNING: protein FASTA not found, exporting with coordinate-only attributes\n")
    fa_meta <- NULL
} else {
    fa_lines <- readLines(fa_file)
    fa_hdrs <- sub("^>", "", grep("^>", fa_lines, value = TRUE))
    tx_info <- strsplit(fa_hdrs, "\\|")
    fa_meta <- data.frame(
        orf_name     = sapply(tx_info, `[`, 1),
        gene_biotype = sapply(tx_info, `[`, 2),
        gene_id      = sapply(tx_info, `[`, 3),
        orf_type     = sapply(tx_info, `[`, 4),
        orf_category = ifelse(lengths(tx_info) >= 5, sapply(tx_info, `[`, 5), "NA"),
        stringsAsFactors = FALSE
    )
    rownames(fa_meta) <- fa_meta[["orf_name"]]
    cat(sprintf("[GTF fix] Parsed %d FASTA entries\n", nrow(fa_meta)))
}

# --- Identify ORF names ------------------------------------------------------
orf_names <- names(g)
if (is.null(orf_names)) {
    orf_names <- paste0("ORFquant_", seq_along(g))
}
unique_orfs <- unique(orf_names)
cat(sprintf("[GTF fix] %d CDS features across %d unique ORFs\n",
            length(g), length(unique_orfs)))

# --- Build per-row metadata vectors (plain R vectors, not GRanges mcols) -----
n <- length(g)

orf_id_vec   <- orf_names
gene_id_vec  <- rep("NA", n)
biotype_vec  <- rep("NA", n)
orf_type_vec <- rep("NA", n)
orf_cat_vec  <- rep("NA", n)

if (!is.null(fa_meta)) {
    matched <- orf_names %in% rownames(fa_meta)
    if (sum(matched) > 0) {
        gene_id_vec[matched]   <- fa_meta[orf_names[matched], "gene_id"]
        biotype_vec[matched]   <- fa_meta[orf_names[matched], "gene_biotype"]
        orf_type_vec[matched]  <- fa_meta[orf_names[matched], "orf_type"]
        orf_cat_vec[matched]   <- fa_meta[orf_names[matched], "orf_category"]
        cat(sprintf("[GTF fix] Metadata assigned to %d / %d features\n", sum(matched), n))
    }
}

# --- Write GTF with base R --------------------------------------------------
gtf_file <- paste0(prefix, "_Detected_ORFs.gtf")
cat(sprintf("[GTF fix] Writing %d features to %s\n", n, gtf_file))

seqname <- as.character(GenomeInfoDb::seqnames(g))
src_vec <- rep("ORFquant", n)
feat    <- rep("CDS", n)
starts  <- GenomicRanges::start(g)
ends    <- GenomicRanges::end(g)
score   <- rep(".", n)
strand  <- as.character(GenomicRanges::strand(g))
frame   <- rep(".", n)

# Build attribute strings in chunks to avoid memory blow-up
CHUNK <- 5000
n_chunks <- ceiling(n / CHUNK)

writeLines("##gtf-version 2", gtf_file)
con <- file(gtf_file, "a")

for (k in seq_len(n_chunks)) {
    idx_start <- (k - 1) * CHUNK + 1
    idx_end   <- min(k * CHUNK, n)
    idx <- idx_start:idx_end

    attrs <- sprintf(
        'gene_id "%s"; transcript_id "%s"; ORF_id "%s"; gene_biotype "%s"; orf_type "%s"; orf_category "%s"; source "ORFquant";',
        gene_id_vec[idx],
        orf_id_vec[idx],
        orf_id_vec[idx],
        biotype_vec[idx],
        orf_type_vec[idx],
        orf_cat_vec[idx]
    )

    lines <- paste(
        seqname[idx], src_vec[idx], feat[idx],
        starts[idx], ends[idx], score[idx], strand[idx], frame[idx],
        attrs,
        sep = "\t"
    )
    writeLines(lines, con)

    if (k %% 10 == 0) {
        cat(sprintf("[GTF fix]   wrote %d / %d features\n", idx_end, n))
    }
}
close(con)
cat("[GTF fix] Done\n")
