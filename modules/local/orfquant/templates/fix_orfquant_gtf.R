suppressMessages({
    library(ORFquant)
    library(rtracklayer)
})

args <- commandArgs(trailingOnly = TRUE)
prefix <- args[1]

cat("[GTF fix] Loading ORFquant results...\n")
load(paste0(prefix, "_final_ORFquant_results"))

if (!exists("ORFquant_results") || length(ORFquant_results[["ORFs_gen"]]) == 0) {
    cat("[GTF fix] No ORFs found, skipping GTF rewrite\n")
    quit(save = "no", status = 0)
}

g <- ORFquant_results[["ORFs_gen"]]

# Parse protein FASTA for per-ORF metadata
fa_file <- paste0(prefix, "_Protein_sequences.fasta")
if (!file.exists(fa_file)) {
    cat("[GTF fix] WARNING: protein FASTA not found, exporting with coordinate-only attributes\n")
    fa_meta <- NULL
} else {
    fa_lines <- readLines(fa_file)
    fa_hdrs <- sub("^>", "", grep("^>", fa_lines, value = TRUE))
    tx_info <- strsplit(fa_hdrs, "\\|")
    fa_meta <- data.frame(
        orf_name  = sapply(tx_info, `[`, 1),
        gene_biotype = sapply(tx_info, `[`, 2),
        gene_id   = sapply(tx_info, `[`, 3),
        orf_type  = sapply(tx_info, `[`, 4),
        orf_category = ifelse(lengths(tx_info) >= 5, sapply(tx_info, `[`, 5), "NA"),
        stringsAsFactors = FALSE
    )
    rownames(fa_meta) <- fa_meta[["orf_name"]]
    cat(sprintf("[GTF fix] Parsed %d FASTA entries\n", nrow(fa_meta)))
}

orf_names <- names(g)
if (is.null(orf_names)) {
    orf_names <- paste0("ORFquant_", seq_along(g))
    names(g) <- orf_names
}
unique_orfs <- unique(orf_names)
cat(sprintf("[GTF fix] %d CDS features across %d unique ORFs\n", length(g), length(unique_orfs)))

n <- length(g)
g[["type"]] <- "CDS"
mcols(g)[["ORF_id"]]    <- orf_names
mcols(g)[["gene_id"]]   <- rep("NA", n)
mcols(g)[["gene_biotype"]] <- rep("NA", n)
mcols(g)[["orf_type"]]  <- rep("NA", n)
mcols(g)[["orf_category"]] <- rep("NA", n)

if (!is.null(fa_meta)) {
    matched <- orf_names %in% rownames(fa_meta)
    if (sum(matched) > 0) {
        mcols(g)[["gene_id"]][matched]      <- fa_meta[orf_names[matched], "gene_id"]
        mcols(g)[["gene_biotype"]][matched] <- fa_meta[orf_names[matched], "gene_biotype"]
        mcols(g)[["orf_type"]][matched]     <- fa_meta[orf_names[matched], "orf_type"]
        mcols(g)[["orf_category"]][matched] <- fa_meta[orf_names[matched], "orf_category"]
        cat(sprintf("[GTF fix] Metadata assigned to %d / %d features\n", sum(matched), n))
    }
}

mcols(g)[["source"]] <- "ORFquant"

gtf_file <- paste0(prefix, "_Detected_ORFs.gtf")
cat(sprintf("[GTF fix] Writing %d features to %s\n", length(g), gtf_file))
rtracklayer::export(g, gtf_file, format = "gtf")
cat("[GTF fix] Done\n")
