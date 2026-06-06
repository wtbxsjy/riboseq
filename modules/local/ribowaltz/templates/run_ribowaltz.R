#!/usr/bin/env Rscript

# riboWaltz P-site analysis for Ribo-seq QC
#
# Bioc 3.20+ compatible version. Uses pre-built cached annotation (RDS)
# to avoid rebuilding TxDb from GTF for each sample (~3 min → <1 sec).
# Annotation cache is built once from the GTF using txdbmaker::
# (GenomicFeatures::makeTxDbFromGFF is defunct in Bioc >= 3.20).

suppressPackageStartupMessages({
    library(riboWaltz)
    library(data.table)
    library(ggplot2)
    library(GenomicRanges)
    library(GenomicFeatures)
    library(GenomicAlignments)
    library(Rsamtools)
    library(GenomeInfoDb)
})

################################################
## Parse Parameters                           ##
################################################

prefix        <- '${prefix}'
bam_file      <- '${bam}'
gtf_file      <- '${gtf_file}'
read_lengths  <- ${read_lengths_r}

cat("=== riboWaltz Analysis ===\n")
cat("Sample:", prefix, "\n")
cat("BAM:", bam_file, "\n")
cat("GTF:", gtf_file, "\n")
cat("Read lengths:", paste(read_lengths, collapse = ", "), "\n")

################################################
## Create Output Directory                    ##
################################################

plot_dir <- paste0(prefix, "_ribowaltz_plots")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

################################################
## Helper: strip version from transcript IDs  ##
################################################

strip_version <- function(x) sub("\\.[0-9]+$", "", x)

################################################
## 1. Load Annotation                         ##
##    Prefer cached RDS, fall back to GTF     ##
################################################

cat("\n[1/7] Loading annotation...\n")

annotation <- tryCatch({
    # Try cached RDS first (built once with txdbmaker, loads in <1 sec).
    # Cache path is passed via ANNOTATION_CACHE env var from the module script.
    cache_file <- Sys.getenv("ANNOTATION_CACHE", unset = "")
    cache_loaded <- FALSE
    if (nchar(cache_file) > 0 && file.exists(cache_file)) {
        cat("Loading cached annotation:", cache_file, "\n")
        ann <- readRDS(cache_file)
        if (is.data.frame(ann) && nrow(ann) > 0 &&
            all(c("transcript", "l_tr", "l_utr5", "l_cds", "l_utr3") %in% names(ann))) {
            cat("Cached annotation loaded:", nrow(ann), "transcripts\n")
            cache_loaded <- TRUE
        } else {
            cat("Cached annotation invalid, rebuilding...\n")
        }
    }
    if (cache_loaded) {
        as.data.table(ann)
    } else {
        # Rebuild annotation from GTF using txdbmaker (Bioc 3.20+ compatible)
        cat("Building annotation from GTF using txdbmaker...\n")
        if (!requireNamespace("txdbmaker", quietly = TRUE)) {
            stop("txdbmaker is required but not installed")
        }
        txdb <- txdbmaker::makeTxDbFromGFF(gtf_file, format = "gtf",
                                            dataSource = "gencode", organism = NA)

        suppressWarnings({
            exon <- GenomicFeatures::exonsBy(txdb, by = "tx", use.names = TRUE)
            utr5 <- GenomicFeatures::fiveUTRsByTranscript(txdb, use.names = TRUE)
            cds  <- GenomicFeatures::cdsBy(txdb, by = "tx", use.names = TRUE)
            utr3 <- GenomicFeatures::threeUTRsByTranscript(txdb, use.names = TRUE)
        })

        exon_dt <- as.data.table(exon[unique(names(exon))])
        utr5_dt <- as.data.table(utr5[unique(names(utr5))])
        cds_dt  <- as.data.table(cds[unique(names(cds))])
        utr3_dt <- as.data.table(utr3[unique(names(utr3))])

        anno_df <- exon_dt[, list(l_tr = sum(width)), by = list(transcript = group_name)]
        l_utr5   <- utr5_dt[, list(l_utr5 = sum(width)), by = list(transcript = group_name)]
        l_cds    <- cds_dt[, list(l_cds = sum(width)), by = list(transcript = group_name)]
        l_utr3   <- utr3_dt[, list(l_utr3 = sum(width)), by = list(transcript = group_name)]

        merge_allx <- function(x, y) merge(x, y, all.x = TRUE)
        anno_df <- Reduce(merge_allx, list(anno_df, l_utr5, l_cds, l_utr3))
        anno_df[is.na(anno_df)] <- 0

        anno_df[, transcript := strip_version(transcript)]
        anno_df <- anno_df[!duplicated(transcript)]

        # Cache for future runs
        saveRDS(anno_df, cache_file, compress = TRUE)
        cat("Annotation built and cached:", nrow(anno_df), "transcripts\n")
        anno_df
    }
}, error = function(e) {
    cat("ERROR creating annotation:", conditionMessage(e), "\n")
    cat("Trying riboWaltz::create_annotation (may fail with Bioc 3.20+)...\n")
    ann <- create_annotation(gtfpath = gtf_file, dataSource = "gencode", organism = NA)
    ann[, transcript := strip_version(transcript)]
    ann[!duplicated(transcript)]
})

if (nrow(annotation) > 0) {
    cat("Sample transcript IDs:", paste(head(annotation$transcript, 5), collapse = ", "), "\n")
}

################################################
## 2. Load BAM Data                           ##
##    Custom loader with version stripping     ##
################################################

cat("\n[2/7] Loading BAM data with transcript version matching...\n")

reads_list <- tryCatch({
    suppressWarnings({
        ga <- GenomicAlignments::readGAlignments(bam_file, use.names = TRUE,
            param = Rsamtools::ScanBamParam(
                what = c("qname", "flag", "rname", "strand", "pos", "qwidth", "mapq"),
                flag = Rsamtools::scanBamFlag(
                    isUnmappedQuery = FALSE,
                    isSecondaryAlignment = FALSE,
                    isNotPassingQualityControls = FALSE
                )))
    })

    # Strip version numbers from transcript names.
    # renameSeqlevels() requires unique seqlevels — impossible after stripping
    # because multiple isoforms (AT1G01020.1, AT1G01020.2) collapse to the
    # same gene-level ID. seqlevels<- also fails because old seqlevels are
    # "in use" by reads.
    # Solution: convert to character early, filter, and build the data.table
    # directly — avoids all GenomicRanges seqlevels/seqnames constraints.
    stripped_names <- strip_version(as.character(seqnames(ga)))
    keep <- stripped_names %in% annotation$transcript
    ga <- ga[keep]
    stripped_names <- stripped_names[keep]

    if (length(ga) == 0) {
        cat("WARNING: No reads mapped to annotated transcripts\n")
        sample_name <- gsub("-", ".", prefix)
        setNames(list(data.table(
            transcript = character(), end5 = integer(), end3 = integer(),
            length = integer(), str = character(),
            cds_start = integer(), cds_stop = integer()
        )), sample_name)
    } else {
        dt <- data.table(
            transcript = stripped_names,
            end5       = ifelse(as.character(strand(ga)) == "+",
                                start(ga), end(ga)),
            end3       = ifelse(as.character(strand(ga)) == "+",
                                end(ga), start(ga)),
            length     = qwidth(ga)
        )
        dt[, str := as.character(strand(ga))]

        # Add cds_start/cds_stop from annotation (required by psite() and
        # length_filter periodicity mode).  Equivalent to what bamtolist does.
        ann_cols <- annotation[, .(transcript, l_utr5, l_cds)]
        dt <- dt[ann_cols, on = 'transcript',
                 c('cds_start', 'cds_stop') := list(i.l_utr5 + 1, i.l_utr5 + i.l_cds)]
        dt[cds_start == 1 & cds_stop == 0, cds_start := 0]

        sample_name <- gsub("-", ".", prefix)
        setNames(list(dt), sample_name)
    }
}, error = function(e) {
    cat("ERROR loading BAM:", conditionMessage(e), "\n")
    quit(status = 1)
})

total_reads <- sum(sapply(reads_list, nrow))
cat("BAM loaded:", total_reads, "reads across", length(reads_list), "sample(s)\n")

if (total_reads == 0) {
    cat("ERROR: No reads loaded from BAM. Check BAM contains transcriptome alignments.\n")
    cat("Annotation transcript IDs (first 5):",
        paste(head(annotation$transcript, 5), collapse = ", "), "\n")
    quit(status = 1)
}

################################################
## 3. Filter by Read Length                   ##
################################################

cat("\n[3/7] Filtering reads by length (", paste(read_lengths, collapse = ", "), ")...\n")

filtered_data <- tryCatch({
    lf_args <- list(data = reads_list, length_filter_mode = "custom")
    if ("length_range" %in% names(formals(riboWaltz::length_filter))) {
        lf_args$length_range <- read_lengths
    } else {
        lf_args$length_filter_vector <- read_lengths
    }
    do.call(length_filter, lf_args)
}, error = function(e) {
    cat("ERROR during length filtering:", conditionMessage(e), "\n")
    quit(status = 1)
})

cat("After filtering:", nrow(filtered_data[[1]]), "reads\n")

if (nrow(filtered_data[[1]]) == 0) {
    cat("WARNING: No reads remain after length filtering. Skipping P-site analysis.\n")
    dt_empty <- data.table(
        length = integer(), total_percentage = numeric(),
        start_percentage = numeric(), around_start = numeric(),
        offset_from_5 = numeric(), offset_from_3 = numeric(),
        corrected_offset_from_5 = numeric(), corrected_offset_from_3 = numeric(),
        sample = character()
    )
    fwrite(dt_empty, paste0(prefix, "_psite_offset.txt"), sep = "\t")
    fwrite(dt_empty, paste0(prefix, "_psite_offset.tsv"), sep = "\t")
    writeLines(
        c(paste0('"', '${task.process}', '":'),
          paste0('    ribowaltz: "', packageVersion("riboWaltz"), '"'),
          paste0('    r-data.table: "', packageVersion("data.table"), '"'),
          paste0('    r-ggplot2: "', packageVersion("ggplot2"), '"')),
        "versions.yml"
    )
    cat("\n=== riboWaltz analysis complete (no reads after filtering) ===\n")
    quit(status = 0)
}

################################################
## 4. Calculate P-site Offsets                ##
################################################

cat("\n[4/7] Calculating P-site offsets...\n")

psite_offset <- tryCatch({
    # NOTE: psite() requires txt_file to contain a directory separator;
    # passing bare filename causes dir.create("") error in riboWaltz.
    psite(
        data = filtered_data,
        flanking = 6,
        extremity = "auto",
        plot = TRUE,
        plot_dir = plot_dir,
        plot_format = "pdf",
        txt = TRUE,
        txt_file = file.path(getwd(), paste0(prefix, "_psite_offset.txt"))
    )
}, error = function(e) {
    cat("WARNING: P-site offset calculation failed:", conditionMessage(e), "\n")
    dt_placeholder <- data.table(
        length = integer(), total_percentage = numeric(),
        start_percentage = numeric(), around_start = numeric(),
        offset_from_5 = numeric(), offset_from_3 = numeric(),
        corrected_offset_from_5 = numeric(), corrected_offset_from_3 = numeric(),
        sample = character()
    )
    fwrite(dt_placeholder, paste0(prefix, "_psite_offset.txt"), sep = "\t")
    return(dt_placeholder)
})

cat("P-site offsets calculated\n")
fwrite(psite_offset, paste0(prefix, "_psite_offset.tsv"), sep = "\t")

if (nrow(psite_offset) > 0) {
    cat("Offset summary:\n")
    print_cols <- intersect(c("length", "corrected_offset_from_5", "corrected_offset_from_3"),
                            names(psite_offset))
    print(psite_offset[, ..print_cols])
}

################################################
## 5. Assign P-sites to Reads                 ##
################################################

cat("\n[5/7] Assigning P-site positions...\n")

psite_data <- tryCatch({
    psite_info(filtered_data, psite_offset)
}, error = function(e) {
    cat("WARNING: P-site assignment failed:", conditionMessage(e), "\n")
    return(filtered_data)
})

################################################
## 6. QC Analyses                             ##
################################################

cat("\n[6/7] Running QC analyses...\n")

# 6a. CDS Coverage
cat("  - CDS coverage...\n")
tryCatch({
    cds_cov <- cds_coverage(psite_data, annotation)
    fwrite(cds_cov, paste0(prefix, "_cds_coverage.tsv"), sep = "\t")
}, error = function(e) {
    cat("  WARNING: CDS coverage failed:", conditionMessage(e), "\n")
})

# 6b. Codon Usage (requires FASTA; skipped if unavailable)
cat("  - Codon usage (skipped - requires BSgenome/FASTA)...\n")

# 6c. Frame Distribution
cat("  - Frame distribution...\n")
tryCatch({
    sample_name <- names(psite_data)[1]
    frame_dist <- frame_psite(psite_data, annotation, sample = sample_name)
    fwrite(frame_dist, paste0(prefix, "_frame_distribution.tsv"), sep = "\t")
}, error = function(e) {
    cat("  WARNING: Frame distribution failed:", conditionMessage(e), "\n")
})

# 6d. Read Length Distribution
cat("  - Read length distribution...\n")
tryCatch({
    rl_plot <- riboWaltz::read_length_plot(data = filtered_data, sample = names(filtered_data)[1])
    ggsave(file.path(plot_dir, paste0(prefix, "_read_length_distribution.pdf")),
           plot = rl_plot, width = 8, height = 6)
}, error = function(e) {
    cat("  WARNING: Read length plot failed:", conditionMessage(e), "\n")
})

# 6e. Metagene Profiles
cat("  - Metagene profiles...\n")
tryCatch({
    sample_name <- names(psite_data)[1]
    if (exists("metaprofile_psite", where = "package:riboWaltz", inherits = FALSE)) {
        mp <- metaprofile_psite(psite_data, annotation, sample = sample_name)
        ggsave(file.path(plot_dir, paste0(prefix, "_metaprofile.pdf")),
               plot = mp, width = 10, height = 6)
    }
    if (exists("metaheatmap_psite", where = "package:riboWaltz", inherits = FALSE)) {
        mh <- metaheatmap_psite(psite_data, annotation, sample = sample_name)
        ggsave(file.path(plot_dir, paste0(prefix, "_metaheatmap.pdf")),
               plot = mh, width = 10, height = 8)
    }
}, error = function(e) {
    cat("  WARNING: Metagene profiles failed:", conditionMessage(e), "\n")
})

# 6f. P-site Region Distribution
# (riboWaltz >=2.0 adds region info via psite_info; percentage_regions removed)
cat("  - Region distribution (from psite_info)...\n")
tryCatch({
    if ("psite_region" %in% names(psite_data[[1]])) {
        region_dist <- psite_data[[1]][, .N, by = psite_region]
        fwrite(region_dist, paste0(prefix, "_region_distribution.tsv"), sep = "\t")
    }
}, error = function(e) {
    cat("  WARNING: Region distribution failed:", conditionMessage(e), "\n")
})

################################################
## Write Versions                             ##
################################################

cat("\nWriting version info...\n")
writeLines(
    c(
        paste0('"', '${task.process}', '":'),
        paste0('    ribowaltz: "', packageVersion("riboWaltz"), '"'),
        paste0('    r-data.table: "', packageVersion("data.table"), '"'),
        paste0('    r-ggplot2: "', packageVersion("ggplot2"), '"')
    ),
    "versions.yml"
)

cat("\n=== riboWaltz analysis complete ===\n")
