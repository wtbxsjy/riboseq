#!/usr/bin/env Rscript

# riboWaltz P-site analysis for Ribo-seq QC
#
# Bioc 3.20 compatible version using txdbmaker
# Handles transcript ID version mismatch between GTF and BAM

suppressPackageStartupMessages({
    library(riboWaltz)
    library(data.table)
    library(ggplot2)
})

################################################
## Parse Parameters                           ##
################################################

prefix        <- '${prefix}'
bam_file      <- '${bam}'
gtf_file      <- '${gtf}'
read_lengths  <- ${read_lengths_r}

cat("=== riboWaltz Analysis ===\\n")
cat("Sample:", prefix, "\\n")
cat("BAM:", bam_file, "\\n")
cat("GTF:", gtf_file, "\\n")
cat("Read lengths:", paste(read_lengths, collapse = ", "), "\\n")

################################################
## Create Output Directory                    ##
################################################

plot_dir <- paste0(prefix, "_ribowaltz_plots")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

################################################
## Helper: strip version from transcript IDs  ##
################################################

strip_version <- function(x) sub("\\\\.[0-9]+$", "", x)

################################################
## 1. Create Annotation from GTF              ##
##    Bioc 3.20 compatible via txdbmaker      ##
################################################

cat("\\n[1/7] Creating annotation from GTF...\\n")

annotation <- tryCatch({
    # Prefer txdbmaker for Bioc 3.20+, fall back to GenomicFeatures
    if (requireNamespace("txdbmaker", quietly = TRUE)) {
        txdb <- txdbmaker::makeTxDbFromGFF(gtf_file, format = "gtf",
                                            dataSource = "gencode", organism = NA)
        cat("Using txdbmaker::makeTxDbFromGFF\\n")
    } else {
        txdb <- GenomicFeatures::makeTxDbFromGFF(gtf_file, format = "gtf",
                                                  dataSource = "gencode", organism = NA)
        cat("Using GenomicFeatures::makeTxDbFromGFF\\n")
    }

    suppressWarnings({
        tx_gr <- GenomicFeatures::transcripts(txdb, columns = c("tx_name", "gene_id", "tx_type"))
    })
    tx_df <- as.data.frame(tx_gr)

    # Extract CDS per transcript
    cds_by_tx <- GenomicFeatures::cdsBy(txdb, by = "tx", use.names = TRUE)

    # Use tx_name if available (matches BAM transcript IDs), otherwise tx_id
    tx_keys <- if (!is.null(tx_df$tx_name)) tx_df$tx_name else as.character(tx_df$group_name)

    tx_width <- tx_df$width
    cds_width <- integer(length(tx_keys))
    cds_start <- integer(length(tx_keys))
    cds_end   <- integer(length(tx_keys))
    has_cds   <- tx_keys %in% names(cds_by_tx)

    for (i in seq_along(tx_keys)) {
        if (has_cds[i]) {
            cds_gr <- cds_by_tx[[tx_keys[i]]]
            cds_width[i] <- sum(width(cds_gr))
            cds_start[i] <- min(start(cds_gr))
            cds_end[i]   <- max(end(cds_gr))
        }
    }

    annotation <- data.table(
        transcript = strip_version(tx_keys),
        l_tr       = tx_width,
        l_utr5     = ifelse(has_cds, cds_start - 1L, tx_width),
        l_cds      = cds_width,
        l_utr3     = ifelse(has_cds, pmax(0L, tx_width - cds_end), 0L),
        gene       = if (!is.null(tx_df$gene_id)) tx_df$gene_id else NA_character_,
        from       = 1L,
        to         = tx_width,
        cds_from   = ifelse(has_cds, as.integer(cds_start), NA_integer_),
        cds_to     = ifelse(has_cds, as.integer(cds_end), NA_integer_)
    )

    # Remove zero-length transcripts and duplicates
    annotation <- annotation[l_tr > 0]
    annotation <- annotation[!duplicated(transcript)]

    annotation
}, error = function(e) {
    cat("ERROR creating annotation:", conditionMessage(e), "\\n")
    cat("Falling back to riboWaltz::create_annotation...\\n")
    tryCatch({
        ann <- create_annotation(gtfpath = gtf_file, dataSource = "gencode", organism = NA)
        ann[, transcript := strip_version(transcript)]
        ann[!duplicated(transcript)]
    }, error = function(e2) {
        cat("FATAL: Cannot create annotation:", conditionMessage(e2), "\\n")
        quit(status = 1)
    })
})

cat("Annotation created:", nrow(annotation), "transcripts\\n")
if (nrow(annotation) > 0) {
    cat("Sample transcript IDs:", paste(head(annotation$transcript, 5), collapse = ", "), "\\n")
}

################################################
## 2. Load BAM Data                           ##
##    Custom loader with version stripping     ##
################################################

cat("\\n[2/7] Loading BAM data with transcript version matching...\\n")

load_bam_single <- function(bam_path, annotation) {
    # Read BAM alignments
    suppressWarnings({
        ga <- GenomicAlignments::readGAlignments(bam_path, use.names = TRUE,
                                                  param = Rsamtools::ScanBamParam(
                                                      what = c("qname", "flag", "rname", "strand",
                                                                "pos", "qwidth", "mapq"),
                                                      flag = Rsamtools::scanBamFlag(
                                                          isUnmappedQuery = FALSE,
                                                          isSecondaryAlignment = FALSE,
                                                          isNotPassingQualityControls = FALSE
                                                      )
                                                  ))
    })

    # Strip version numbers from transcript names in BAM reads
    new_seqlevels <- strip_version(seqlevels(ga))
    names(new_seqlevels) <- seqlevels(ga)
    ga <- GenomeInfoDb::renameSeqlevels(ga, new_seqlevels)

    # Filter to transcripts in annotation
    ga <- ga[seqnames(ga) %in% annotation$transcript]

    if (length(ga) == 0) {
        cat("WARNING: No reads mapped to annotated transcripts\\n")
        return(data.table(
            transcript = character(), end5 = integer(), end3 = integer(),
            length = integer(), str = character()
        ))
    }

    # Create data.table in riboWaltz-compatible format
    dt <- data.table(
        transcript = as.character(seqnames(ga)),
        end5       = ifelse(as.character(strand(ga)) == "+",
                            start(ga), end(ga)),
        end3       = ifelse(as.character(strand(ga)) == "+",
                            end(ga), start(ga)),
        length     = qwidth(ga)
    )
    dt[, str := as.character(strand(ga))]

    return(dt)
}

reads_list <- tryCatch({
    dt <- load_bam_single(bam_file, annotation)
    sample_name <- gsub("-", ".", prefix)
    setnames(dt, "transcript", "transcript")
    set(list(dt), names = sample_name)
}, error = function(e) {
    cat("ERROR loading BAM:", conditionMessage(e), "\\n")
    cat("Falling back to riboWaltz::bamtolist...\\n")
    tryCatch({
        bam_dir <- dirname(bam_file)
        tmp_dir <- file.path(tempdir(), "rw_bam")
        dir.create(tmp_dir, showWarnings = FALSE, recursive = TRUE)
        base_bam <- basename(bam_file)
        file.symlink(bam_file, file.path(tmp_dir, base_bam))
        bai_file <- paste0(bam_file, ".bai")
        if (file.exists(bai_file)) {
            file.symlink(bai_file, file.path(tmp_dir, paste0(base_bam, ".bai")))
        }
        rl <- bamtolist(bamfolder = tmp_dir, annotation = annotation)
        unlink(tmp_dir, recursive = TRUE)
        for (i in seq_along(rl)) {
            if (nrow(rl[[i]]) > 0) {
                rl[[i]][, transcript := strip_version(transcript)]
            }
        }
        rl
    }, error = function(e2) {
        cat("FATAL loading BAM:", conditionMessage(e2), "\\n")
        quit(status = 1)
    })
})

total_reads <- sum(sapply(reads_list, nrow))
cat("BAM loaded:", total_reads, "reads across", length(reads_list), "sample(s)\\n")

if (total_reads == 0) {
    cat("ERROR: No reads loaded from BAM. Check that BAM contains transcriptome alignments.\\n")
    cat("Annotation transcript IDs (first 5):", paste(head(annotation$transcript, 5), collapse = ", "), "\\n")
    quit(status = 1)
}

################################################
## 3. Filter by Read Length                   ##
################################################

cat("\\n[3/7] Filtering reads by length (", paste(read_lengths, collapse = ", "), ")...\\n")

filtered_data <- tryCatch({
    # riboWaltz >=2.0 uses 'length_range', older versions use 'length_filter_vector'
    lf_args <- list(data = reads_list, length_filter_mode = "custom")
    if ("length_range" %in% names(formals(riboWaltz::length_filter))) {
        lf_args$length_range <- read_lengths
    } else {
        lf_args$length_filter_vector <- read_lengths
    }
    do.call(length_filter, lf_args)
}, error = function(e) {
    cat("ERROR during length filtering:", conditionMessage(e), "\\n")
    quit(status = 1)
})

cat("After filtering:", nrow(filtered_data[[1]]), "reads\\n")

if (nrow(filtered_data[[1]]) == 0) {
    cat("WARNING: No reads remain after length filtering. Skipping P-site analysis.\\n")
    cat("Creating minimal output files...\\n")
    dt_empty <- data.table(
        length = integer(), total_percentage = numeric(),
        start_percentage = numeric(), around_start = numeric(),
        offset_from_5 = numeric(), offset_from_3 = numeric(),
        corrected_offset_from_5 = numeric(), corrected_offset_from_3 = numeric(),
        sample = character()
    )
    fwrite(dt_empty, paste0(prefix, "_psite_offset.txt"), sep = "\\t")
    fwrite(dt_empty, paste0(prefix, "_psite_offset.tsv"), sep = "\\t")

    # Write version info and exit
    writeLines(
        c(
            paste0('"', '${task.process}', '":'),
            paste0('    ribowaltz: "', packageVersion("riboWaltz"), '"'),
            paste0('    r-data.table: "', packageVersion("data.table"), '"'),
            paste0('    r-ggplot2: "', packageVersion("ggplot2"), '"')
        ),
        "versions.yml"
    )
    cat("\\n=== riboWaltz analysis complete (no reads after filtering) ===\\n")
    quit(status = 0)
}

################################################
## 4. Calculate P-site Offsets                ##
################################################

cat("\\n[4/7] Calculating P-site offsets...\\n")

psite_offset <- tryCatch({
    psite(
        data = filtered_data,
        flanking = 6,
        extremity = "auto",
        plot = TRUE,
        plot_dir = plot_dir,
        plot_format = "pdf",
        txt = TRUE,
        txt_file = paste0(prefix, "_psite_offset.txt")
    )
}, error = function(e) {
    cat("WARNING: P-site offset calculation failed:", conditionMessage(e), "\\n")
    # Create minimal placeholder
    dt_placeholder <- data.table(
        length = integer(), total_percentage = numeric(),
        start_percentage = numeric(), around_start = numeric(),
        offset_from_5 = numeric(), offset_from_3 = numeric(),
        corrected_offset_from_5 = numeric(), corrected_offset_from_3 = numeric(),
        sample = character()
    )
    fwrite(dt_placeholder, paste0(prefix, "_psite_offset.txt"), sep = "\\t")
    return(dt_placeholder)
})

cat("P-site offsets calculated\\n")
fwrite(psite_offset, paste0(prefix, "_psite_offset.tsv"), sep = "\\t")

# Show offset summary
if (nrow(psite_offset) > 0) {
    cat("Offset summary:\\n")
    print_cols <- intersect(c("length", "corrected_offset_from_5", "corrected_offset_from_3"),
                            names(psite_offset))
    print(psite_offset[, ..print_cols])
}

################################################
## 5. Assign P-sites to Reads                 ##
################################################

cat("\\n[5/7] Assigning P-site positions...\\n")

psite_data <- tryCatch({
    psite_info(filtered_data, psite_offset)
}, error = function(e) {
    cat("WARNING: P-site assignment failed:", conditionMessage(e), "\\n")
    return(filtered_data)
})

################################################
## 6. QC Analyses                             ##
################################################

cat("\\n[6/7] Running QC analyses...\\n")

# 6a. CDS Coverage
cat("  - CDS coverage...\\n")
tryCatch({
    cds_cov <- cds_coverage(psite_data, annotation)
    fwrite(cds_cov, paste0(prefix, "_cds_coverage.tsv"), sep = "\\t")
}, error = function(e) {
    cat("  WARNING: CDS coverage failed:", conditionMessage(e), "\\n")
})

# 6b. Codon Usage
cat("  - Codon usage...\\n")
tryCatch({
    codon_usage <- codon_usage_psite(psite_data, annotation)
    if (nrow(codon_usage) > 0) {
        fwrite(codon_usage, paste0(prefix, "_codon_usage.tsv"), sep = "\\t")
    }
}, error = function(e) {
    cat("  WARNING: Codon usage failed:", conditionMessage(e), "\\n")
})

# 6c. Frame Distribution
cat("  - Frame distribution...\\n")
tryCatch({
    frame_dist <- frame_psite(psite_data, annotation)
    fwrite(frame_dist, paste0(prefix, "_frame_distribution.tsv"), sep = "\\t")
}, error = function(e) {
    cat("  WARNING: Frame distribution failed:", conditionMessage(e), "\\n")
})

# 6d. Read Length Distribution
cat("  - Read length distribution...\\n")
tryCatch({
    rl_plot <- read_length_plot(data = filtered_data, sample = names(filtered_data)[1])
    ggsave(file.path(plot_dir, paste0(prefix, "_read_length_distribution.pdf")),
           plot = rl_plot, width = 8, height = 6)
}, error = function(e) {
    cat("  WARNING: Read length plot failed:", conditionMessage(e), "\\n")
})

# 6e. Metagene Profiles
cat("  - Metagene profiles...\\n")
tryCatch({
    if (exists("metaplots")) {
        metaplots_output <- metaplots(psite_data, annotation, sample = names(psite_data)[1])
        ggsave(file.path(plot_dir, paste0(prefix, "_metaplots.pdf")),
               plot = metaplots_output, width = 10, height = 8)
    } else {
        # Try metaprofile + metaheatmap separately (riboWaltz >=2.0)
        if (exists("metaprofile_psite")) {
            mp <- metaprofile_psite(psite_data, annotation, sample = names(psite_data)[1])
            ggsave(file.path(plot_dir, paste0(prefix, "_metaprofile.pdf")),
                   plot = mp, width = 10, height = 6)
        }
        if (exists("metaheatmap_psite")) {
            mh <- metaheatmap_psite(psite_data, annotation, sample = names(psite_data)[1])
            ggsave(file.path(plot_dir, paste0(prefix, "_metaheatmap.pdf")),
                   plot = mh, width = 10, height = 8)
        }
    }
}, error = function(e) {
    cat("  WARNING: Metagene profiles failed:", conditionMessage(e), "\\n")
})

# 6f. P-site Region Distribution
cat("  - P-site region distribution...\\n")
tryCatch({
    if (exists("percentage_regions")) {
        pct_regions <- percentage_regions(psite_data, annotation)
        fwrite(pct_regions, paste0(prefix, "_region_distribution.tsv"), sep = "\\t")
    }
}, error = function(e) {
    cat("  WARNING: Region distribution failed:", conditionMessage(e), "\\n")
})

################################################
## Write Versions                             ##
################################################

cat("\\nWriting version info...\\n")
writeLines(
    c(
        paste0('"', '${task.process}', '":'),
        paste0('    ribowaltz: "', packageVersion("riboWaltz"), '"'),
        paste0('    r-data.table: "', packageVersion("data.table"), '"'),
        paste0('    r-ggplot2: "', packageVersion("ggplot2"), '"')
    ),
    "versions.yml"
)

cat("\\n=== riboWaltz analysis complete ===\\n")
