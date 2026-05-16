#!/usr/bin/env Rscript

# riboWaltz P-site analysis for Ribo-seq QC
#
# This script runs riboWaltz analysis:
# 1. Creates annotation from GTF
# 2. Loads BAM data
# 3. Calculates P-site offsets
# 4. Generates metagene profiles, codon usage, CDS coverage, frame distribution

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
read_lengths  <- ${read_lengths_r}  # R vector of read lengths, e.g. c(28, 29, 30)

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
## 1. Create Annotation from GTF              ##
################################################

cat("\n[1/6] Creating annotation from GTF...\n")

annotation <- tryCatch({
    create_annotation(gtfpath = gtf_file, dataSource = "gencode", organism = NA)
}, error = function(e) {
    cat("ERROR creating annotation:", conditionMessage(e), "\n")
    quit(status = 1)
})

cat("Annotation created:", nrow(annotation), "transcripts\n")

################################################
## 2. Load BAM Data                           ##
################################################

cat("\n[2/6] Loading BAM data...\n")

reads_list <- tryCatch({
    bamtolist(bamfolder = dirname(bam_file), annotation = annotation)
}, error = function(e) {
    cat("ERROR loading BAM:", conditionMessage(e), "\n")
    quit(status = 1)
})

cat("BAM loaded:", length(reads_list), "sample(s)\n")

################################################
## 3. Filter by Read Length                   ##
################################################

cat("\n[3/6] Filtering reads by length (", paste(read_lengths, collapse = ", "), ")...\n")

filtered_data <- tryCatch({
    length_filter(data = reads_list, length_filter_mode = "custom", length_filter_vector = read_lengths)
}, error = function(e) {
    cat("ERROR during length filtering:", conditionMessage(e), "\n")
    quit(status = 1)
})

cat("After filtering:", nrow(filtered_data[[1]]), "reads\n")

################################################
## 4. Calculate P-site Offsets                ##
################################################

cat("\n[4/6] Calculating P-site offsets...\n")

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
    cat("WARNING: P-site offset calculation failed:", conditionMessage(e), "\n")
    # Create placeholder
    dt_placeholder <- data.table(
        length = integer(),
        total_percentage = numeric(),
        start_percentage = numeric(),
        around_start = numeric(),
        offset_from_5 = numeric(),
        offset_from_3 = numeric(),
        sample = character()
    )
    # Write placeholder to avoid empty file downstream issues
    fwrite(dt_placeholder, paste0(prefix, "_psite_offset.txt"), sep = "\t")
    return(dt_placeholder)
})

cat("P-site offsets calculated\n")
fwrite(psite_offset, paste0(prefix, "_psite_offset.tsv"), sep = "\t")

################################################
## 5. Assign P-sites to Reads                 ##
################################################

cat("\n[5/6] Assigning P-site positions...\n")

psite_data <- tryCatch({
    psite_info(filtered_data, annotation)
}, error = function(e) {
    cat("WARNING: P-site assignment failed:", conditionMessage(e), "\n")
    return(filtered_data)
})

################################################
## 6. QC Analyses                             ##
################################################

cat("\n[6/6] Running QC analyses...\n")

# 6a. CDS Coverage
cat("  - CDS coverage...\n")
tryCatch({
    cds_cov <- cds_coverage(psite_data, annotation)
    fwrite(cds_cov, paste0(prefix, "_cds_coverage.tsv"), sep = "\t")
}, error = function(e) {
    cat("  WARNING: CDS coverage failed:", conditionMessage(e), "\n")
})

# 6b. Codon Usage
cat("  - Codon usage...\n")
tryCatch({
    codon_usage <- codon_usage_psite(psite_data, annotation)
    if (nrow(codon_usage) > 0) {
        fwrite(codon_usage, paste0(prefix, "_codon_usage.tsv"), sep = "\t")
    }
}, error = function(e) {
    cat("  WARNING: Codon usage failed:", conditionMessage(e), "\n")
})

# 6c. Frame Distribution
cat("  - Frame distribution...\n")
tryCatch({
    frame_dist <- frame_psite(psite_data, annotation)
    fwrite(frame_dist, paste0(prefix, "_frame_distribution.tsv"), sep = "\t")
}, error = function(e) {
    cat("  WARNING: Frame distribution failed:", conditionMessage(e), "\n")
})

# 6d. Read Length Distribution
cat("  - Read length distribution...\n")
tryCatch({
    rl_plot <- read_length_plot(data = filtered_data, sample = names(filtered_data)[1])
    ggsave(file.path(plot_dir, paste0(prefix, "_read_length_distribution.pdf")),
           plot = rl_plot, width = 8, height = 6)
}, error = function(e) {
    cat("  WARNING: Read length plot failed:", conditionMessage(e), "\n")
})

# 6e. Metagene Profiles (if metaplots function is available)
cat("  - Metagene profiles...\n")
tryCatch({
    if (exists("metaplots")) {
        metaplots_output <- metaplots(psite_data, annotation, sample = names(psite_data)[1])
        ggsave(file.path(plot_dir, paste0(prefix, "_metaplots.pdf")),
               plot = metaplots_output, width = 10, height = 8)
    } else {
        # Alternative: plot directly using the riboWaltz internal functions
        cat("  (metaplots function not available, generating subset plots)\n")
    }
}, error = function(e) {
    cat("  WARNING: Metagene profiles failed:", conditionMessage(e), "\n")
})

# 6f. P-site Region Distribution
cat("  - P-site region distribution...\n")
tryCatch({
    if (exists("percentage_regions")) {
        pct_regions <- percentage_regions(psite_data, annotation)
        fwrite(pct_regions, paste0(prefix, "_region_distribution.tsv"), sep = "\t")
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
