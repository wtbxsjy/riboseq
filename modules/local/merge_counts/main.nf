process MERGE_COUNTS {
    tag "merge_counts"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/r-data.table:latest' :
        'community.wave.seqera.io/library/r-data.table:latest' }"

    input:
    path count_files          // all *counts.tsv files from featureCounts
    path sample_sheet_csv     // sample sheet: sample,type,treatment,...

    output:
    path "merged_counts.tsv"  , emit: counts
    path "sample_sheet.csv"   , emit: sample_sheet
    path "versions.yml"       , emit: versions

    script:
    """
    #!/bin/bash
    # Merge featureCounts output files into a single count matrix (genes x samples)
    # featureCounts output: Geneid \t Chr \t Start \t End \t Strand \t Length \t count
    # We extract Geneid + count column (last column) from each file and merge on Geneid

    # Build R script for merging
    cat <<'RSCRIPT' > merge_counts.R
    # Get input files from glob pattern or space-separated path list
    input_str <- "${count_files}"
    if (grepl(" ", input_str)) {
        count_files <- strsplit(input_str, " ")[[1]]
        count_files <- count_files[file.exists(count_files)]
    } else {
        count_files <- Sys.glob(input_str)
    }
    if (length(count_files) == 0) {
        count_files <- list.files(pattern = "_counts\\\\.tsv\$", recursive = TRUE, full.names = TRUE)
    }

    cat("Found", length(count_files), "count files\\n")

    # Read and merge
    mat_list <- lapply(count_files, function(f) {
        # Extract sample name from filename
        sample_name <- gsub("_counts\\\\.tsv\$", "", basename(f))
        cat("  Reading:", sample_name, "\\n")

        x <- read.delim(f, header = TRUE, comment.char = "#", stringsAsFactors = FALSE)
        # Keep only Geneid and count (last column)
        count_col <- ncol(x)
        out <- x[, c("Geneid", colnames(x)[count_col]), drop = FALSE]
        colnames(out)[2] <- sample_name
        out
    })

    # Merge all on Geneid
    merged <- mat_list[[1]]
    if (length(mat_list) > 1) {
        for (i in 2:length(mat_list)) {
            merged <- merge(merged, mat_list[[i]], by = "Geneid", all = TRUE)
        }
    }
    merged[is.na(merged)] <- 0
    colnames(merged)[1] <- "gene_id"

    write.table(merged, "merged_counts.tsv", sep = "\\t", row.names = FALSE, quote = FALSE)
    cat("Merged matrix:", nrow(merged), "genes x", ncol(merged) - 1, "samples\\n")

    # Copy sample sheet
    file.copy("${sample_sheet_csv}", "sample_sheet.csv", overwrite = TRUE)

    # Write versions
    writeLines(c(
        '"${task.process}":',
        '    r-base: "4.0"'
    ), "versions.yml")
    quit(save = "no")
RSCRIPT

    Rscript merge_counts.R
    """

    stub:
    """
    touch merged_counts.tsv
    touch sample_sheet.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: "4.0"
    END_VERSIONS
    """
}
