/*
 * Extract read length and P-site offset (cutoff) from RiboseQC P_sites_calcs output
 * 
 * This module extracts rows where max_coverage is TRUE from the RiboseQC
 * P_sites_calcs file, keeping read_length, cutoff, and comp columns.
 * The resulting file is used as input for prepare_for_ORFquant to improve
 * P-site offset accuracy.
 * 
 * Output format (tab-separated):
 *   read_length  cutoff  comp
 *   28           12      nucl
 *   27           12      nucl
 *   ...
 */

process EXTRACT_RL_CUTOFF {
    tag "$meta.id"
    label 'process_low'

    // Use same container as RiboseQC (contains R environment)
    conda "../../riboseqc/analysis/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/riboseqc:1.1--r36_1' :
        'quay.io/biocontainers/riboseqc:1.1--r36_1' }"

    input:
    tuple val(meta), path(psites_calcs)  // *_P_sites_calcs file from RiboseQC

    output:
    tuple val(meta), path("*_rl_cutoff.tsv"), emit: rl_cutoff
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/usr/bin/env Rscript

    # Read P_sites_calcs file
    # The file is typically an RDS file from RiboseQC containing a data frame
    psites_data <- tryCatch({
        readRDS("${psites_calcs}")
    }, error = function(e) {
        # If not RDS, try reading as tab-delimited text
        read.table("${psites_calcs}", header = TRUE, sep = "\\t", stringsAsFactors = FALSE)
    })

    # If the data is a list (from RDS), extract the P_sites_calcs data frame
    if (is.list(psites_data) && !is.data.frame(psites_data)) {
        # RiboseQC stores results in a list structure
        # Look for P_sites_calcs or similar component
        if ("P_sites_calcs" %in% names(psites_data)) {
            psites_data <- psites_data[["P_sites_calcs"]]
        } else if ("read_stats" %in% names(psites_data)) {
            # Alternative structure in some RiboseQC versions
            psites_data <- psites_data[["read_stats"]]
        } else {
            # Try to find the correct component
            for (name in names(psites_data)) {
                if (is.data.frame(psites_data[[name]]) && 
                    "max_coverage" %in% colnames(psites_data[[name]])) {
                    psites_data <- psites_data[[name]]
                    break
                }
            }
        }
    }

    # Check and normalize column names
    # The file may have different column name variations
    
    # Handle read_length variations
    if ("rl" %in% colnames(psites_data) && !"read_length" %in% colnames(psites_data)) {
        colnames(psites_data)[colnames(psites_data) == "rl"] <- "read_length"
    }
    
    # Handle cutoff variations (might be called 'offset' in some versions)
    if ("offset" %in% colnames(psites_data) && !"cutoff" %in% colnames(psites_data)) {
        colnames(psites_data)[colnames(psites_data) == "offset"] <- "cutoff"
    }
    
    # Handle compartment variations - ORFquant expects 'comp' not 'compartment'
    if ("compartment" %in% colnames(psites_data) && !"comp" %in% colnames(psites_data)) {
        colnames(psites_data)[colnames(psites_data) == "compartment"] <- "comp"
    }

    # Check required columns exist
    required_cols <- c("read_length", "cutoff", "max_coverage")
    missing_cols <- setdiff(required_cols, colnames(psites_data))
    if (length(missing_cols) > 0) {
        stop(paste("Missing required columns:", paste(missing_cols, collapse = ", "), 
                   "\\nAvailable columns:", paste(colnames(psites_data), collapse = ", ")))
    }

    # Debug: Print data structure
    cat("Data structure:\\n")
    cat("Rows:", nrow(psites_data), "\\n")
    cat("Columns:", paste(colnames(psites_data), collapse = ", "), "\\n")
    cat("max_coverage class:", class(psites_data[["max_coverage"]]), "\\n")
    cat("max_coverage unique values:", paste(unique(psites_data[["max_coverage"]]), collapse = ", "), "\\n")
    
    # Filter rows where max_coverage is TRUE
    # Handle both logical (TRUE/FALSE) and numeric (1/0) representations
    if (is.logical(psites_data[["max_coverage"]])) {
        filtered_data <- psites_data[psites_data[["max_coverage"]] == TRUE, ]
    } else if (is.numeric(psites_data[["max_coverage"]])) {
        # Treat 1 as TRUE, others as FALSE
        filtered_data <- psites_data[psites_data[["max_coverage"]] == 1, ]
    } else {
        # Try to coerce to logical
        filtered_data <- psites_data[as.logical(psites_data[["max_coverage"]]), ]
    }

    if (nrow(filtered_data) == 0) {
        # Fallback: If no max_coverage rows, use all unique read_length/cutoff pairs
        cat("WARNING: No rows with max_coverage == TRUE found. Using all unique read_length/cutoff pairs.\\n")
        
        # Remove duplicates based on read_length
        filtered_data <- psites_data[!duplicated(psites_data[["read_length"]]), ]
        
        if (nrow(filtered_data) == 0) {
            stop("No valid data found in P_sites_calcs file")
        }
    } else {
        cat("Found", nrow(filtered_data), "rows with max_coverage indicator\\n")
    }

    # Select columns: read_length, cutoff, and comp (if exists)
    # ORFquant expects: read_length, cutoff, comp
    output_cols <- c("read_length", "cutoff")
    if ("comp" %in% colnames(filtered_data)) {
        output_cols <- c(output_cols, "comp")
    } else {
        # Add default compartment 'nucl' if not present
        filtered_data[["comp"]] <- "nucl"
        output_cols <- c(output_cols, "comp")
    }

    output_data <- filtered_data[, output_cols, drop = FALSE]

    # Sort by read_length for consistency
    output_data <- output_data[order(output_data[["read_length"]]), ]

    # Write output in format expected by ORFquant prepare_for_ORFquant()
    write.table(
        output_data,
        file = "${prefix}_rl_cutoff.tsv",
        sep = "\\t",
        row.names = FALSE,
        quote = FALSE
    )

    cat("Extracted", nrow(output_data), "read length / cutoff pairs for ORFquant\\n")

    # Write versions
    writeLines(
        c(
            '"${task.process}":',
            paste0('    r-base: "', R.version[["major"]], ".", R.version[["minor"]], '"')
        ),
        "versions.yml"
    )
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo -e "read_length\\tcutoff\\tcomp" > ${prefix}_rl_cutoff.tsv
    echo -e "28\\t12\\tnucl" >> ${prefix}_rl_cutoff.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: "4.1.0"
    END_VERSIONS
    """
}
