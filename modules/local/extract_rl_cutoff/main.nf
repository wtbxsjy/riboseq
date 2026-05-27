/*
 * Extract read length and P-site offset (cutoff) from RiboseQC P_sites_calcs output
 *
 * This module extracts rows where max_coverage is TRUE from the RiboseQC
 * P_sites_calcs file, keeping read_length, cutoff, and comp columns.
 * The resulting file is used as input for prepare_for_ORFquant to improve
 * P-site offset accuracy.
 *
 * Fallback: When RiboseQC P_sites_calcs is empty/invalid, riboWaltz-derived
 * offsets (converted from psite_offset.tsv) are used if available.
 * Hardcoded defaults (28-32 → 12) are the last resort.
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

    conda "../../riboseqc/analysis/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/riboseqc:1.1--r36_1' :
        'quay.io/biocontainers/riboseqc:1.1--r36_1' }"

    input:
    tuple val(meta), path(psites_calcs)                         // RiboseQC P_sites_calcs (primary)
    tuple val(meta2), path(ribowaltz_psite), val(use_rw)        // riboWaltz psite_offset.tsv (optional fallback)
                                                                 // use_rw: boolean flag from params

    output:
    tuple val(meta), path("*_rl_cutoff.tsv"), emit: rl_cutoff
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // Determine if riboWaltz fallback file exists
    def rw_file = (ribowaltz_psite && ribowaltz_psite.name && ribowaltz_psite.name != 'null' && use_rw)
        ? ribowaltz_psite.name
        : ''
    """
    #!/bin/bash
    set -euo pipefail

    cat <<'RSCRIPT' > extract_rl.R
    # Parse args: input_file [ribowaltz_file] output_file
    args <- commandArgs(trailingOnly = TRUE)
    input_file <- args[1]
    rw_file <- if (length(args) >= 3 && nchar(args[2]) > 0) args[2] else NA_character_
    output_file <- args[length(args)]

    create_default_output <- function(output_file, reason) {
        cat("WARNING:", reason, "\\n")
        cat("Creating default rl_cutoff with common ribosome footprint lengths\\n")
        default_data <- data.frame(
            read_length = c(28, 29, 30, 31, 32),
            cutoff = c(12, 12, 12, 12, 12),
            comp = c("nucl", "nucl", "nucl", "nucl", "nucl"),
            stringsAsFactors = FALSE
        )
        write.table(default_data, file = output_file, sep = "\\t", row.names = FALSE, quote = FALSE)
        cat("Created default rl_cutoff with", nrow(default_data), "read lengths\\n")
        quit(status = 0)
    }

    # Try to load riboWaltz offsets as fallback
    load_ribowaltz_offsets <- function(rw_path) {
        if (is.na(rw_path) || !file.exists(rw_path)) return(NULL)
        finfo <- file.info(rw_path)
        if (is.na(finfo[["size"]]) || finfo[["size"]] == 0) return(NULL)
        cat("Loading riboWaltz offsets from:", rw_path, "\\n")
        rw <- tryCatch({
            read.table(rw_path, header = TRUE, sep = "\\t", stringsAsFactors = FALSE)
        }, error = function(e) {
            cat("Failed to read riboWaltz file:", conditionMessage(e), "\\n")
            return(NULL)
        })
        if (is.null(rw) || nrow(rw) == 0) return(NULL)

        # riboWaltz columns: length, corrected_offset_from_5, corrected_offset_from_3, sample
        # Convert to rl_cutoff format: read_length, cutoff, comp
        if (!"length" %in% colnames(rw)) {
            cat("riboWaltz file missing 'length' column\\n")
            return(NULL)
        }
        # Use corrected_offset_from_5 if available, otherwise offset_from_5
        offset_col <- if ("corrected_offset_from_5" %in% colnames(rw)) {
            "corrected_offset_from_5"
        } else if ("offset_from_5" %in% colnames(rw)) {
            "offset_from_5"
        } else {
            cat("riboWaltz file missing offset column\\n")
            return(NULL)
        }

        result <- data.frame(
            read_length = rw[["length"]],
            cutoff = rw[[offset_col]],
            comp = "nucl",
            stringsAsFactors = FALSE
        )
        # Remove rows with NA cutoff
        result <- result[!is.na(result[["cutoff"]]), ]
        result <- result[order(result[["read_length"]]), ]
        cat("Converted", nrow(result), "riboWaltz offsets to rl_cutoff format\\n")
        return(result)
    }

    # Check if RiboseQC input is valid
    input_valid <- TRUE
    if (!file.exists(input_file)) {
        cat("Input file not found:", input_file, "\\n")
        input_valid <- FALSE
    } else {
        finfo <- file.info(input_file)
        if (is.na(finfo[["size"]]) || finfo[["size"]] == 0) {
            cat("Input file is empty:", input_file, "\\n")
            input_valid <- FALSE
        } else if (finfo[["size"]] < 10) {
            cat("Input file too small (", finfo[["size"]], " bytes), likely placeholder/corrupted:", input_file, "\\n")
            input_valid <- FALSE
        }
    }

    # Check for placeholder content in first line
    if (input_valid) {
        first_line <- tryCatch({
            readLines(input_file, n = 1, warn = FALSE)
        }, error = function(e) { "" })
        if (length(first_line) == 0 || nchar(first_line) == 0 ||
            grepl("^# No P-site data|^# Placeholder|placeholder|insufficient|failed",
                  first_line, ignore.case = TRUE)) {
            cat("Input file indicates P-site calculation failed:", first_line, "\\n")
            input_valid <- FALSE
        }
    }

    # If RiboseQC data is invalid, try riboWaltz fallback
    if (!input_valid) {
        cat("RiboseQC P_sites_calcs is invalid, trying riboWaltz fallback...\\n")
        rw_data <- load_ribowaltz_offsets(rw_file)
        if (!is.null(rw_data) && nrow(rw_data) > 0) {
            write.table(rw_data, file = output_file, sep = "\\t", row.names = FALSE, quote = FALSE)
            cat("Used riboWaltz-derived offsets (", nrow(rw_data), "read lengths)\\n")
            quit(status = 0)
        }
        create_default_output(output_file, "No valid RiboseQC P_sites_calcs and no riboWaltz fallback available")
    }

    cat("Reading P_sites_calcs file:", input_file, "\\n")
    cat("File size:", finfo[["size"]], "bytes\\n")

    # Read P_sites_calcs file (RDS format from RiboseQC)
    psites_data <- tryCatch({
        data <- readRDS(input_file)
        cat("Successfully read RDS file\\n")
        data
    }, error = function(e) {
        cat("Failed to read as RDS:", conditionMessage(e), "\\n")
        cat("Attempting to read as tab-delimited text...\\n")
        tryCatch({
            read.table(input_file, header = TRUE, sep = "\\t", stringsAsFactors = FALSE)
        }, error = function(e2) {
            cat("Failed to read file as RDS or TSV.\\n")
            # Try riboWaltz fallback
            rw_data <- load_ribowaltz_offsets(rw_file)
            if (!is.null(rw_data) && nrow(rw_data) > 0) {
                write.table(rw_data, file = output_file, sep = "\\t", row.names = FALSE, quote = FALSE)
                cat("Used riboWaltz-derived offsets after RDS/TSV read failure\\n")
                quit(status = 0)
            }
            create_default_output(output_file, "Unable to parse RiboseQC file")
        })
    })

    cat("Data object class:", class(psites_data), "\\n")
    if (is.list(psites_data) && !is.data.frame(psites_data)) {
        cat("List names:", paste(names(psites_data), collapse = ", "), "\\n")
    }

    # If the data is a list, extract the P_sites_calcs data frame
    if (is.list(psites_data) && !is.data.frame(psites_data)) {
        cat("Data is a list, attempting to extract P_sites_calcs data frame...\\n")
        if ("P_sites_calcs" %in% names(psites_data)) {
            psites_data <- psites_data[["P_sites_calcs"]]
            cat("Extracted P_sites_calcs component\\n")
        } else if ("read_stats" %in% names(psites_data)) {
            psites_data <- psites_data[["read_stats"]]
            cat("Extracted read_stats component\\n")
        } else {
            cat("Searching for data frame with max_coverage column...\\n")
            found <- FALSE
            for (name in names(psites_data)) {
                obj <- psites_data[[name]]
                if (is.data.frame(obj) && "max_coverage" %in% colnames(obj)) {
                    psites_data <- obj
                    cat("Found suitable data frame in component:", name, "\\n")
                    found <- TRUE
                    break
                }
            }
            if (!found) {
                cat("ERROR: Could not find P_sites_calcs data frame\\n")
                # Try riboWaltz fallback
                rw_data <- load_ribowaltz_offsets(rw_file)
                if (!is.null(rw_data) && nrow(rw_data) > 0) {
                    write.table(rw_data, file = output_file, sep = "\\t", row.names = FALSE, quote = FALSE)
                    cat("Used riboWaltz-derived offsets after component search failure\\n")
                    quit(status = 0)
                }
                create_default_output(output_file, "P_sites_calcs data frame not found in RDS file")
            }
        }
    }

    if (!is.data.frame(psites_data)) {
        rw_data <- load_ribowaltz_offsets(rw_file)
        if (!is.null(rw_data) && nrow(rw_data) > 0) {
            write.table(rw_data, file = output_file, sep = "\\t", row.names = FALSE, quote = FALSE)
            cat("Used riboWaltz-derived offsets (P_sites_calcs not a data frame)\\n")
            quit(status = 0)
        }
        create_default_output(output_file, paste("Expected data frame, but got:", class(psites_data)))
    }

    cat("Final data frame:", nrow(psites_data), "rows x", ncol(psites_data), "columns\\n")
    cat("Column names:", paste(colnames(psites_data), collapse = ", "), "\\n")

    # Normalize column names
    if ("rl" %in% colnames(psites_data) && !"read_length" %in% colnames(psites_data)) {
        colnames(psites_data)[colnames(psites_data) == "rl"] <- "read_length"
    }
    if ("offset" %in% colnames(psites_data) && !"cutoff" %in% colnames(psites_data)) {
        colnames(psites_data)[colnames(psites_data) == "offset"] <- "cutoff"
    }
    if ("compartment" %in% colnames(psites_data) && !"comp" %in% colnames(psites_data)) {
        colnames(psites_data)[colnames(psites_data) == "compartment"] <- "comp"
    }

    # Check required columns
    required_cols <- c("read_length", "cutoff", "max_coverage")
    missing_cols <- setdiff(required_cols, colnames(psites_data))
    if (length(missing_cols) > 0) {
        rw_data <- load_ribowaltz_offsets(rw_file)
        if (!is.null(rw_data) && nrow(rw_data) > 0) {
            write.table(rw_data, file = output_file, sep = "\\t", row.names = FALSE, quote = FALSE)
            cat("Used riboWaltz-derived offsets (missing columns in P_sites_calcs)\\n")
            quit(status = 0)
        }
        create_default_output(output_file, paste("Missing required columns:", paste(missing_cols, collapse = ", "),
                   "\\nAvailable columns:", paste(colnames(psites_data), collapse = ", ")))
    }

    cat("max_coverage class:", class(psites_data[["max_coverage"]]), "\\n")
    cat("max_coverage unique values:", paste(unique(psites_data[["max_coverage"]]), collapse = ", "), "\\n")

    # Filter rows where max_coverage is TRUE
    if (is.logical(psites_data[["max_coverage"]])) {
        filtered_data <- psites_data[psites_data[["max_coverage"]] == TRUE, ]
    } else if (is.numeric(psites_data[["max_coverage"]])) {
        filtered_data <- psites_data[psites_data[["max_coverage"]] == 1, ]
    } else {
        filtered_data <- psites_data[as.logical(psites_data[["max_coverage"]]), ]
    }

    if (nrow(filtered_data) == 0) {
        cat("WARNING: No rows with max_coverage indicator. Using all unique read_length/cutoff pairs.\\n")
        filtered_data <- psites_data[!duplicated(psites_data[["read_length"]]), ]
        if (nrow(filtered_data) == 0) {
            rw_data <- load_ribowaltz_offsets(rw_file)
            if (!is.null(rw_data) && nrow(rw_data) > 0) {
                write.table(rw_data, file = output_file, sep = "\\t", row.names = FALSE, quote = FALSE)
                cat("Used riboWaltz-derived offsets (no valid data in P_sites_calcs)\\n")
                quit(status = 0)
            }
            create_default_output(output_file, "No valid data found in P_sites_calcs file")
        }
    } else {
        cat("Found", nrow(filtered_data), "rows with max_coverage indicator\\n")
    }

    # Select output columns
    output_cols <- c("read_length", "cutoff")
    if ("comp" %in% colnames(filtered_data)) {
        output_cols <- c(output_cols, "comp")
    } else {
        filtered_data[["comp"]] <- "nucl"
        output_cols <- c(output_cols, "comp")
    }

    output_data <- filtered_data[, output_cols, drop = FALSE]
    output_data <- output_data[order(output_data[["read_length"]]), ]

    write.table(output_data, file = output_file, sep = "\\t", row.names = FALSE, quote = FALSE)
    cat("Extracted", nrow(output_data), "read length / cutoff pairs for ORFquant\\n")
RSCRIPT

    # Run R script with arguments
    Rscript extract_rl.R "${psites_calcs}" "${rw_file}" "${prefix}_rl_cutoff.tsv"

    # Write versions
    R_VERSION=\$(Rscript -e 'cat(paste(R.version[["major"]], R.version[["minor"]], sep="."))')
    cat > versions.yml <<EOF
"${task.process}":
    r-base: "\$R_VERSION"
EOF
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
