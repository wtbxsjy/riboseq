/*
 * Extract read length and P-site offset (cutoff) for ORFquant
 *
 * Priority: riboWaltz > RiboseQC > hardcoded defaults (28-32 -> 12)
 *
 * riboWaltz is preferred because its offsets are more consistent across
 * samples and cover a broader range of read lengths. RiboseQC serves as
 * fallback with its frame_preference-based quality metric. Hardcoded
 * defaults are the last resort.
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
    tuple val(meta), path(psites_calcs)                         // RiboseQC P_sites_calcs (fallback)
    tuple val(meta2), path(ribowaltz_psite), val(use_rw)        // riboWaltz psite_offset.tsv (primary)

    output:
    tuple val(meta), path("*_rl_cutoff.tsv"), emit: rl_cutoff
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def rw_file = (ribowaltz_psite && ribowaltz_psite.name && ribowaltz_psite.name != 'null' && use_rw)
        ? ribowaltz_psite.name
        : ''
    """
    #!/bin/bash
    set -euo pipefail

    cat <<'RSCRIPT' > extract_rl.R
    # Parse args: ribowaltz_file input_file output_file
    args <- commandArgs(trailingOnly = TRUE)
    rw_file    <- if (length(args) >= 3 && nchar(args[1]) > 0) args[1] else NA_character_
    input_file <- args[2]   # RiboseQC P_sites_calcs (fallback)
    output_file <- args[length(args)]

    # ---- helpers ----

    create_default_output <- function(output_file, reason) {
        cat("WARNING:", reason, "\n")
        cat("Using hardcoded defaults (28-32 nt -> cutoff 12)\n")
        default_data <- data.frame(
            read_length = c(28, 29, 30, 31, 32),
            cutoff = c(12, 12, 12, 12, 12),
            comp = "nucl",
            stringsAsFactors = FALSE
        )
        write.table(default_data, file = output_file, sep = "\t", row.names = FALSE, quote = FALSE)
        quit(status = 0)
    }

    # ---- PRIMARY: load riboWaltz offsets ----
    load_ribowaltz_offsets <- function(rw_path) {
        if (is.na(rw_path) || !file.exists(rw_path)) return(NULL)
        finfo <- file.info(rw_path)
        if (is.na(finfo[["size"]]) || finfo[["size"]] == 0) return(NULL)
        cat("Loading riboWaltz offsets from:", rw_path, "\n")
        rw <- tryCatch({
            read.table(rw_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
        }, error = function(e) {
            cat("  Failed:", conditionMessage(e), "\n")
            return(NULL)
        })
        if (is.null(rw) || nrow(rw) == 0) return(NULL)

        if (!"length" %in% colnames(rw)) {
            cat("  riboWaltz file missing 'length' column\n")
            return(NULL)
        }
        offset_col <- if ("corrected_offset_from_5" %in% colnames(rw)) {
            "corrected_offset_from_5"
        } else if ("offset_from_5" %in% colnames(rw)) {
            "offset_from_5"
        } else {
            cat("  riboWaltz file missing offset column\n")
            return(NULL)
        }

        result <- data.frame(
            read_length = rw[["length"]],
            cutoff = rw[[offset_col]],
            comp = "nucl",
            stringsAsFactors = FALSE
        )
        result <- result[!is.na(result[["cutoff"]]), ]
        result <- result[order(result[["read_length"]]), ]
        cat("  Converted", nrow(result), "riboWaltz offsets\n")
        return(result)
    }

    # ---- FALLBACK: load RiboseQC P_sites_calcs ----
    load_riboseqc_offsets <- function(pc_file) {
        if (!file.exists(pc_file)) {
            cat("RiboseQC P_sites_calcs not found:", pc_file, "\n")
            return(NULL)
        }
        finfo <- file.info(pc_file)
        if (is.na(finfo[["size"]]) || finfo[["size"]] < 10) {
            cat("RiboseQC P_sites_calcs too small / empty\n")
            return(NULL)
        }
        first_line <- tryCatch(readLines(pc_file, n = 1, warn = FALSE), error = function(e) "")
        if (length(first_line) == 0 || nchar(first_line) == 0 ||
            grepl("^# No P-site|^# Placeholder|placeholder|insufficient|failed",
                  first_line, ignore.case = TRUE)) {
            cat("RiboseQC P_sites_calcs indicates failure:", first_line, "\n")
            return(NULL)
        }

        cat("Loading RiboseQC P_sites_calcs:", pc_file, "\n")

        # Read as RDS or TSV
        psites <- tryCatch({
            readRDS(pc_file)
        }, error = function(e) {
            cat("  Not RDS, trying TSV...\n")
            tryCatch({
                read.table(pc_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
            }, error = function(e2) {
                cat("  Failed to parse RiboseQC file\n")
                return(NULL)
            })
        })
        if (is.null(psites)) return(NULL)

        # Extract data.frame from list if needed
        if (is.list(psites) && !is.data.frame(psites)) {
            if ("P_sites_calcs" %in% names(psites)) {
                psites <- psites[["P_sites_calcs"]]
            } else {
                for (nm in names(psites)) {
                    if (is.data.frame(psites[[nm]]) && "max_coverage" %in% colnames(psites[[nm]])) {
                        psites <- psites[[nm]]
                        break
                    }
                }
            }
        }
        if (!is.data.frame(psites)) {
            cat("  Could not extract P_sites_calcs data.frame\n")
            return(NULL)
        }

        # Normalize column names
        if ("rl" %in% colnames(psites) && !"read_length" %in% colnames(psites))
            colnames(psites)[colnames(psites) == "rl"] <- "read_length"
        if ("offset" %in% colnames(psites) && !"cutoff" %in% colnames(psites))
            colnames(psites)[colnames(psites) == "offset"] <- "cutoff"
        if ("compartment" %in% colnames(psites) && !"comp" %in% colnames(psites))
            colnames(psites)[colnames(psites) == "compartment"] <- "comp"

        required <- c("read_length", "cutoff", "max_coverage")
        if (!all(required %in% colnames(psites))) {
            cat("  Missing required columns:", paste(setdiff(required, colnames(psites)), collapse=", "), "\n")
            return(NULL)
        }

        # Filter max_coverage == TRUE rows
        if (is.logical(psites[["max_coverage"]])) {
            best <- psites[psites[["max_coverage"]] == TRUE, ]
        } else if (is.numeric(psites[["max_coverage"]])) {
            best <- psites[psites[["max_coverage"]] == 1, ]
        } else {
            best <- psites[as.logical(psites[["max_coverage"]]), ]
        }

        if (nrow(best) == 0) {
            cat("  No max_coverage rows, using all unique pairs\n")
            best <- psites[!duplicated(psites[["read_length"]]), ]
        }

        output_cols <- c("read_length", "cutoff")
        if ("comp" %in% colnames(best)) output_cols <- c(output_cols, "comp")
        result <- best[, output_cols, drop = FALSE]
        if (!"comp" %in% names(result)) result[["comp"]] <- "nucl"
        result <- result[order(result[["read_length"]]), ]
        cat("  Extracted", nrow(result), "RiboseQC offsets\n")
        return(result)
    }

    # ============================================================
    # Main: try riboWaltz first, then RiboseQC, then defaults
    # ============================================================

    cat("=== EXTRACT_RL_CUTOFF ===\n")
    cat("riboWaltz file:", if (is.na(rw_file)) "(none)" else rw_file, "\n")
    cat("RiboseQC file :", input_file, "\n")

    # Priority 1: riboWaltz
    offsets <- load_ribowaltz_offsets(rw_file)
    if (!is.null(offsets) && nrow(offsets) > 0) {
        cat(">>> Using riboWaltz offsets (primary)\n")
        write.table(offsets, file = output_file, sep = "\t", row.names = FALSE, quote = FALSE)
        cat("  Wrote", nrow(offsets), "read_length / cutoff pairs\n")
        quit(status = 0)
    }

    # Priority 2: RiboseQC
    offsets <- load_riboseqc_offsets(input_file)
    if (!is.null(offsets) && nrow(offsets) > 0) {
        cat(">>> Using RiboseQC offsets (fallback)\n")
        write.table(offsets, file = output_file, sep = "\t", row.names = FALSE, quote = FALSE)
        cat("  Wrote", nrow(offsets), "read_length / cutoff pairs\n")
        quit(status = 0)
    }

    # Priority 3: hardcoded defaults
    create_default_output(output_file, "Neither riboWaltz nor RiboseQC offsets available")
RSCRIPT

    # Run R script: args = rw_file, psites_calcs, output
    Rscript extract_rl.R "${rw_file}" "${psites_calcs}" "${prefix}_rl_cutoff.tsv"

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
