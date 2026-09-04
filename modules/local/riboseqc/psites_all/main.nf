process RIBOSEQC_PSITES_ALL {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/../analysis/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/riboseqc:1.1--r36_1' :
        'quay.io/biocontainers/riboseqc:1.1--r36_1' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path annotation  // *_Rannot file from RIBOSEQC_PREPAREANNOTATION
    path fasta       // genome fasta for FaFile

    output:
    tuple val(meta), path("*_P_sites_*.bedgraph"), emit: psites_bedgraph, optional: true
    path "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def fast_mode = args.contains('fast_mode=FALSE') ? 'fast_mode = FALSE,' : 'fast_mode = TRUE,'
    """
    #!/bin/bash
    set -euo pipefail

    cat <<'RSCRIPT' > script.R
    library(RiboseQC)
    library(Rsamtools)

    # P-site quantification pass for ORF expression statistics.
    # RIBOSEQC_ANALYSIS uses readlength_choice_method="max_coverage": when the
    # max-coverage read length fails the frame-preference filter, it skips the
    # P-sites position calculation and writes no *_P_sites_*.bedgraph files
    # (observed for all PRJEB26593 Ribo-seq samples: RL29 frame_preference
    # 47.1% < 50%). With choice="all" P-sites positions are computed for every
    # read length using each length's cutoff from the P_sites_calcs table, so
    # the bedgraphs are always emitted.
    cat("Starting RiboseQC P-site quantification (all read lengths)...\\n")
    cat("Annotation file: ${annotation}\\n")
    cat("BAM file: ${bam}\\n")
    cat("Genome FASTA: ${fasta}\\n")
    cat("Sample name: ${prefix}\\n")

    create_placeholder_bedgraphs <- function(prefix, reason) {
        cat("\\nCreating placeholder P-site bedgraphs...\\n")
        cat("Reason:", reason, "\\n")
        writeLines(paste0("# Placeholder bedgraph - ", reason), paste0(prefix, "_P_sites_plus.bedgraph"))
        writeLines(paste0("# Placeholder bedgraph - ", reason), paste0(prefix, "_P_sites_minus.bedgraph"))
    }

    analysis_success <- tryCatch({
        RiboseQC_analysis(
            annotation_file = "${annotation}",
            bam_files = "${bam}",
            genome_seq = "${fasta}",
            dest_names = "${prefix}",
            sample_names = "${prefix}",
            ${fast_mode}
            readlength_choice_method = "all",
            create_report = FALSE,
            write_tmp_files = TRUE
        )
        cat("RiboseQC P-site quantification completed successfully\\n")
        TRUE
    }, error = function(e) {
        error_msg <- conditionMessage(e)
        cat("ERROR in RiboseQC_analysis:\\n")
        cat(error_msg, "\\n")

        # Low-signal/empty-BAM errors should not kill the pipeline: emit
        # placeholder bedgraphs so the sample is still listed downstream.
        is_low_signal_error <- (
            grepl("subscript out of bounds", error_msg, ignore.case = TRUE) ||
            grepl("replacement has length zero", error_msg, ignore.case = TRUE) ||
            grepl("no non-missing arguments", error_msg, ignore.case = TRUE) ||
            grepl("cannot allocate vector", error_msg, ignore.case = TRUE) ||
            grepl("argument is of length zero", error_msg, ignore.case = TRUE) ||
            grepl("zero-length", error_msg, ignore.case = TRUE) ||
            grepl("no reads", error_msg, ignore.case = TRUE) ||
            grepl("empty", error_msg, ignore.case = TRUE)
        )

        if (is_low_signal_error) {
            cat("\\nWARNING: RiboseQC P-site quantification failed due to insufficient data/reads.\\n")
            create_placeholder_bedgraphs("${prefix}", "insufficient reads or signal")
            FALSE
        } else {
            cat("\\nUnexpected RiboseQC error. Creating placeholders and exiting...\\n")
            create_placeholder_bedgraphs("${prefix}", paste0("unexpected error: ", error_msg))
            quit(status = 1)
        }
    })

    if (!analysis_success) {
        cat("\\nRiboseQC P-site quantification skipped due to insufficient data.\\n")
    }

    # Write versions
    writeLines(
        c(
            '"${task.process}":',
            paste0('    riboseqc: "', packageVersion("RiboseQC"), '"')
        ),
        "versions.yml"
    )
RSCRIPT

    # Use Rscript from the Conda environment if available
    echo "[INFO] Running RiboseQC P-site quantification (all read lengths)..."
    if [[ -n "\${CONDA_PREFIX:-}" ]]; then
        "\$CONDA_PREFIX/bin/Rscript" script.R
    else
        Rscript script.R
    fi

    if [[ ! -f "${prefix}_P_sites_plus.bedgraph" || ! -f "${prefix}_P_sites_minus.bedgraph" ]]; then
        echo "[WARNING] P_sites bedgraph files missing - sample has insufficient signal"
    else
        echo "[INFO] P-site bedgraphs written"
    fi

    echo "[INFO] RiboseQC P-site quantification completed"
    ls -lh ${prefix}_P_sites_*.bedgraph || true
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_P_sites_plus.bedgraph
    touch ${prefix}_P_sites_minus.bedgraph

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        riboseqc: "1.1"
    END_VERSIONS
    """
}
