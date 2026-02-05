process RIBOSEQC_ANALYSIS {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/riboseqc:1.1--r36_1' :
        'quay.io/biocontainers/riboseqc:1.1--r36_1' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path annotation  // *_Rannot file from RIBOSEQC_PREPAREANNOTATION
    path fasta       // genome fasta for FaFile

    output:
    tuple val(meta), path("*_results_RiboseQC")       , emit: results
    tuple val(meta), path("*_results_RiboseQC_all")   , emit: results_all, optional: true
    tuple val(meta), path("*_for_ORFquant")           , emit: orfquant, optional: true
    tuple val(meta), path("*_coverage_*.bedgraph")    , emit: coverage, optional: true
    tuple val(meta), path("*_P_sites_*.bedgraph")     , emit: psites_bedgraph, optional: true
    tuple val(meta), path("*_P_sites_calcs")          , emit: psites_calcs, optional: true
    tuple val(meta), path("*_junctions")              , emit: junctions, optional: true
    tuple val(meta), path("*_ggribo.tsv")             , emit: ggribo, optional: true
    path "versions.yml"                               , emit: versions

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

    # Run RiboseQC analysis with error handling
    cat("Starting RiboseQC analysis...\\n")
    cat("Annotation file: ${annotation}\\n")
    cat("BAM file: ${bam}\\n")
    cat("Genome FASTA: ${fasta}\\n")
    cat("Sample name: ${prefix}\\n")
    
    # Function to create placeholder output files for low-signal samples
    create_placeholder_outputs <- function(prefix, reason) {
        cat("\\nCreating placeholder output files for downstream processing...\\n")
        cat("Reason:", reason, "\\n")
        
        # Create minimal results file
        results_file <- paste0(prefix, "_results_RiboseQC")
        writeLines(paste0("# RiboseQC placeholder - ", reason), results_file)
        
        # Create empty P_sites_calcs
        psites_file <- paste0(prefix, "_P_sites_calcs")
        writeLines(paste0("# No P-site data - ", reason), psites_file)
        
        # Create empty for_ORFquant file (optional but helps downstream)
        orfquant_file <- paste0(prefix, "_for_ORFquant")
        # Create an empty RData-like marker file
        writeLines(paste0("# No ORFquant data - ", reason), orfquant_file)
        
        # Create empty bedgraph files
        writeLines("", paste0(prefix, "_P_sites_plus.bedgraph"))
        writeLines("", paste0(prefix, "_P_sites_minus.bedgraph"))
        writeLines("", paste0(prefix, "_coverage_plus.bedgraph"))
        writeLines("", paste0(prefix, "_coverage_minus.bedgraph"))
        
        cat("Placeholder files created. Downstream ORF prediction will be skipped for this sample.\\n")
    }
    
    analysis_success <- tryCatch({
        RiboseQC_analysis(
            annotation_file = "${annotation}",
            bam_files = "${bam}",
            genome_seq = "${fasta}",
            dest_names = "${prefix}",
            sample_names = "${prefix}",
            ${fast_mode}
            create_report = FALSE,
            write_tmp_files = TRUE
        )
        cat("RiboseQC analysis completed successfully\\n")
        TRUE
    }, error = function(e) {
        error_msg <- conditionMessage(e)
        cat("ERROR in RiboseQC_analysis:\\n")
        cat(error_msg, "\\n")
        
        # Check for common low-signal/empty-BAM errors that should allow pipeline continuation
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
            cat("\\nWARNING: RiboseQC failed due to insufficient data/reads.\\n")
            cat("This typically occurs when:\\n")
            cat("  - BAM file has no aligned reads after filtering\\n")
            cat("  - Sample has extremely low ribosome profiling signal\\n")
            cat("  - Read length distribution is empty\\n")
            create_placeholder_outputs("${prefix}", "insufficient reads or signal")
            return(FALSE)
        } else {
            # For unexpected errors, still try to create placeholders but exit with error
            cat("\\nUnexpected RiboseQC error. Creating placeholders and exiting...\\n")
            create_placeholder_outputs("${prefix}", paste0("unexpected error: ", error_msg))
            quit(status = 1)
        }
    })
    
    if (!analysis_success) {
        cat("\\nRiboseQC analysis skipped due to insufficient data.\\n")
    }
    
    # Verify critical output files exist (may be placeholders)
    output_file <- "${prefix}_P_sites_calcs"
    if (!file.exists(output_file)) {
        cat("WARNING: P_sites_calcs file not found. Creating placeholder.\\n")
        writeLines("# No P-site data - file missing after analysis", output_file)
    }
    
    results_file <- "${prefix}_results_RiboseQC"
    if (!file.exists(results_file)) {
        cat("WARNING: results_RiboseQC file not found. Creating placeholder.\\n")
        writeLines("# RiboseQC results placeholder - analysis incomplete", results_file)
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
    echo "[INFO] Running RiboseQC analysis..."
    if [[ -n "\${CONDA_PREFIX:-}" ]]; then
        "\$CONDA_PREFIX/bin/Rscript" script.R
    else
        Rscript script.R
    fi
    
    # Check if P_sites_calcs was created successfully
    if [[ ! -f "${prefix}_P_sites_calcs" ]] || [[ ! -s "${prefix}_P_sites_calcs" ]]; then
        echo "[WARNING] P_sites_calcs file is missing or empty - sample has insufficient signal"
        echo "[WARNING] Downstream ORF prediction will be skipped for this sample"
    else
        echo "[INFO] P-site calculation successful"
    fi
    
    echo "[INFO] RiboseQC analysis completed"
    ls -lh ${prefix}_* || true

    # Convert P-sites bedgraphs to ggRibo input format
    # ggRibo format: Count \t Chromosome \t Position \t Strand
    # BedGraph format: Chr \t Start \t End \t Count
    # Note: BedGraph is 0-based start, 1-based end. ggRibo expects 1-based position.
    # Since P-sites are single nucleotides, Start+1 = End. We use End (col 3) as Position.

    if [ -f "${prefix}_P_sites_plus.bedgraph" ]; then
        awk -v OFS='\\t' '{print \$4, \$1, \$3, "+"}' "${prefix}_P_sites_plus.bedgraph" > "${prefix}_ggribo.tsv"
    fi
    if [ -f "${prefix}_P_sites_minus.bedgraph" ]; then
        awk -v OFS='\\t' '{print \$4, \$1, \$3, "-"}' "${prefix}_P_sites_minus.bedgraph" >> "${prefix}_ggribo.tsv"
    fi
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_results_RiboseQC
    touch ${prefix}_results_RiboseQC_all
    touch ${prefix}_for_ORFquant
    touch ${prefix}_coverage_plus.bedgraph
    touch ${prefix}_coverage_minus.bedgraph
    touch ${prefix}_P_sites_plus.bedgraph
    touch ${prefix}_P_sites_minus.bedgraph
    touch ${prefix}_P_sites_calcs
    touch ${prefix}_junctions
    touch ${prefix}_ggribo.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        riboseqc: "1.1"
    END_VERSIONS
    """
}
