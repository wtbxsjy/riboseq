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
    
    tryCatch({
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
    }, error = function(e) {
        cat("ERROR in RiboseQC_analysis:\\n")
        cat(conditionMessage(e), "\\n")
        quit(status = 1)
    })
    
    # Verify critical output files exist and have content
    output_file <- "${prefix}_P_sites_calcs"
    if (!file.exists(output_file)) {
        cat("ERROR: Required output file not found:", output_file, "\\n")
        quit(status = 1)
    }
    finfo <- file.info(output_file)
    if (is.na(finfo[["size"]]) || finfo[["size"]] < 10) {
        cat("ERROR: Output file is empty or too small:", output_file, "(", finfo[["size"]], "bytes)\\n")
        quit(status = 1)
    }
    cat("Verified output:", output_file, "(", finfo[["size"]], "bytes)\\n")

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
    if [[ -n "\$CONDA_PREFIX" ]]; then
        "\$CONDA_PREFIX/bin/Rscript" script.R
    else
        Rscript script.R
    fi
    
    # Check if critical files were created
    if [[ ! -f "${prefix}_P_sites_calcs" ]] || [[ ! -s "${prefix}_P_sites_calcs" ]]; then
        echo "[ERROR] P_sites_calcs file is missing or empty"
        exit 1
    fi
    
    echo "[INFO] RiboseQC analysis completed"
    ls -lh ${prefix}_*

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
