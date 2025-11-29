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
    path "versions.yml"                               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def fast_mode = args.contains('fast_mode=FALSE') ? 'fast_mode = FALSE,' : 'fast_mode = TRUE,'
    """
    #!/bin/bash

    cat <<EOF > script.R
    library(RiboseQC)
    library(Rsamtools)

    # Run RiboseQC analysis
    # Note: genome_seq should be a file path string, not an FaFile object
    # The function will internally create FaFile and FaFile_Circ objects
    # HTML report generation is disabled due to compatibility issues with
    # rmarkdown in containerized environments (RiboseQC 1.1)
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

    # Write versions
    writeLines(
        c(
            '"${task.process}":',
            paste0('    riboseqc: "', packageVersion("RiboseQC"), '"')
        ),
        "versions.yml"
    )
    EOF

    # Use Rscript from the Conda environment if available
    if [[ -n "\$CONDA_PREFIX" ]]; then
        "\$CONDA_PREFIX/bin/Rscript" script.R
    else
        Rscript script.R
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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        riboseqc: "1.1"
    END_VERSIONS
    """
}
