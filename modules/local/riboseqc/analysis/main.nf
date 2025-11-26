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
    tuple val(meta), path("*_results_RiboseQC")    , emit: results
    tuple val(meta), path("*_RiboseQC_report.html"), emit: html, optional: true
    tuple val(meta), path("*_for_ORFquant")        , emit: orfquant, optional: true
    path "versions.yml"                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def create_report = args.contains('create_report=FALSE') ? '' : 'create_report = TRUE,'
    def fast_mode = args.contains('fast_mode=FALSE') ? 'fast_mode = FALSE,' : 'fast_mode = TRUE,'
    """
    #!/usr/bin/env Rscript

    library(RiboseQC)
    library(Rsamtools)

    # Run RiboseQC analysis
    # Note: genome_seq should be a file path string, not an FaFile object
    # The function will internally create FaFile and FaFile_Circ objects
    RiboseQC_analysis(
        annotation_file = "${annotation}",
        bam_files = "${bam}",
        genome_seq = "${fasta}",
        dest_names = "${prefix}",
        sample_names = "${prefix}",
        report_file = "${prefix}_RiboseQC_report.html",
        ${fast_mode}
        ${create_report}
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
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_results_RiboseQC
    touch ${prefix}_RiboseQC_report.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        riboseqc: "1.1"
    END_VERSIONS
    """
}
