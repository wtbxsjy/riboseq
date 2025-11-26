process RIBOSEQC_PREPAREANNOTATION {
    tag "${gtf.baseName}"
    label 'process_medium'

    conda "bioconda::riboseqc=1.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/riboseqc:1.1--r36_1' :
        'quay.io/biocontainers/riboseqc:1.1--r36_1' }"

    input:
    path gtf
    path fasta

    output:
    path "*_Rannot"    , emit: annotation
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = gtf.baseName
    """
    #!/usr/bin/env Rscript

    library(RiboseQC)

    # Prepare annotation files
    # Uses FASTA file directly (forge_BSgenome=FALSE to avoid package installation)
    prepare_annotation_files(
        annotation_directory = ".",
        genome_seq = "${fasta}",
        gtf_file = "${gtf}",
        scientific_name = "Genome.annotation",
        annotation_name = "custom",
        export_bed_tables_TxDb = FALSE,
        forge_BSgenome = FALSE,
        create_TxDb = TRUE
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
    def prefix = gtf.baseName
    """
    touch ${prefix}_Rannot

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        riboseqc: "1.1"
    END_VERSIONS
    """
}
