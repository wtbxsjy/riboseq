process RIBOSEQC_PREPAREANNOTATION {
    tag "${gtf.baseName}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
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
    #!/bin/bash

    cat <<EOF > script.R
    library(RiboseQC)

    # Prepare annotation files
    # forge_BSgenome=TRUE builds a temporary BSgenome package from FASTA
    # so ORFquant has access to genomic sequences for gene-level analysis.
    # Required for non-model organisms that lack pre-built BSgenome packages.
    prepare_annotation_files(
        annotation_directory = ".",
        genome_seq = "${fasta}",
        gtf_file = "${gtf}",
        scientific_name = "Genome.annotation",
        annotation_name = "custom",
        export_bed_tables_TxDb = FALSE,
        forge_BSgenome = TRUE,
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
    EOF

    # Use Rscript from the Conda environment if available
    if [[ -n "\$CONDA_PREFIX" ]]; then
        "\$CONDA_PREFIX/bin/Rscript" script.R
    else
        Rscript script.R
    fi
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
