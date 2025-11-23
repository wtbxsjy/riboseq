process BOWTIE_BUILD {
    tag "$fasta"
    label 'process_high'

    conda "bioconda::bowtie=1.3.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bowtie:1.3.1--py38h828cd81_9' :
        'quay.io/biocontainers/bowtie:1.3.1--py38h828cd81_9' }"

    input:
    tuple val(meta), path(fasta)

    output:
    path "bowtie", emit: index
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    mkdir bowtie
    bowtie-build $args --threads $task.cpus $fasta bowtie/${fasta.baseName}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie: \$(echo \$(bowtie --version 2>&1) | sed 's/^.*bowtie-align-s version //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    """
    mkdir bowtie
    touch bowtie/${fasta.baseName}.1.ebwt
    touch versions.yml
    """
}
