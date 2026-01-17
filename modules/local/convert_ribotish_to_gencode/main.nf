process CONVERT_RIBOTISH_TO_GENCODE {
    tag "${meta.id}"
    label 'process_low'

    conda "conda-forge::python=3.9 conda-forge::biopython=1.81"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.81' :
        'biocontainers/biopython:1.81' }"

    input:
    tuple val(meta), path(ribotish_predict), path(ribotish_quality)
    path fasta
    path gtf

    output:
    tuple val(meta), path("*.gencode.fa"), emit: fasta
    tuple val(meta), path("*.gencode.bed"), emit: bed
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def study_id = meta.study_id ?: meta.id
    def min_length = params.gencode_orf_min_length ?: 16

    """
    ribotish_to_gencode.py \\
        --predict ${ribotish_predict} \\
        --fasta ${fasta} \\
        --study_id ${study_id} \\
        --output_prefix ${prefix} \\
        --min_length ${min_length} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
        biopython: \$(python3 -c "import Bio; print(Bio.__version__)" 2>/dev/null || echo "1.81")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.gencode.fa
    touch ${prefix}.gencode.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "3.9.0"
        biopython: "1.81"
    END_VERSIONS
    """
}
