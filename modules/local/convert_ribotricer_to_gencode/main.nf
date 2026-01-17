process CONVERT_RIBOTRICER_TO_GENCODE {
    tag "${meta.id}"
    label 'process_low'

    conda "conda-forge::python=3.9 conda-forge::biopython=1.81 conda-forge::pandas=1.3.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-8849acf39a43cdd6c839a369a74c0adc823e2f91:ab110436faf952a33575c64dd74615a84011450b-0' :
        'biocontainers/mulled-v2-8849acf39a43cdd6c839a369a74c0adc823e2f91:ab110436faf952a33575c64dd74615a84011450b-0' }"

    input:
    tuple val(meta), path(ribotricer_orfs)
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
    ribotricer_to_gencode.py \\
        --tsv ${ribotricer_orfs} \\
        --fasta ${fasta} \\
        --study_id ${study_id} \\
        --output_prefix ${prefix} \\
        --min_length ${min_length} \\
        --min_phase_score 0.5 \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
        biopython: \$(python3 -c "import Bio; print(Bio.__version__)" 2>/dev/null || echo "1.81")
        pandas: \$(python3 -c "import pandas; print(pandas.__version__)" 2>/dev/null || echo "1.3.5")
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
        pandas: "1.3.5"
    END_VERSIONS
    """
}
