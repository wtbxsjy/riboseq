process UNIFY_ORF_PREDICTIONS {
    def prefix = (params.unify_orf_predictions_prefix ?: 'unified_orfs').tokenize('/').last()

    tag "${prefix}"
    label 'process_medium'

    publishDir "${params.outdir}/orf_unification", mode: params.publish_dir_mode

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9' :
        'python:3.9' }"

    input:
    tuple val(ribotish_files), val(ribotricer_files), val(orfquant_files), path(all_inputs)
    path gtf
    path fasta
    path unify_script

    output:
    path "${prefix}.metadata.tsv", emit: metadata
    path "${prefix}.bed"         , emit: bed
    path "${prefix}.gtf"         , emit: gtf
    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def min_len = params.unify_orf_min_len ?: 10
    def extra_args = params.extra_unify_orf_predictions_args ?: ''
    def ribotish_arg = (ribotish_files && ribotish_files.size() > 0) ? "--ribotish ${ribotish_files.collect{ it.getFileName().toString() }.join(' ')}" : ''
    def ribotricer_arg = (ribotricer_files && ribotricer_files.size() > 0) ? "--ribotricer ${ribotricer_files.collect{ it.getFileName().toString() }.join(' ')}" : ''
    def orfquant_arg = (orfquant_files && orfquant_files.size() > 0) ? "--orfquant ${orfquant_files.collect{ it.getFileName().toString() }.join(' ')}" : ''
    """
    set -euo pipefail

    export PYTHONUSERBASE="$PWD/.pylibs"
    export PATH="$PYTHONUSERBASE/bin:$PATH"

    pip install --user --quiet --no-warn-script-location biopython pyfaidx

    python3 ${unify_script} \
        --gtf ${gtf} \
        --fasta ${fasta} \
        --output ${prefix} \
        --min_len ${min_len} \
        ${ribotish_arg} \
        ${ribotricer_arg} \
        ${orfquant_arg} \
        ${extra_args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
        biopython: \$(python3 -c "import Bio; print(Bio.__version__)" 2>/dev/null || echo "unknown")
        pyfaidx: \$(python3 -c "import pyfaidx; print(pyfaidx.__version__)" 2>/dev/null || echo "unknown")
    END_VERSIONS
    """

    stub:
    """
    touch ${prefix}.metadata.tsv
    touch ${prefix}.bed
    touch ${prefix}.gtf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "3.9.0"
        biopython: "1.81"
        pyfaidx: "0.7"
    END_VERSIONS
    """
}
