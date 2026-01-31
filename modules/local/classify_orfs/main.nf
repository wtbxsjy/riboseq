process CLASSIFY_ORFS_GENCODE {
    def output_dir = (params.orf_classify_output_dir ?: 'orf_classification').tokenize('/').last()

    tag "${output_dir}"
    label 'process_medium'

    publishDir "${params.outdir}/orf_classification", mode: params.publish_dir_mode

    conda "${moduleDir}/environment_gencode.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-8849acf39a43cdd6c839a369a74c0adc823e2f91:ab110436faf952a33575c64dd74615a84011450b-0' :
        'quay.io/biocontainers/mulled-v2-8849acf39a43cdd6c839a369a74c0adc823e2f91:ab110436faf952a33575c64dd74615a84011450b-0' }"

    input:
    path unified_bed
    path unified_metadata
    val input_prefix
    path classify_wrapper
    path class_orf_dir
    path gencode_orf_dir
    path ensembl_dir

    output:
    path "${output_dir}/gencode_results.*", emit: results
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def extra_args = params.extra_orf_classify_args ?: ''
    """
    set -euo pipefail

    mkdir -p ${output_dir}

    python3 ${classify_wrapper} \
        --mode gencode \
        --input ${input_prefix} \
        --output_dir ${output_dir} \
        --ensembl_dir ${ensembl_dir} \
        --cpus ${task.cpus} \
        ${extra_args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p ${output_dir}
    touch ${output_dir}/gencode_results.orfs.gtf
    touch ${output_dir}/gencode_results.orfs.out

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "3.9.0"
    END_VERSIONS
    """
}

process CLASSIFY_ORFS_ORFQUANT {
    def output_dir = (params.orf_classify_output_dir ?: 'orf_classification').tokenize('/').last()

    tag "${output_dir}"
    label 'process_medium'

    publishDir "${params.outdir}/orf_classification", mode: params.publish_dir_mode

    conda "${moduleDir}/environment_orfquant.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/orfquant:1.1.0--r40_1' :
        'quay.io/biocontainers/orfquant:1.1.0--r40_1' }"

    input:
    path unified_gtf
    val input_prefix
    path classify_wrapper
    path class_orf_dir
    path ref_gtf

    output:
    path "${output_dir}/orfquant_classification.tsv", emit: results
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def extra_args = params.extra_orf_classify_args ?: ''
    """
    set -euo pipefail

    mkdir -p ${output_dir}

    python3 ${classify_wrapper} \
        --mode orfquant \
        --input ${input_prefix} \
        --output_dir ${output_dir} \
        --gtf ${ref_gtf} \
        --cpus ${task.cpus} \
        ${extra_args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript --version 2>&1 | sed -n '1p' | awk '{print \$3}')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p ${output_dir}
    touch ${output_dir}/orfquant_classification.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: "4.0"
    END_VERSIONS
    """
}

process CLASSIFY_ORFS_ORF_TYPE {
    def output_dir = (params.orf_classify_output_dir ?: 'orf_classification').tokenize('/').last()

    tag "${output_dir}"
    label 'process_medium'

    publishDir "${params.outdir}/orf_classification", mode: params.publish_dir_mode

    conda "${moduleDir}/environment_orf_type.yml"
    // Use unified container if provided, otherwise use biopython container
    container "${ params.unify_orf_container ?
        params.unify_orf_container :
        (workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
            'https://depot.galaxyproject.org/singularity/biopython:1.79' :
            'quay.io/biocontainers/biopython:1.79') }"

    input:
    path unified_metadata
    val input_prefix
    path classify_wrapper
    path class_orf_dir
    path ref_gtf

    output:
    path "${output_dir}/orftype_classification.tsv", emit: results
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def extra_args = params.extra_orf_classify_args ?: ''
    """
    set -euo pipefail

    mkdir -p ${output_dir}

    python3 ${classify_wrapper} \
        --mode orf_type \
        --input ${input_prefix} \
        --output_dir ${output_dir} \
        --gtf ${ref_gtf} \
        --cpus ${task.cpus} \
        ${extra_args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p ${output_dir}
    touch ${output_dir}/orftype_classification.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "3.9.0"
    END_VERSIONS
    """
}
