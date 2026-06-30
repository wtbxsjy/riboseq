process EXPRESSION_QUANT {
    tag "expression_quant"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        (params.expression_quant_container ?: 'https://depot.galaxyproject.org/singularity/python:3.11') :
        'python:3.11-slim' }"

    input:
    path expression_summary    // from UNIFY (pre-computed per-sample P-site expression)
    path expression_rpkm_tpm   // from UNIFY (pre-computed per-sample RPKM/TPM)
    path orf_confidence        // ORF confidence TSV from ORF_QC (optional)
    path format_script         // format_expression_output.py from bin/

    output:
    path "*_expression_summary.tsv"                 , emit: expression
    path "*_expression_rpkm_tpm.tsv"                , emit: rpkm_tpm
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "expression_quant"
    def min_ocs = params.expression_quant_min_ocs ?: 0.0
    def max_orfs = params.expression_quant_max_orfs ?: 0
    // ORF confidence may not be available (e.g. ORF_QC hasn't completed yet)
    def _orf_conf = orf_confidence instanceof List ? (orf_confidence.isEmpty() ? null : orf_confidence[0]) : (orf_confidence.name != 'NO_FILE' ? orf_confidence : null)
    def orf_conf_opt = _orf_conf ? "--orf-confidence ${_orf_conf}" : ''
    """
    #!/bin/bash
    set -euo pipefail

    echo "=== ORF Expression Quantification (reformat) ==="
    echo "Prefix: ${prefix}"
    echo "Min OCS: ${min_ocs}"

    # Expression stats already computed during UNIFY_ORF_PREDICTIONS
    # This process only reformats and optionally applies OCS filtering.
    python3 ${format_script} \\
        --expression-summary ${expression_summary} \\
        --expression-rpkm-tpm ${expression_rpkm_tpm} \\
        --min-ocs ${min_ocs} \\
        --max-orfs ${max_orfs} \\
        --output-summary "${prefix}_expression_summary.tsv" \\
        --output-rpkm-tpm "${prefix}_expression_rpkm_tpm.tsv" \\
        ${orf_conf_opt}

    echo "--- Done ---"
    ls -lh ${prefix}_*

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | cut -d' ' -f2)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "expression_quant"
    """
    touch ${prefix}_expression_summary.tsv
    touch ${prefix}_expression_rpkm_tpm.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "3.11"
    END_VERSIONS
    """
}
