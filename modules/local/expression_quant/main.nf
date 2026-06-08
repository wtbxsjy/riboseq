process EXPRESSION_QUANT {
    tag "$prefix"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        (params.expression_quant_container ?: 'https://depot.galaxyproject.org/singularity/python:3.11') :
        'python:3.11-slim' }"

    input:
    path unified_meta                               // unified ORF metadata
    path unified_bed                                // unified ORF BED
    path orf_confidence                             // ORF confidence TSV (from ORF_QC)
    tuple val(sample_list), path(psites_bedgraph)   // all P-site bedgraphs from RiboseQC
    tuple val(sample_list_cov), path(coverage_bedgraph) // all coverage bedgraphs

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
    def workers = params.expression_quant_workers ?: 4
    """
    #!/bin/bash
    set -euo pipefail

    # Install dependencies if needed
    pip install --quiet pandas numpy 2>/dev/null || true

    echo "=== ORF Expression Quantification ==="
    echo "Prefix: ${prefix}"
    echo "Min OCS: ${min_ocs}"
    echo "Workers: ${workers}"

    # Collect sample names from bedgraph files
    ls *_P_sites_plus.bedgraph 2>/dev/null | sed 's/_P_sites_plus.bedgraph//' > sample_list.txt || true
    n_samples=\$(wc -l < sample_list.txt || echo 0)
    echo "Samples with P-site data: \$n_samples"

    if [ "\$n_samples" -eq 0 ]; then
        echo "WARNING: No P-site bedgraph files found. Creating placeholder output."
        echo -e "orf_id\\tchrom\\tstart\\tend\\tstrand\\ttotal_reads\\tn_expressed_samples" > "${prefix}_expression_summary.tsv"
        echo -e "orf_id\\tchrom\\tstart\\tend\\tstrand\\torf_length\\torf_length_kb" > "${prefix}_expression_rpkm_tpm.tsv"
    else
        # Phase 1: P-site expression quantification
        echo "--- Phase 1: P-site expression quantification ---"
        quantify_orf_expression.py \\
            --orf-meta ${unified_meta} \\
            --orf-confidence ${orf_confidence} \\
            --psites-dir . \\
            --sample-pattern "*_P_sites_plus.bedgraph" \\
            --output "${prefix}_expression_summary.tsv" \\
            --min-ocs ${min_ocs} \\
            --max-orfs ${max_orfs} \\
            --workers ${workers}

        # Phase 2: RPKM/TPM calculation
        echo "--- Phase 2: RPKM/TPM calculation ---"
        calc_orf_rpkm_tpm.py \\
            --expression "${prefix}_expression_summary.tsv" \\
            --coverage-dir . \\
            --sample-pattern "*_coverage_plus.bedgraph" \\
            --output "${prefix}_expression_rpkm_tpm.tsv" \\
            --workers ${workers}

        echo "--- Done ---"
        ls -lh ${prefix}_*
    fi

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
