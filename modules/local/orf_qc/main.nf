process ORF_QC {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        (params.orf_qc_container ?: 'https://depot.galaxyproject.org/singularity/python:3.11') :
        'python:3.11-slim' }"

    input:
    // Required: unified ORFs (post-UNIFY_ORF_PREDICTIONS)
    tuple val(meta), path(unified_bed), path(unified_meta)

    // All ORF prediction tool outputs — each optional
    tuple val(meta2), path(ribocode_collapsed)       // RiboCode _collapsed.txt
    tuple val(meta3), path(riboseqc_psites)           // RiboseQC _P_sites_calcs
    tuple val(meta4), path(ribowaltz_offsets)         // riboWaltz _psite_offset.tsv
    tuple val(meta5), path(ribowaltz_frames)          // riboWaltz _frame_distribution.tsv
    tuple val(meta6), path(ribotricer_tsv)            // Ribotricer _translating_ORFs.tsv
    tuple val(meta7), path(ribotish_pred)             // Ribo-TISH _pred.txt
    tuple val(meta8), path(ribotish_para)             // Ribo-TISH .para.py
    tuple val(meta9), path(price_tsv)                 // PRICE .orfs.tsv
    tuple val(meta10), path(rpbp_bayes)               // rp-bp bayes-factors.bed.gz
    tuple val(meta11), path(orfquant_gtf)             // ORFquant _Detected_ORFs.gtf.gz

    output:
    tuple val(meta), path("${prefix}_qc_report.html")    , emit: report
    tuple val(meta), path("${prefix}_orf_confidence.tsv"), emit: confidence
    tuple val(meta), path("${prefix}_tool_agreement.tsv"), emit: agreement
    tuple val(meta), path("${prefix}_psite_harmonized.tsv"), emit: psite_harmonized
    tuple val(meta), path("${prefix}_sample_flags.json") , emit: flags
    tuple val(meta), path("${prefix}_qc_metrics.tsv")    , emit: metrics
    path "versions.yml"                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    """
    #!/bin/bash
    set -euo pipefail

    echo "=== ORF QC Module ==="
    echo "Sample: ${prefix}"

    # Collect all available tool output files into a list
    # We use a file list approach so the QC scripts can auto-detect
    rm -f input_files.txt

    for f in \\
        ${ribocode_collapsed} \\
        ${riboseqc_psites} \\
        ${ribowaltz_offsets} \\
        ${ribowaltz_frames} \\
        ${ribotricer_tsv} \\
        ${ribotish_pred} \\
        ${ribotish_para} \\
        ${price_tsv} \\
        ${rpbp_bayes} \\
        ${orfquant_gtf} \\
        ; do
        if [ -n "\$f" ] && [ "\$f" != "null" ] && [ -f "\$f" ]; then
            echo "\$f" >> input_files.txt
        fi
    done

    n_files=\$(wc -l < input_files.txt || echo 0)
    echo "Available tool output files: \$n_files"

    if [ "\$n_files" -eq 0 ]; then
        echo "WARNING: No tool output files found. Creating placeholder QC outputs."
        echo -e "orf_id\\tgene_id\\tgene_name\\tchrom\\torf_type\\tocs\\ttier\\ts_translation\\ts_agreement\\ts_coverage\\ts_periodicity\\ts_readlevel\\tdetecting_tools\\tn_detecting" > "${prefix}_orf_confidence.tsv"
        echo -e "tool_a\\ttool_b\\tjaccard\\toverlap_count\\ta_total\\tb_total" > "${prefix}_tool_agreement.tsv"
        echo -e "read_length\\tn_tools\\tconsensus_offset\\tmax_delta\\tflag" > "${prefix}_psite_harmonized.tsv"
        echo '{"flags":[],"summary":{"error":"No tool output files found"}}' > "${prefix}_sample_flags.json"
        echo -e "metric\\tvalue\\nsample_id\\t${prefix}" > "${prefix}_qc_metrics.tsv"

        # Generate minimal placeholder HTML report
        cat > "${prefix}_qc_report.html" <<'EOF'
    <!DOCTYPE html>
    <html><head><title>ORF QC — No Data</title></head>
    <body style="font-family:sans-serif;padding:40px;text-align:center">
    <h1>ORF QC Report</h1>
    <p style="color:#dc3545">⚠ No ORF prediction tool outputs were found.</p>
    <p>All tools may have been skipped or failed. Check pipeline logs for details.</p>
    </body></html>
    EOF
    else
        # Install Python dependencies if needed
        pip install --quiet pandas numpy scipy plotly 2>/dev/null || true

        # Phase 1: Extract metrics from all tools
        echo "--- Phase 1: Extracting metrics ---"
        extract_orf_qc_metrics.py \\
            --file-list input_files.txt \\
            --output tool_data.json

        # Phase 2: Harmonize P-site offsets and periodicity
        echo "--- Phase 2: Harmonizing metrics ---"
        harmonize_orf_qc.py \\
            --input tool_data.json \\
            --output-prefix "${prefix}"

        # Phase 3+4: Cross-tool comparison + OCS scoring
        echo "--- Phase 3+4: Cross-tool comparison and confidence scoring ---"
        compare_orf_tools.py \\
            --tool-data tool_data.json \\
            --unified-meta ${unified_meta} \\
            --periodicity "${prefix}_periodicity.json" \\
            --output-prefix "${prefix}"

        # Phase 5: Generate HTML report
        echo "--- Phase 5: Generating report ---"
        generate_orf_qc_report.py \\
            --psite "${prefix}_psite_harmonized.tsv" \\
            --confidence "${prefix}_orf_confidence.tsv" \\
            --periodicity "${prefix}_periodicity.json" \\
            --agreement "${prefix}_tool_agreement.json" \\
            --flags "${prefix}_sample_flags.json" \\
            --sample-id "${prefix}" \\
            --output "${prefix}_qc_report.html"
    fi

    # Write versions
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | cut -d' ' -f2)
        bedtools: \$(bedtools --version 2>&1 | cut -d' ' -f2 || echo "unknown")
    END_VERSIONS

    echo "=== ORF QC Module complete ==="
    ls -lh ${prefix}_*
    """
}
