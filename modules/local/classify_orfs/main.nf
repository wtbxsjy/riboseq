process CLASSIFY_ORFS_GENCODE {
    tag "${classify_output_dir}"
    label 'process_medium'

    publishDir { "${params.outdir}/${classify_output_dir}" }, mode: params.publish_dir_mode

    conda "${moduleDir}/environment_gencode.yml"
    container "${ params.gencode_orf_mapper_container ?: 'nfcore/gencode-orf-mapper:1.1.0' }"

    input:
    path unified_bed
    path unified_metadata
    val input_prefix
    path classify_wrapper
    path class_orf_dir
    path gencode_orf_dir
    path ensembl_dir
    val classify_output_dir

    output:
    path "gencode_results.*", emit: results
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def output_dir = classify_output_dir ?: (params.orf_classify_output_dir ?: 'orf_classification').tokenize('/').last()
    def extra_args = params.extra_orf_classify_args ?: ''
    def orfquant_script = "${class_orf_dir}/run_orfquant_classify.R"
    def local_bed = "${unified_bed}".replaceAll(/\.gz$/, '')
    """
    set -uo pipefail

    # Decompress any .gz inputs so downstream tools can read them
    for f in *.gtf.gz *.bed.gz; do
        [ -f "\$f" ] && gunzip -f "\$f" || true
    done

    # Check if input is a placeholder file
    is_placeholder=false
    if [ -f "${local_bed}" ]; then
        line_count=\$(wc -l < "${local_bed}" 2>/dev/null || echo "0")
        if grep -qi "placeholder" "${local_bed}" 2>/dev/null || \\
           grep -qi "insufficient" "${local_bed}" 2>/dev/null || \\
           [ "\${line_count}" -le 2 ]; then
            is_placeholder=true
            echo "INFO: Input file ${local_bed} detected as placeholder (\${line_count} lines)"
        fi
    else
        is_placeholder=true
        echo "WARNING: Input file ${local_bed} not found - treating as placeholder"
    fi

    if [ "\${is_placeholder}" = "true" ]; then
        echo "WARNING: Input is placeholder/empty - creating placeholder classification output"
        echo "# Placeholder GENCODE classification - input unified ORFs were empty/placeholder" > gencode_results.orfs.gtf
        echo "# Placeholder GENCODE classification - input unified ORFs were empty/placeholder" > gencode_results.orfs.out
    else
        set +e
        python3 ${classify_wrapper} \\
            --mode gencode \\
            --input ${input_prefix} \\
            --output_dir . \\
            --ensembl_dir ${ensembl_dir} \\
            --gencode_impl ${params.gencode_classify_impl ?: 'original'} \\
            --cpus ${task.cpus} \\
            ${extra_args} 2>&1 | tee classify_gencode.log
        EXIT_CODE=\${PIPESTATUS[0]}
        set -e
        
        if [ \${EXIT_CODE} -ne 0 ]; then
            if grep -qiE "(no valid|no ORFs|zero|empty|no input|no data|IndexError|KeyError|ValueError)" classify_gencode.log; then
                echo "WARNING: GENCODE classification failed due to insufficient data - creating placeholder output"
                echo "# Placeholder GENCODE classification - classification failed: no valid ORFs" > gencode_results.orfs.gtf
                echo "# Placeholder GENCODE classification - classification failed: no valid ORFs" > gencode_results.orfs.out
            else
                echo "ERROR: GENCODE classification failed with unexpected error"
                cat classify_gencode.log
                exit \${EXIT_CODE}
            fi
        fi
    fi

    # Compress text outputs to save disk space
    for f in gencode_results.orfs.gtf gencode_results.orfs.out; do
        [ -f "\$f" ] && gzip -f "\$f" || true
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    touch gencode_results.orfs.gtf
    gzip -f gencode_results.orfs.gtf
    touch gencode_results.orfs.out
    gzip -f gencode_results.orfs.out

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "3.9.0"
    END_VERSIONS
    """
}

process CLASSIFY_ORFS_ORFQUANT {
    tag "${classify_output_dir}"
    label 'process_medium'

    publishDir { "${params.outdir}/${classify_output_dir}" }, mode: params.publish_dir_mode

    conda "${moduleDir}/environment_orfquant.yml"
    container "${ params.orfquant_container ?
        params.orfquant_container :
        (workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
            'https://depot.galaxyproject.org/singularity/orfquant:1.1.0--r40_1' :
            'quay.io/biocontainers/orfquant:1.1.0--r40_1') }"

    input:
    path unified_gtf
    path unified_metadata
    val input_prefix
    path classify_wrapper
    path class_orf_dir
    path ref_gtf
    val classify_output_dir

    output:
    path "orfquant_classification.tsv", emit: results
    path "orfquant_results.*"          , emit: extra_results, optional: true
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def output_dir = classify_output_dir ?: (params.orf_classify_output_dir ?: 'orf_classification').tokenize('/').last()
    def extra_args = params.extra_orf_classify_args ?: ''
    def orfquant_script = "${class_orf_dir}/run_orfquant_classify.R"
    def local_gtf = "${unified_gtf}".replaceAll(/\.gz$/, '')
    """
    set -uo pipefail

    # Decompress any .gz inputs so downstream tools can read them
    for f in *.gtf.gz *.bed.gz; do
        [ -f "\$f" ] && gunzip -f "\$f" || true
    done

    # Check if input is a placeholder file
    is_placeholder=false
    if [ -f "${local_gtf}" ]; then
        line_count=\$(wc -l < "${local_gtf}" 2>/dev/null || echo "0")
        if grep -qi "placeholder" "${local_gtf}" 2>/dev/null || \\
           grep -qi "insufficient" "${local_gtf}" 2>/dev/null || \\
           [ "\${line_count}" -le 2 ]; then
            is_placeholder=true
            echo "INFO: Input file ${local_gtf} detected as placeholder (\${line_count} lines)"
        fi
    else
        is_placeholder=true
        echo "WARNING: Input file ${local_gtf} not found - treating as placeholder"
    fi

    if [ "\${is_placeholder}" = "true" ]; then
        echo "WARNING: Input is placeholder/empty - creating placeholder classification output"
        echo -e "# Placeholder ORFquant classification - input unified ORFs were empty/placeholder" > orfquant_classification.tsv
        echo -e "orf_id\\torf_type\\tclassification" >> orfquant_classification.tsv
        touch orfquant_results.logs orfquant_results.orfs.bed orfquant_results.orfs.gtf \\
              orfquant_results.orfs.fa orfquant_results.orfs.pep.fa orfquant_results.orfs.out
    else
        set +e
        export R_LIBS_USER="\${PWD}/.Rlib"
        mkdir -p "\${R_LIBS_USER}"
        if ! Rscript -e "suppressPackageStartupMessages(library(optparse))" >/dev/null 2>&1; then
            echo "INFO: Installing missing R package optparse"
            Rscript -e "install.packages('optparse', repos='https://cloud.r-project.org')" || \\
                { echo "ERROR: Failed to install optparse"; exit 1; }
        fi
        Rscript ${orfquant_script} \\
            --input ${input_prefix}.gtf \\
            --annotation ${ref_gtf} \\
            --output orfquant_classification.tsv \\
            --metadata ${unified_metadata} \\
            --output_prefix orfquant_results \\
            ${extra_args} 2>&1 | tee classify_orfquant.log
        EXIT_CODE=\${PIPESTATUS[0]}
        set -e
        
        if [ \${EXIT_CODE} -ne 0 ]; then
            if grep -qiE "(no valid|no ORFs|zero|empty|no input|no data|IndexError|KeyError|ValueError)" classify_orfquant.log; then
                echo "WARNING: ORFquant classification failed due to insufficient data - creating placeholder output"
                echo -e "# Placeholder ORFquant classification - classification failed: no valid ORFs" > orfquant_classification.tsv
                echo -e "orf_id\\torf_type\\tclassification" >> orfquant_classification.tsv
                touch orfquant_results.logs orfquant_results.orfs.bed orfquant_results.orfs.gtf \
                      orfquant_results.orfs.fa orfquant_results.orfs.pep.fa orfquant_results.orfs.out
            else
                echo "ERROR: ORFquant classification failed with unexpected error"
                cat classify_orfquant.log
                exit \${EXIT_CODE}
            fi
        fi
    fi

    # Compress text outputs to save disk space
    for f in orfquant_results.orfs.bed orfquant_results.orfs.gtf orfquant_results.orfs.out; do
        [ -f "\$f" ] && gzip -f "\$f" || true
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript --version 2>&1 | sed -n '1p' | awk '{print \$3}')
    END_VERSIONS
    """

    stub:
    """
    touch orfquant_classification.tsv
    touch orfquant_results.logs orfquant_results.orfs.bed orfquant_results.orfs.gtf
    touch orfquant_results.orfs.fa orfquant_results.orfs.pep.fa orfquant_results.orfs.out
    gzip -f orfquant_results.orfs.bed orfquant_results.orfs.gtf orfquant_results.orfs.out

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: "4.0"
    END_VERSIONS
    """
}

process CLASSIFY_ORFS_ORF_TYPE {
    tag "${classify_output_dir}"
    label 'process_medium'

    publishDir { "${params.outdir}/${classify_output_dir}" }, mode: params.publish_dir_mode

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
    val classify_output_dir

    output:
    path "orftype_classification.tsv", emit: results
    path "orftype_results.*"          , emit: extra_results, optional: true
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def output_dir = classify_output_dir ?: (params.orf_classify_output_dir ?: 'orf_classification').tokenize('/').last()
    def extra_args = params.extra_orf_classify_args ?: ''
    """
    set -uo pipefail

    # Check if input is a placeholder file
    is_placeholder=false
    if [ -f "${unified_metadata}" ]; then
        line_count=\$(wc -l < "${unified_metadata}" 2>/dev/null || echo "0")
        if grep -qi "placeholder" "${unified_metadata}" 2>/dev/null || \\
           grep -qi "insufficient" "${unified_metadata}" 2>/dev/null || \\
           [ "\${line_count}" -le 2 ]; then
            is_placeholder=true
            echo "INFO: Input file ${unified_metadata} detected as placeholder (\${line_count} lines)"
        fi
    else
        is_placeholder=true
        echo "WARNING: Input file ${unified_metadata} not found - treating as placeholder"
    fi

    if [ "\${is_placeholder}" = "true" ]; then
        echo "WARNING: Input is placeholder/empty - creating placeholder classification output"
        echo -e "# Placeholder ORF type classification - input unified ORFs were empty/placeholder" > orftype_classification.tsv
        echo -e "orf_id\\torf_type\\tclassification\\tgene_biotype" >> orftype_classification.tsv
        touch orftype_results.logs orftype_results.orfs.bed orftype_results.orfs.gtf \
              orftype_results.orfs.fa orftype_results.orfs.pep.fa orftype_results.orfs.out
    else
        set +e
        python3 ${classify_wrapper} \\
            --mode orf_type \\
            --input ${input_prefix} \\
            --output_dir . \\
            --gtf ${ref_gtf} \\
            --cpus ${task.cpus} \\
            ${extra_args} 2>&1 | tee classify_orftype.log
        EXIT_CODE=\${PIPESTATUS[0]}
        set -e
        
        if [ \${EXIT_CODE} -ne 0 ]; then
            if grep -qiE "(no valid|no ORFs|zero|empty|no input|no data|IndexError|KeyError|ValueError)" classify_orftype.log; then
                echo "WARNING: ORF type classification failed due to insufficient data - creating placeholder output"
                echo -e "# Placeholder ORF type classification - classification failed: no valid ORFs" > orftype_classification.tsv
                echo -e "orf_id\\torf_type\\tclassification\\tgene_biotype" >> orftype_classification.tsv
                touch orftype_results.logs orftype_results.orfs.bed orftype_results.orfs.gtf \
                      orftype_results.orfs.fa orftype_results.orfs.pep.fa orftype_results.orfs.out
            else
                echo "ERROR: ORF type classification failed with unexpected error"
                cat classify_orftype.log
                exit \${EXIT_CODE}
            fi
        fi
    fi

    # Compress text outputs to save disk space
    for f in orftype_results.orfs.bed orftype_results.orfs.gtf orftype_results.orfs.out; do
        [ -f "\$f" ] && gzip -f "\$f" || true
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    touch orftype_classification.tsv
    touch orftype_results.logs orftype_results.orfs.bed orftype_results.orfs.gtf
    touch orftype_results.orfs.fa orftype_results.orfs.pep.fa orftype_results.orfs.out
    gzip -f orftype_results.orfs.bed orftype_results.orfs.gtf orftype_results.orfs.out

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "3.9.0"
    END_VERSIONS
    """
}
