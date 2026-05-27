process UNIFY_ORF_PREDICTIONS {
    tag "${(params.unify_orf_predictions_prefix ?: 'unified_orfs').tokenize('/').last()}"
    label 'process_medium'

    publishDir "${params.outdir}/orf_unification", mode: params.publish_dir_mode

    conda "${moduleDir}/environment.yml"
    // Use custom container if provided, otherwise use biopython container
    container "${ params.unify_orf_container ?
        params.unify_orf_container :
        (workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
            'https://depot.galaxyproject.org/singularity/biopython:1.79' :
            'quay.io/biocontainers/biopython:1.79') }"

    input:
    tuple val(ribotish_files), val(ribotricer_files), val(ribocode_files), val(orfquant_files), path(all_inputs)
    path gtf
    path fasta
    path unify_script
    path run_orf_script
    path psites_bedgraph, stageAs: 'bedgraph/*'  // RiboseQC P-site bedgraph files (optional)
    val sample_list       // List of sample names for bedgraph files

    output:
    path "*.metadata.tsv", emit: metadata
    path "*.bed.gz"      , emit: bed
    path "*.gtf.gz"      , emit: gtf
    path "*.stats.txt"   , emit: stats
    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = (params.unify_orf_predictions_prefix ?: 'unified_orfs').tokenize('/').last()
    def min_len = params.unify_orf_min_len ?: 6
    def extra_args = params.extra_unify_orf_predictions_args ?: ''
    // Advanced merging parameters
    def frame_merge_min_overlap = params.unify_orf_frame_merge_min_overlap ?: 0.9
    def no_frame_merge = params.unify_orf_no_frame_merge ? "--no-frame-merge" : ''
    def seq_cluster = params.unify_orf_seq_cluster ? "--seq-cluster" : ''
    // Files are already staged via all_inputs, just need filenames properly quoted
    def ribotish_arg = (ribotish_files && ribotish_files instanceof List && ribotish_files.size() > 0) ?
        "--ribotish ${ribotish_files.collect{ "\"${it}\"" }.join(' ')}" : ''
    def ribotricer_arg = (ribotricer_files && ribotricer_files instanceof List && ribotricer_files.size() > 0) ?
        "--ribotricer ${ribotricer_files.collect{ "\"${it}\"" }.join(' ')}" : ''
    def ribocode_arg = (ribocode_files && ribocode_files instanceof List && ribocode_files.size() > 0) ?
        "--ribocode ${ribocode_files.collect{ "\"${it}\"" }.join(' ')}" : ''
    def orfquant_arg = (orfquant_files && orfquant_files instanceof List && orfquant_files.size() > 0) ?
        "--orfquant ${orfquant_files.collect{ "\"${it}\"" }.join(' ')}" : ''
    // Bedgraph arguments for P-site statistics from RiboseQC
    // psites_bedgraph can be a list or path, check if it has files
    def has_bedgraph = (psites_bedgraph instanceof List && psites_bedgraph.size() > 0) ||
                       (psites_bedgraph && !(psites_bedgraph instanceof List) && psites_bedgraph.name != 'NO_FILE')
    def bedgraph_arg = has_bedgraph ? "--bedgraph-dir bedgraph" : ''
    def sample_arg = (sample_list && sample_list instanceof List && sample_list.size() > 0) ? "--sample-list ${sample_list.join(',')}" : ''

    // Generate list of input files for bash to iterate over
    def input_files_list = (all_inputs instanceof List) ? all_inputs.collect{ it.name }.join(' ') : (all_inputs ? all_inputs.name : '')
    """
    set -uo pipefail

    # Check if dependencies are already available (custom container)
    if python3 -c "import Bio; import pyfaidx" 2>/dev/null; then
        echo "Dependencies already available in container"
    else
        # Setup user-local Python package directory to avoid permission issues
        export PYTHONUSERBASE="\$PWD/.pylibs"
        export PATH="\$PYTHONUSERBASE/bin:\$PATH"
        export PYTHONPATH="\$PYTHONUSERBASE/lib/python3.9/site-packages:\${PYTHONPATH:-}"
        export PIP_NO_CACHE_DIR=1
        mkdir -p "\$PYTHONUSERBASE"

        echo "Installing Python dependencies..."
        pip install --user --no-cache-dir --no-warn-script-location pyfaidx biopython 2>&1 || {
            echo "pip install failed, trying with python -m pip..."
            python3 -m pip install --user --no-cache-dir pyfaidx biopython
        }
    fi
    
    # Verify installation
    echo "Verifying dependencies..."
    python3 -c "import Bio; import pyfaidx; print('Dependencies OK: Bio=' + Bio.__version__ + ', pyfaidx=' + pyfaidx.__version__)"

    # Optionally install orfont extras (pandas + duckdb) for optimized path
    if [ -f "${run_orf_script}" ]; then
        if ! python3 -c "import pandas; import duckdb" 2>/dev/null; then
            echo "Installing orfont extras (pandas + duckdb)..."
            pip install --user --no-cache-dir --no-warn-script-location pandas duckdb 2>&1 || \
                python3 -m pip install --user --no-cache-dir pandas duckdb 2>&1 || \
                echo "WARNING: Could not install pandas/duckdb - falling back to original path"
        fi
    fi

    # Decompress any .gz input files so validation tools can read them
    for f in *.gtf.gz *.bed.gz; do
        [ -f "\$f" ] && gunzip -f "\$f" || true
    done

    # Check if all input files are empty/placeholder files
    # Input files list generated from Nextflow: ${input_files_list}
    has_valid_input=false
    input_files_to_check="${input_files_list}"

    if [ -n "\${input_files_to_check}" ]; then
        for input_file in \${input_files_to_check}; do
            check_file="\${input_file}"
            # After decompression, the .gz file may be gone; check uncompressed name
            [ -f "\${check_file}" ] || check_file="\${input_file%.gz}"
            if [ -f "\${check_file}" ]; then
                # Check if file has actual content (not just header or placeholder)
                line_count=\$(wc -l < "\${check_file}" 2>/dev/null || echo "0")
                if [ "\${line_count}" -gt 2 ]; then
                    # Check it's not a placeholder file
                    if ! grep -qi "placeholder" "\${check_file}" 2>/dev/null && \\
                       ! grep -qi "# Empty" "\${check_file}" 2>/dev/null && \\
                       ! grep -qi "insufficient" "\${check_file}" 2>/dev/null; then
                        has_valid_input=true
                        echo "Found valid input file: \${check_file} (\${line_count} lines)"
                        break
                    else
                        echo "Skipping placeholder file: \${check_file}"
                    fi
                else
                    echo "Skipping empty/small file: \${check_file} (\${line_count} lines)"
                fi
            fi
        done
    else
        echo "WARNING: No input files provided"
    fi

    if [ "\${has_valid_input}" = "false" ]; then
        echo "WARNING: All input files are empty or placeholder files - creating placeholder unified output"
        
        # Create placeholder metadata
        echo "# Placeholder unified ORF metadata - no valid input data from ORF prediction tools" > ${prefix}.metadata.tsv
        echo -e "orf_id\\torf_name\\torf_type\\tchrom\\tstrand\\tstart\\tend\\tsource\\tgene_id\\ttranscript_id\\torf_length\\tscore\\tsource_count" >> ${prefix}.metadata.tsv
        
        # Create placeholder BED file
        echo "# Placeholder unified ORF BED - no valid input data from ORF prediction tools" > ${prefix}.bed
        echo -e "#chrom\\tchromStart\\tchromEnd\\tname\\tscore\\tstrand\\tthickStart\\tthickEnd\\titemRgb\\tblockCount\\tblockSizes\\tblockStarts" >> ${prefix}.bed
        
        # Create placeholder GTF file  
        echo "# Placeholder unified ORF GTF - no valid input data from ORF prediction tools" > ${prefix}.gtf
        
        # Create stats file
        {
            echo "=== ORF Unification Statistics ==="
            echo "WARNING: No valid input files found"
            echo "All input files were empty or placeholder files"
            echo "No ORF unification was performed"
        } > ${prefix}.stats.txt
        
        echo "WARNING: UNIFY_ORF_PREDICTIONS created placeholder outputs due to insufficient input data"
    else
        set +e

        # Choose optimized (orfont) or original path
        if [ -f "${run_orf_script}" ] && python3 -c "import pandas; import duckdb" 2>/dev/null; then
            echo "Using optimized orfont path (run_orf.py)"
            python3 ${run_orf_script} unify \\
                --gtf ${gtf} \\
                --fasta ${fasta} \\
                --output ${prefix} \\
                --min_len ${min_len} \\
                --threads ${task.cpus} \\
                --frame-merge-min-overlap ${frame_merge_min_overlap} \\
                ${no_frame_merge} \\
                ${seq_cluster} \\
                ${ribotish_arg} \\
                ${ribotricer_arg} \\
                ${ribocode_arg} \\
                ${orfquant_arg} \\
                ${bedgraph_arg} \\
                ${sample_arg} \\
                ${extra_args} 2>&1 | tee unify_orf.log
        else
            echo "Using original unify path"
            python3 ${unify_script} \\
                --gtf ${gtf} \\
                --fasta ${fasta} \\
                --output ${prefix} \\
                --min_len ${min_len} \\
                --threads ${task.cpus} \\
                --frame-merge-min-overlap ${frame_merge_min_overlap} \\
                ${no_frame_merge} \\
                ${seq_cluster} \\
                ${ribotish_arg} \\
                ${ribotricer_arg} \\
                ${ribocode_arg} \\
                ${orfquant_arg} \\
                ${bedgraph_arg} \\
                ${sample_arg} \\
                ${extra_args} 2>&1 | tee unify_orf.log
        fi
        EXIT_CODE=\${PIPESTATUS[0]}
        set -e

        # Extract statistics from the log and save to stats file
        {
            echo "=== ORF Unification Statistics ==="
            grep -E "^===|^By Tool:|^By Sample:|^Final|^After|^Total|^Note:|^  (ribotish|ribotricer|ribocode|orfquant|RiboCode|[A-Za-z0-9_])" unify_orf.log || true
        } > ${prefix}.stats.txt

        if [ \${EXIT_CODE} -ne 0 ]; then
            # Check for recoverable errors (no valid ORFs found, empty results, etc.)
            if grep -qiE "(no valid|no ORFs|zero ORFs|empty|no predictions|0 ORFs)" unify_orf.log || \\
               grep -qiE "(ValueError|KeyError|IndexError).*empty" unify_orf.log || \\
               grep -qiE "(cannot|failed to).*merge|parse" unify_orf.log; then
                echo "WARNING: UNIFY_ORF_PREDICTIONS failed due to insufficient/invalid ORF data - creating placeholder files"

                # Create placeholder metadata
                echo "# Placeholder unified ORF metadata - unification failed: insufficient valid ORF predictions" > ${prefix}.metadata.tsv
                echo -e "orf_id\\torf_name\\torf_type\\tchrom\\tstrand\\tstart\\tend\\tsource\\tgene_id\\ttranscript_id\\torf_length\\tscore\\tsource_count" >> ${prefix}.metadata.tsv

                # Create placeholder BED file
                echo "# Placeholder unified ORF BED - unification failed: insufficient valid ORF predictions" > ${prefix}.bed
                echo -e "#chrom\\tchromStart\\tchromEnd\\tname\\tscore\\tstrand\\tthickStart\\tthickEnd\\titemRgb\\tblockCount\\tblockSizes\\tblockStarts" >> ${prefix}.bed

                # Create placeholder GTF file
                echo "# Placeholder unified ORF GTF - unification failed: insufficient valid ORF predictions" > ${prefix}.gtf

                # Add error note to stats
                echo "ERROR: ORF unification failed due to insufficient valid ORF data" >> ${prefix}.stats.txt
            else
                echo "ERROR: UNIFY_ORF_PREDICTIONS failed with unexpected error"
                cat unify_orf.log
                exit \${EXIT_CODE}
            fi
        fi
    fi

    # Compress text outputs to save disk space
    gzip -f ${prefix}.bed ${prefix}.gtf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
        biopython: \$(python3 -c "import Bio; print(Bio.__version__)" 2>/dev/null || echo "unknown")
        pyfaidx: \$(python3 -c "import pyfaidx; print(pyfaidx.__version__)" 2>/dev/null || echo "unknown")
        pandas: \$(python3 -c "import pandas; print(pandas.__version__)" 2>/dev/null || echo "N/A")
        duckdb: \$(python3 -c "import duckdb; print(duckdb.__version__)" 2>/dev/null || echo "N/A")
    END_VERSIONS
    """

    stub:
    def prefix = (params.unify_orf_predictions_prefix ?: 'unified_orfs').tokenize('/').last()
    """
    touch ${prefix}.metadata.tsv
    touch ${prefix}.bed
    gzip -f ${prefix}.bed
    touch ${prefix}.gtf
    gzip -f ${prefix}.gtf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "3.9.0"
        biopython: "1.81"
        pyfaidx: "0.7"
        pandas: "N/A"
        duckdb: "N/A"
    END_VERSIONS
    """
}

// Per-tool variant: exact-match dedup only, no cross-tool frame-aware merging.
// Outputs per-tool files for each tool that ran, plus a combined exact-dedup set.
// Use this when skip_unify_orf_predictions = true to enable per-tool classification.
process UNIFY_ORF_PREDICTIONS_PER_TOOL {
    tag "${(params.unify_orf_predictions_prefix ?: 'unified_orfs').tokenize('/').last()}"
    label 'process_medium'

    publishDir "${params.outdir}/orf_unification/per_tool", mode: params.publish_dir_mode

    conda "${moduleDir}/environment.yml"
    container "${ params.unify_orf_container ?
        params.unify_orf_container :
        (workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
            'https://depot.galaxyproject.org/singularity/biopython:1.79' :
            'quay.io/biocontainers/biopython:1.79') }"

    input:
    tuple val(ribotish_files), val(ribotricer_files), val(ribocode_files), val(orfquant_files), path(all_inputs)
    path gtf
    path fasta
    path unify_script
    path run_orf_script
    path psites_bedgraph, stageAs: 'bedgraph/*'
    val sample_list

    output:
    path "*_ribotish.metadata.tsv"  , emit: ribotish_metadata  , optional: true
    path "*_ribotish.bed.gz"        , emit: ribotish_bed       , optional: true
    path "*_ribotish.gtf.gz"        , emit: ribotish_gtf       , optional: true
    path "*_ribotricer.metadata.tsv", emit: ribotricer_metadata, optional: true
    path "*_ribotricer.bed.gz"      , emit: ribotricer_bed     , optional: true
    path "*_ribotricer.gtf.gz"      , emit: ribotricer_gtf     , optional: true
    path "*_ribocode.metadata.tsv"  , emit: ribocode_metadata  , optional: true
    path "*_ribocode.bed.gz"        , emit: ribocode_bed       , optional: true
    path "*_ribocode.gtf.gz"        , emit: ribocode_gtf       , optional: true
    path "*_orfquant.metadata.tsv"  , emit: orfquant_metadata  , optional: true
    path "*_orfquant.bed.gz"        , emit: orfquant_bed       , optional: true
    path "*_orfquant.gtf.gz"        , emit: orfquant_gtf       , optional: true
    path "${(params.unify_orf_predictions_prefix ?: 'unified_orfs').tokenize('/').last()}.metadata.tsv", emit: combined_metadata
    path "${(params.unify_orf_predictions_prefix ?: 'unified_orfs').tokenize('/').last()}.bed.gz"      , emit: combined_bed
    path "${(params.unify_orf_predictions_prefix ?: 'unified_orfs').tokenize('/').last()}.gtf.gz"      , emit: combined_gtf
    path "${(params.unify_orf_predictions_prefix ?: 'unified_orfs').tokenize('/').last()}.stats.txt"   , emit: stats
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = (params.unify_orf_predictions_prefix ?: 'unified_orfs').tokenize('/').last()
    def min_len = params.unify_orf_min_len ?: 6
    def extra_args = params.extra_unify_orf_predictions_args ?: ''
    def ribotish_arg = (ribotish_files && ribotish_files instanceof List && ribotish_files.size() > 0) ?
        "--ribotish ${ribotish_files.collect{ "\"${it}\"" }.join(' ')}" : ''
    def ribotricer_arg = (ribotricer_files && ribotricer_files instanceof List && ribotricer_files.size() > 0) ?
        "--ribotricer ${ribotricer_files.collect{ "\"${it}\"" }.join(' ')}" : ''
    def ribocode_arg = (ribocode_files && ribocode_files instanceof List && ribocode_files.size() > 0) ?
        "--ribocode ${ribocode_files.collect{ "\"${it}\"" }.join(' ')}" : ''
    def orfquant_arg = (orfquant_files && orfquant_files instanceof List && orfquant_files.size() > 0) ?
        "--orfquant ${orfquant_files.collect{ "\"${it}\"" }.join(' ')}" : ''
    def has_bedgraph = (psites_bedgraph instanceof List && psites_bedgraph.size() > 0) ||
                       (psites_bedgraph && !(psites_bedgraph instanceof List) && psites_bedgraph.name != 'NO_FILE')
    def bedgraph_arg = has_bedgraph ? "--bedgraph-dir bedgraph" : ''
    def sample_arg = (sample_list && sample_list instanceof List && sample_list.size() > 0) ? "--sample-list ${sample_list.join(',')}" : ''
    def input_files_list = (all_inputs instanceof List) ? all_inputs.collect{ it.name }.join(' ') : (all_inputs ? all_inputs.name : '')
    """
    set -uo pipefail

    if python3 -c "import Bio; import pyfaidx" 2>/dev/null; then
        echo "Dependencies already available in container"
    else
        export PYTHONUSERBASE="\$PWD/.pylibs"
        export PATH="\$PYTHONUSERBASE/bin:\$PATH"
        export PYTHONPATH="\$PYTHONUSERBASE/lib/python3.9/site-packages:\${PYTHONPATH:-}"
        export PIP_NO_CACHE_DIR=1
        mkdir -p "\$PYTHONUSERBASE"
        pip install --user --no-cache-dir --no-warn-script-location pyfaidx biopython 2>&1 || \
            python3 -m pip install --user --no-cache-dir pyfaidx biopython
    fi
    python3 -c "import Bio; import pyfaidx; print('Dependencies OK')"

    # Optionally install orfont extras (pandas + duckdb) for optimized path
    if [ -f "${run_orf_script}" ]; then
        if ! python3 -c "import pandas; import duckdb" 2>/dev/null; then
            echo "Installing orfont extras (pandas + duckdb)..."
            pip install --user --no-cache-dir --no-warn-script-location pandas duckdb 2>&1 || \
                python3 -m pip install --user --no-cache-dir pandas duckdb 2>&1 || \
                echo "WARNING: Could not install pandas/duckdb - falling back to original path"
        fi
    fi

    # Decompress any .gz input files for validation
    for f in *.gtf.gz *.bed.gz; do
        [ -f "\$f" ] && gunzip -f "\$f" || true
    done

    has_valid_input=false
    for input_file in ${input_files_list}; do
        check_file="\${input_file}"
        [ -f "\${check_file}" ] || check_file="\${input_file%.gz}"
        if [ -f "\${check_file}" ]; then
            line_count=\$(wc -l < "\${check_file}" 2>/dev/null || echo "0")
            if [ "\${line_count}" -gt 2 ] && \
               ! grep -qi "placeholder" "\${check_file}" 2>/dev/null && \
               ! grep -qi "# Empty" "\${check_file}" 2>/dev/null && \
               ! grep -qi "insufficient" "\${check_file}" 2>/dev/null; then
                has_valid_input=true
                break
            fi
        fi
    done

    if [ "\${has_valid_input}" = "false" ]; then
        echo "WARNING: All input files are empty/placeholder - creating placeholder outputs"
        for f in ${prefix}_ribotish ${prefix}_ribotricer ${prefix}_ribocode ${prefix}_orfquant; do
            echo "# Placeholder (no valid input)" > "\${f}.metadata.tsv"
            echo "# Placeholder (no valid input)" > "\${f}.bed"
            echo "# Placeholder (no valid input)" > "\${f}.gtf"
        done
        echo "# Placeholder (no valid input)" > ${prefix}.metadata.tsv
        echo "# Placeholder (no valid input)" > ${prefix}.bed
        echo "# Placeholder (no valid input)" > ${prefix}.gtf
        echo "WARNING: No valid ORF input data" > ${prefix}.stats.txt
    else
        set +e

        # Choose optimized (orfont) or original path
        if [ -f "${run_orf_script}" ] && python3 -c "import pandas; import duckdb" 2>/dev/null; then
            echo "Using optimized orfont path (run_orf.py)"
            python3 ${run_orf_script} unify \\
                --gtf ${gtf} \\
                --fasta ${fasta} \\
                --output ${prefix} \\
                --per-tool-output ${prefix} \\
                --min_len ${min_len} \\
                --threads ${task.cpus} \\
                --no-frame-merge \\
                ${ribotish_arg} \\
                ${ribotricer_arg} \\
                ${ribocode_arg} \\
                ${orfquant_arg} \\
                ${bedgraph_arg} \\
                ${sample_arg} \\
                ${extra_args} 2>&1 | tee unify_orf.log
        else
            echo "Using original unify path"
            python3 ${unify_script} \\
                --gtf ${gtf} \\
                --fasta ${fasta} \\
                --output ${prefix} \\
                --per-tool-output ${prefix} \\
                --min_len ${min_len} \\
                --threads ${task.cpus} \\
                --no-frame-merge \\
                ${ribotish_arg} \\
                ${ribotricer_arg} \\
                ${ribocode_arg} \\
                ${orfquant_arg} \\
                ${bedgraph_arg} \\
                ${sample_arg} \\
                ${extra_args} 2>&1 | tee unify_orf.log
        fi
        EXIT_CODE=\${PIPESTATUS[0]}
        set -e

        {
            echo "=== Per-Tool ORF Unification Statistics ==="
            grep -E "^===|^By Tool:|^By Sample:|^Final|^After|^Total|^Note:|^  (ribotish|ribotricer|ribocode|orfquant|Ribo-TISH|Ribotricer|RiboCode|ORFquant|[A-Za-z0-9_])" unify_orf.log || true
        } > ${prefix}.stats.txt

        if [ \${EXIT_CODE} -ne 0 ]; then
            if grep -qiE "(no valid|no ORFs|zero ORFs|empty|no predictions|0 ORFs|ValueError|KeyError|IndexError)" unify_orf.log; then
                echo "WARNING: ORF unification failed due to insufficient data - creating placeholder files"
                for f in ${prefix}_ribotish ${prefix}_ribotricer ${prefix}_ribocode ${prefix}_orfquant; do
                    [ ! -f "\${f}.metadata.tsv" ] && echo "# Placeholder" > "\${f}.metadata.tsv"
                    [ ! -f "\${f}.bed" ]          && echo "# Placeholder" > "\${f}.bed"
                    [ ! -f "\${f}.gtf" ]          && echo "# Placeholder" > "\${f}.gtf"
                done
                [ ! -f "${prefix}.metadata.tsv" ] && echo "# Placeholder" > ${prefix}.metadata.tsv
                [ ! -f "${prefix}.bed" ]           && echo "# Placeholder" > ${prefix}.bed
                [ ! -f "${prefix}.gtf" ]           && echo "# Placeholder" > ${prefix}.gtf
                echo "ERROR: unification failed" >> ${prefix}.stats.txt
            else
                echo "ERROR: UNIFY_ORF_PREDICTIONS_PER_TOOL failed with unexpected error"
                cat unify_orf.log
                exit \${EXIT_CODE}
            fi
        fi
    fi

    # Compress text outputs to save disk space
    for f in ${prefix}_ribotish.bed ${prefix}_ribotish.gtf \
             ${prefix}_ribotricer.bed ${prefix}_ribotricer.gtf \
             ${prefix}_ribocode.bed ${prefix}_ribocode.gtf \
             ${prefix}_orfquant.bed ${prefix}_orfquant.gtf \
             ${prefix}.bed ${prefix}.gtf; do
        [ -f "\$f" ] && gzip -f "\$f" || true
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
        biopython: \$(python3 -c "import Bio; print(Bio.__version__)" 2>/dev/null || echo "unknown")
        pyfaidx: \$(python3 -c "import pyfaidx; print(pyfaidx.__version__)" 2>/dev/null || echo "unknown")
        pandas: \$(python3 -c "import pandas; print(pandas.__version__)" 2>/dev/null || echo "N/A")
        duckdb: \$(python3 -c "import duckdb; print(duckdb.__version__)" 2>/dev/null || echo "N/A")
    END_VERSIONS
    """

    stub:
    def prefix = (params.unify_orf_predictions_prefix ?: 'unified_orfs').tokenize('/').last()
    """
    touch ${prefix}_ribotish.metadata.tsv ${prefix}_ribotish.bed ${prefix}_ribotish.gtf
    touch ${prefix}_ribotricer.metadata.tsv ${prefix}_ribotricer.bed ${prefix}_ribotricer.gtf
    touch ${prefix}_ribocode.metadata.tsv ${prefix}_ribocode.bed ${prefix}_ribocode.gtf
    touch ${prefix}_orfquant.metadata.tsv ${prefix}_orfquant.bed ${prefix}_orfquant.gtf
    touch ${prefix}.metadata.tsv ${prefix}.bed ${prefix}.gtf ${prefix}.stats.txt
    gzip -f ${prefix}_ribotish.bed ${prefix}_ribotish.gtf \
            ${prefix}_ribotricer.bed ${prefix}_ribotricer.gtf \
            ${prefix}_ribocode.bed ${prefix}_ribocode.gtf \
            ${prefix}_orfquant.bed ${prefix}_orfquant.gtf \
            ${prefix}.bed ${prefix}.gtf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "3.9.0"
        biopython: "1.81"
        pyfaidx: "0.7"
    END_VERSIONS
    """
}
