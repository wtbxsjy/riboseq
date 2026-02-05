process UNIFY_ORF_PREDICTIONS {
    def prefix = (params.unify_orf_predictions_prefix ?: 'unified_orfs').tokenize('/').last()

    tag "${prefix}"
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
    tuple val(ribotish_files), val(ribotricer_files), val(orfquant_files), path(all_inputs)
    path gtf
    path fasta
    path unify_script
    path psites_bedgraph, stageAs: 'bedgraph/*'  // RiboseQC P-site bedgraph files (optional)
    val sample_list       // List of sample names for bedgraph files

    output:
    path "${prefix}.metadata.tsv", emit: metadata
    path "${prefix}.bed"         , emit: bed
    path "${prefix}.gtf"         , emit: gtf
    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def min_len = params.unify_orf_min_len ?: 24
    def extra_args = params.extra_unify_orf_predictions_args ?: ''
    // Advanced merging parameters
    def merge_tolerance = params.unify_orf_merge_tolerance ?: 3
    def min_overlap = params.unify_orf_min_overlap ?: 0.5
    def no_frame_merge = params.unify_orf_no_frame_merge ? "--no-frame-merge" : ''
    def no_overlap_group = params.unify_orf_no_overlap_group ? "--no-overlap-group" : ''
    // Files are already staged via all_inputs, just need filenames properly quoted
    def ribotish_arg = (ribotish_files && ribotish_files instanceof List && ribotish_files.size() > 0) ? 
        "--ribotish ${ribotish_files.collect{ "\"${it}\"" }.join(' ')}" : ''
    def ribotricer_arg = (ribotricer_files && ribotricer_files instanceof List && ribotricer_files.size() > 0) ? 
        "--ribotricer ${ribotricer_files.collect{ "\"${it}\"" }.join(' ')}" : ''
    def orfquant_arg = (orfquant_files && orfquant_files instanceof List && orfquant_files.size() > 0) ? 
        "--orfquant ${orfquant_files.collect{ "\"${it}\"" }.join(' ')}" : ''
    // Bedgraph arguments for P-site statistics from RiboseQC
    // psites_bedgraph can be a list or path, check if it has files
    def has_bedgraph = (psites_bedgraph instanceof List && psites_bedgraph.size() > 0) || 
                       (psites_bedgraph && !(psites_bedgraph instanceof List) && psites_bedgraph.name != 'NO_FILE')
    def bedgraph_arg = has_bedgraph ? "--bedgraph-dir bedgraph" : ''
    def sample_arg = (sample_list && sample_list instanceof List && sample_list.size() > 0) ? "--sample-list ${sample_list.join(',')}" : ''
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

    # Check if all input files are empty/placeholder files
    has_valid_input=false
    for input_file in \${all_inputs:-}; do
        if [ -f "\${input_file}" ]; then
            # Check if file has actual content (not just header or placeholder)
            line_count=\$(wc -l < "\${input_file}" || echo "0")
            if [ "\${line_count}" -gt 1 ]; then
                # Check it's not a placeholder file
                if ! grep -q "# Placeholder" "\${input_file}" 2>/dev/null && \\
                   ! grep -q "# Empty" "\${input_file}" 2>/dev/null && \\
                   ! grep -q "placeholder" "\${input_file}" 2>/dev/null; then
                    has_valid_input=true
                    break
                fi
            fi
        fi
    done

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
        
        echo "WARNING: UNIFY_ORF_PREDICTIONS created placeholder outputs due to insufficient input data"
    else
        set +e
        python3 ${unify_script} \\
            --gtf ${gtf} \\
            --fasta ${fasta} \\
            --output ${prefix} \\
            --min_len ${min_len} \\
            --threads ${task.cpus} \\
            --merge-tolerance ${merge_tolerance} \\
            --min-overlap ${min_overlap} \\
            ${no_frame_merge} \\
            ${no_overlap_group} \\
            ${ribotish_arg} \\
            ${ribotricer_arg} \\
            ${orfquant_arg} \\
            ${bedgraph_arg} \\
            ${sample_arg} \\
            ${extra_args} 2>&1 | tee unify_orf.log
        EXIT_CODE=\${PIPESTATUS[0]}
        set -e
        
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
            else
                echo "ERROR: UNIFY_ORF_PREDICTIONS failed with unexpected error"
                cat unify_orf.log
                exit \${EXIT_CODE}
            fi
        fi
    fi

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
