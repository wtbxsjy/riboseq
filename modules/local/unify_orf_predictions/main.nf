process UNIFY_ORF_PREDICTIONS {
    def prefix = (params.unify_orf_predictions_prefix ?: 'unified_orfs').tokenize('/').last()

    tag "${prefix}"
    label 'process_medium'

    publishDir "${params.outdir}/orf_unification", mode: params.publish_dir_mode

    conda "${moduleDir}/environment.yml"
    // Use biopython container which includes pip and allows installing pyfaidx
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.79' :
        'quay.io/biocontainers/biopython:1.79' }"

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
    // Files are already staged via all_inputs, just need filenames properly quoted
    def ribotish_arg = (ribotish_files && ribotish_files instanceof List && ribotish_files.size() > 0) ? 
        "--ribotish ${ribotish_files.collect{ "\"${it}\"" }.join(' ')}" : ''
    def ribotricer_arg = (ribotricer_files && ribotricer_files instanceof List && ribotricer_files.size() > 0) ? 
        "--ribotricer ${ribotricer_files.collect{ "\"${it}\"" }.join(' ')}" : ''
    def orfquant_arg = (orfquant_files && orfquant_files instanceof List && orfquant_files.size() > 0) ? 
        "--orfquant ${orfquant_files.collect{ "\"${it}\"" }.join(' ')}" : ''
    """
    set -euo pipefail

    # Setup user-local Python package directory to avoid permission issues
    export PYTHONUSERBASE="\$PWD/.pylibs"
    export PATH="\$PYTHONUSERBASE/bin:\$PATH"
    export PYTHONPATH="\$PYTHONUSERBASE/lib/python3.9/site-packages:\${PYTHONPATH:-}"
    # Disable pip cache entirely to avoid permission issues in containers
    export PIP_NO_CACHE_DIR=1
    mkdir -p "\$PYTHONUSERBASE"

    # Install pyfaidx (biopython should be in container)
    echo "Installing Python dependencies..."
    pip install --user --no-cache-dir --no-warn-script-location pyfaidx 2>&1 || {
        echo "pip install failed, trying with python -m pip..."
        python3 -m pip install --user --no-cache-dir pyfaidx 2>&1 || true
    }
    
    # If biopython is not in container, install it too
    python3 -c "import Bio" 2>/dev/null || {
        echo "Installing biopython..."
        pip install --user --no-cache-dir --no-warn-script-location biopython 2>&1 || python3 -m pip install --user --no-cache-dir biopython
    }
    
    # Verify installation
    echo "Verifying dependencies..."
    python3 -c "import Bio; import pyfaidx; print('Dependencies OK: Bio=' + Bio.__version__ + ', pyfaidx=' + pyfaidx.__version__)"

    python3 ${unify_script} \\
        --gtf ${gtf} \\
        --fasta ${fasta} \\
        --output ${prefix} \\
        --min_len ${min_len} \\
        ${ribotish_arg} \\
        ${ribotricer_arg} \\
        ${orfquant_arg} \\
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
