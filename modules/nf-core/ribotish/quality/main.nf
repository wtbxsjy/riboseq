process RIBOTISH_QUALITY {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ribotish:0.2.7--pyhdfd78af_0':
        'biocontainers/ribotish:0.2.7--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(gtf)

    output:
    tuple val(meta), path("*.txt")    , emit: distribution
    tuple val(meta), path("*.pdf")    , emit: pdf
    tuple val(meta), path("*.para.py"), emit: offset
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    set +e
    ribotish quality \\
        -b $bam \\
        -g $gtf \\
        -o ${prefix}_qual.txt \\
        -f ${prefix}_qual.pdf \\
        -r ${prefix}.para.py \\
        -p $task.cpus \\
        $args 2>&1 | tee ribotish_quality.log
    
    EXIT_CODE=\${PIPESTATUS[0]}
    set -e
    
    # Check for low signal/quality errors that should allow pipeline continuation
    if [ \$EXIT_CODE -ne 0 ]; then
        # Check for known low-quality sample errors (no reads, insufficient signal)
        if grep -qiE "(no reads found|Counted reads: 0|no valid|insufficient|empty)" ribotish_quality.log; then
            echo "WARNING: Ribotish quality failed due to insufficient reads/signal - creating placeholder files"
            
            # Create placeholder quality distribution file
            echo "# Ribotish quality placeholder - insufficient reads/signal for sample ${prefix}" > ${prefix}_qual.txt
            echo "# No read length distribution available" >> ${prefix}_qual.txt
            echo "read_length\tcount\tproportion" >> ${prefix}_qual.txt
            
            # Create placeholder PDF (empty file, as we can't generate a real PDF)
            echo "# Placeholder - no quality plot generated due to insufficient reads" > ${prefix}_qual.pdf
            
            # Create placeholder offset parameter file with default values
            cat > ${prefix}.para.py << 'PYEOF'
# Ribotish quality placeholder - insufficient reads/signal
# Using default ribosome footprint parameters
# These are typical values for ribosome profiling data
offdict = {28: 12, 29: 12, 30: 12, 31: 12, 32: 12}
PYEOF
            
            echo "Placeholder files created - downstream analysis will use default parameters"
        else
            echo "ERROR: Ribotish quality failed with unexpected error:"
            cat ribotish_quality.log
            exit \$EXIT_CODE
        fi
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ribotish: \$(ribotish --version | sed 's/ribotish //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_qual.txt ${prefix}_qual.pdf ${prefix}.para.py

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ribotish: \$(ribotish --version | sed 's/ribotish //')
    END_VERSIONS
    """
}
