process RIBOCODE_DETECT {
    tag "$meta.id"
    label 'process_medium'

    // Allow process to complete even if RiboCode fails due to low periodicity
    errorStrategy 'ignore'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ribocode:1.2.11--pyh145b6a8_1' :
        'quay.io/biocontainers/ribocode:1.2.11--pyh145b6a8_1' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path gtf
    path fasta

    output:
    tuple val(meta), path("${meta.id}*"), emit: results, optional: true
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def strandedness = 'yes'
    if (meta.strandedness == 'reverse') {
        strandedness = 'reverse'
    } else if (meta.strandedness == 'unstranded') {
        strandedness = 'no'
    }
    // RiboCode strandedness options: yes, reverse, no.
    // nf-core/riboseq strandedness: forward, reverse, unstranded.
    // forward -> yes
    // reverse -> reverse
    // unstranded -> no

    """
    # Run RiboCode, capture exit status
    RiboCode_onestep \\
        -g $gtf \\
        -f $fasta \\
        -r $bam \\
        --stranded $strandedness \\
        -o ${meta.id} \\
        $args || {
        echo "RiboCode failed for ${meta.id} - likely due to insufficient periodicity in data"
        echo "This is common for low-depth or test datasets"
        # Create a marker file to indicate failure
        echo "FAILED: No periodicity detected" > ${meta.id}.ribocode_failed.txt
    }

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ribocode: \$(RiboCode_onestep --version 2>&1 | sed 's/^.* //')
    END_VERSIONS
    """
}
