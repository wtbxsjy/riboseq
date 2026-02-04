process RIBOTRICER_DETECTORFS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ribotricer:1.3.3--pyhdfd78af_0':
        'biocontainers/ribotricer:1.3.3--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(candidate_orfs)

    output:
    tuple val(meta), path('*_protocol.txt')             , emit: protocol, optional: true
    tuple val(meta), path('*_bam_summary.txt')          , emit: bam_summary
    tuple val(meta), path('*_read_length_dist.pdf')     , emit: read_length_dist
    tuple val(meta), path('*_metagene_profiles_5p.tsv') , emit: metagene_profile_5p
    tuple val(meta), path('*_metagene_profiles_3p.tsv') , emit: metagene_profile_3p
    tuple val(meta), path('*_metagene_plots.pdf')       , emit: metagene_plots
    tuple val(meta), path('*_psite_offsets.txt')        , emit: psite_offsets, optional: true
    tuple val(meta), path('*_pos.wig')                  , emit: pos_wig
    tuple val(meta), path('*_neg.wig')                  , emit: neg_wig
    tuple val(meta), path('*_translating_ORFs.tsv')     , emit: orfs
    path "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def strandedness_cmd = ''

    switch(meta.strandedness) {
        case "forward":
            strandedness_cmd = "--stranded yes"
            break
        case "reverse":
            strandedness_cmd = "--stranded reverse"
            break
        //
        // Specifying unstranded seems broken - see
        // https://github.com/smithlabcode/ribotricer/issues/153. Leaving it
        // undefined works, though ribotricer may incorrectly infer
        // strandednesss?
        //
        //case "unstranded":
        //    strandedness_cmd = "--stranded no"
        //    break
    }
    """
    # Wrap ribotricer in error handling for low-quality samples
    set +e
    ribotricer detect-orfs \\
        --bam $bam \\
        --ribotricer_index $candidate_orfs \\
        --prefix $prefix \\
        $strandedness_cmd \\
        $args 2>&1 | tee ribotricer.log
    
    EXIT_CODE=\${PIPESTATUS[0]}
    set -e
    
    # Check for low signal indicators
    if [ \$EXIT_CODE -ne 0 ]; then
        if grep -q "no periodic read length found" ribotricer.log || \
           grep -q "WARNING.*periodic" ribotricer.log; then
            echo "WARNING: Ribotricer failed due to insufficient periodic signal - creating placeholder files"
            
            # Create placeholder files with informative headers
            echo "# No ribotricer results - insufficient periodic signal detected" > ${prefix}_protocol.txt
            echo "# No ribotricer results - insufficient periodic signal detected" > ${prefix}_bam_summary.txt
            echo "# No ribotricer results - insufficient periodic signal detected" > ${prefix}_metagene_profiles_5p.tsv
            echo "# No ribotricer results - insufficient periodic signal detected" > ${prefix}_metagene_profiles_3p.tsv
            echo "# No ribotricer results - insufficient periodic signal detected" > ${prefix}_psite_offsets.txt
            echo "# No ribotricer results - insufficient periodic signal detected" > ${prefix}_pos.wig
            echo "# No ribotricer results - insufficient periodic signal detected" > ${prefix}_neg.wig
            echo "# No ribotricer results - insufficient periodic signal detected" > ${prefix}_translating_ORFs.tsv
            
            # Create empty PDFs (using touch as placeholder)
            touch ${prefix}_read_length_dist.pdf
            touch ${prefix}_metagene_plots.pdf
        else
            echo "ERROR: Ribotricer failed with unexpected error"
            exit \$EXIT_CODE
        fi
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ribotricer: \$(ribotricer --version 2>&1 | grep ribotricer | sed '1!d ; s/ribotricer, version //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_protocol.txt
    touch ${prefix}_bam_summary.txt
    touch ${prefix}_read_length_dist.pdf
    touch ${prefix}_metagene_profiles_5p.tsv
    touch ${prefix}_metagene_profiles_3p.tsv
    touch ${prefix}_metagene_plots.pdf
    touch ${prefix}_psite_offsets.txt
    touch ${prefix}_pos.wig
    touch ${prefix}_neg.wig
    touch ${prefix}_translating_ORFs.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ribotricer: \$(ribotricer --version 2>&1 | grep ribotricer | sed '1!d ; s/ribotricer, version //')
    END_VERSIONS
    """
}
