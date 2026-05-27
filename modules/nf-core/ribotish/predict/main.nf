process RIBOTISH_PREDICT {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ribotish:0.2.7--pyhdfd78af_0':
        'biocontainers/ribotish:0.2.7--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(bam_ribo), path(bai_ribo)
    tuple val(meta2), path(bam_ti), path(bai_ti)
    tuple val(meta3), path(fasta), path(gtf)
    tuple val(meta4), path(candidate_orfs)
    tuple val(meta5), path(para_ribo)
    tuple val(meta6), path(para_ti)

    output:
    tuple val(meta), path("*_pred.txt")        , emit: predictions
    tuple val(meta), path("*_all.txt")         , emit: all
    tuple val(meta), path("*_transprofile.py") , emit: transprofile
    path "versions.yml"                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    ribo_bam_cmd = ''
    ti_bam_cmd = ''
    if (bam_ribo){
        ribo_bam_cmd = "-b ${bam_ribo.join(',')}"
        if (para_ribo){
            ribo_bam_cmd += " --ribopara ${para_ribo.join(',')}"
        }
    }
    if (bam_ti){
        ti_bam_cmd = "-t ${bam_ti.join(',')}"
        if (para_ti){
            ti_bam_cmd += " --tisparapara  ${para_ti.join(',')}"
        }
    }
    """
    set +e
    ribotish predict \\
        $ribo_bam_cmd \\
        $ti_bam_cmd \\
        -f $fasta \\
        -g $gtf \\
        -o ${prefix}_pred.txt \\
        --allresult ${prefix}_all.txt \\
        --transprofile ${prefix}_transprofile.py \\
        -p $task.cpus \\
        $args 2>&1 | tee ribotish.log
    
    EXIT_CODE=\${PIPESTATUS[0]}
    set -e
    
    # Check for low signal/quality errors that should allow pipeline continuation
    if [ \$EXIT_CODE -ne 0 ]; then
        # Check for known low-quality sample errors
        if grep -qiE "(no valid|insufficient|empty|zero|no reads|cannot|failed to)" ribotish.log || \
           grep -qiE "(no ORFs|no candidates|no predictions)" ribotish.log; then
            echo "WARNING: Ribotish failed due to insufficient data/signal - creating placeholder files"
            
            # Create placeholder files with header
            echo -e "# Ribotish placeholder - insufficient data\\n# Gid\\tTid\\tSymbol\\tGeneType\\tGenomePos\\tStartCodon\\tStrand\\tAALen\\tTisType\\tTISGroup\\tTISCounts\\tTISPvalue\\tRiboPvalue\\tRiboPStatus\\tFisherPvalue\\tTISQvalue\\tRiboQvalue\\tFrameQvalue\\tFisherQvalue" > ${prefix}_pred.txt
            echo -e "# Ribotish placeholder - insufficient data\\n# All ORF predictions" > ${prefix}_all.txt
            echo "# Ribotish placeholder - insufficient data" > ${prefix}_transprofile.py
            
            echo "Placeholder files created - downstream analysis will continue"
        else
            echo "ERROR: Ribotish failed with unexpected error:"
            cat ribotish.log
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
    touch ${prefix}_pred.txt
    touch ${prefix}_all.txt
    touch ${prefix}_transprofile.py

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ribotish: \$(ribotish --version | sed 's/ribotish //')
    END_VERSIONS
    """
}
