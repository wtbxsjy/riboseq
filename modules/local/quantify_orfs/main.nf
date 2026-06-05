process QUANTIFY_ORFS {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/subread:2.1.1--h577a1d6_0' :
        'quay.io/biocontainers/subread:2.1.1--h577a1d6_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path annotation_bed  // ORF annotation in BED12 format for featureCounts

    output:
    tuple val(meta), path("*_counts.tsv"), emit: counts
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def paired = meta.single_end ? '' : '-p --countReadPairs'
    """
    # Convert BED to SAF format (subread 2.1.1 lacks BED support)
    awk 'BEGIN{OFS="\\t"; print "GeneID\\tChr\\tStart\\tEnd\\tStrand"}
         {print \$4, \$1, \$2+1, \$3, \$6}' ${annotation_bed} > annotation.saf

    featureCounts \\
        -a annotation.saf \\
        -o ${prefix}_counts.tsv \\
        -F SAF \\
        -t exon \\
        -g GeneID \\
        -s 0 \\
        -T ${task.cpus} \\
        --minOverlap 1 \\
        ${paired} \\
        ${bam}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        subread: \$(featureCounts -v 2>&1 | head -1 | sed 's/featureCounts //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_counts.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        subread: "2.0.6"
    END_VERSIONS
    """
}
