process QUANTIFY_ORFS {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/subread:latest' :
        'community.wave.seqera.io/library/subread:latest' }"

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
    """
    featureCounts \\
        -a ${annotation_bed} \\
        -o ${prefix}_counts.tsv \\
        -F BED \\
        -t exon \\
        -g gene_id \\
        -s 0 \\
        -T ${task.cpus} \\
        --minOverlap 1 \\
        -R BAM \\
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
