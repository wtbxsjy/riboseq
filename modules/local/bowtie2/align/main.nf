process BOWTIE2_ALIGN {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::bowtie2=2.4.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bowtie2:2.4.4--py39hbb4e92a_0' :
        'quay.io/biocontainers/bowtie2:2.4.4--py39hbb4e92a_0' }"

    input:
    tuple val(meta), path(reads)
    path index

    output:
    tuple val(meta), path("*.fastq.gz"), emit: clean_reads
    tuple val(meta), path("*.log")     , emit: log
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def input_reads = ''
    if (meta.single_end) {
        def read_list = reads instanceof Path ? [ reads ] : reads
        input_reads = read_list.collect { fastq -> "-U \"${fastq}\"" }.join(' ')
    } else {
        input_reads = "-1 \"${reads[0]}\" -2 \"${reads[1]}\""
    }

    def unaligned = meta.single_end ? "--un-gz ${prefix}.clean.fastq.gz" : "--un-conc-gz ${prefix}.clean.%.fastq.gz"

    """
    INDEX=`find -L ${index} -name "*.1.bt2" | sed 's/.1.bt2//' | head -n 1`
    [ -z "\$INDEX" ] && INDEX=`find -L ${index} -name "*.1.bt2l" | sed 's/.1.bt2l//' | head -n 1`

    bowtie2 \\
        $args \\
        --threads $task.cpus \\
        -x \$INDEX \\
        $input_reads \\
        $unaligned \\
        -S /dev/null \\
        2> ${prefix}.bowtie2.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie2: \$(echo \$(bowtie2 --version 2>&1) | sed 's/^.*bowtie2-align-s version //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (meta.single_end) {
        """
        touch ${prefix}.clean.fastq.gz
        touch ${prefix}.bowtie2.log
        touch versions.yml
        """
    } else {
        """
        touch ${prefix}.clean_1.fastq.gz
        touch ${prefix}.clean_2.fastq.gz
        touch ${prefix}.bowtie2.log
        touch versions.yml
        """
    }
}
