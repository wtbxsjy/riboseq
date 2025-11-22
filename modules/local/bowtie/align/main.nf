process BOWTIE_ALIGN {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::bowtie=1.3.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bowtie:1.3.1--py38h1b792b2_3' :
        'quay.io/biocontainers/bowtie:1.3.1--py38h1b792b2_3' }"

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

    def read_args = ''
    if (meta.single_end) {
        def read_list = reads instanceof Path ? [ reads ] : reads
        read_args = read_list.collect { fastq -> "\"${fastq}\"" }.join(' ')
    } else {
        read_args = "-1 \"${reads[0]}\" -2 \"${reads[1]}\""
    }

    def unaligned_flag = meta.single_end ? "--un ${prefix}.clean.fastq" : "--un-conc ${prefix}.clean.fastq"

    """
    INDEX=`find -L ${index} -name "*.1.ebwt" | sed 's/.1.ebwt//'`

    bowtie \\
        $args \\
        -q \\
        -S \\
        -p $task.cpus \\
        $unaligned_flag \\
        \$INDEX \\
        $read_args \\
        > /dev/null \\
        2> ${prefix}.bowtie.log

    if [ -f ${prefix}.clean.fastq ]; then
        gzip -c ${prefix}.clean.fastq > ${prefix}.clean.fastq.gz
        rm ${prefix}.clean.fastq
    fi

    if [ -f ${prefix}.clean.fastq.1 ]; then
        gzip -c ${prefix}.clean.fastq.1 > ${prefix}.clean_1.fastq.gz
        gzip -c ${prefix}.clean.fastq.2 > ${prefix}.clean_2.fastq.gz
        rm ${prefix}.clean.fastq.1 ${prefix}.clean.fastq.2
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie: \$(echo \$(bowtie --version 2>&1) | sed 's/^.*bowtie-align-s version //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (meta.single_end) {
        """
        touch ${prefix}.clean.fastq.gz
        touch ${prefix}.bowtie.log
        touch versions.yml
        """
    } else {
        """
        touch ${prefix}.clean_1.fastq.gz
        touch ${prefix}.clean_2.fastq.gz
        touch ${prefix}.bowtie.log
        touch versions.yml
        """
    }
}
