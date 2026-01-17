process GENCODE_ORF_MAPPER {
    tag "${project_id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/YOUR_TOOL:YOUR_VERSION' :
        'nfcore/gencode-orf-mapper:1.1.0' }"

    input:
    path orfs_fasta
    path orfs_bed
    path ensembl_dir
    val project_id

    output:
    path "*.orfs.fa"       , emit: unified_fasta
    path "*.orfs.bed"      , emit: unified_bed
    path "*.orfs.gtf"      , emit: unified_gtf
    path "*.orfs.out"      , emit: unified_table
    path "*.altmapped"     , optional: true, emit: altmapped
    path "*.unmapped"      , optional: true, emit: unmapped
    path "versions.yml"    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${project_id}"

    // Parameters with defaults
    def min_len = params.gencode_orf_min_length ?: 16
    def max_len = params.gencode_orf_max_length ?: 999999
    def collapse_thr = params.gencode_collapse_threshold ?: 0.9
    def collapse_method = params.gencode_collapse_method ?: 'longest_string'
    def add_cds = params.gencode_add_cds ? 'yes' : 'no'

    """
    # Run GENCODE ORF mapper
    python3 /opt/gencode-riboseqORFs/ORF_mapper_to_GENCODE_v1.1.py \\
        -d ${ensembl_dir} \\
        -f ${orfs_fasta} \\
        -b ${orfs_bed} \\
        -o ${prefix} \\
        -l ${min_len} \\
        -L ${max_len} \\
        -c ${collapse_thr} \\
        -m ${collapse_method} \\
        -C ${add_cds} \\
        ${args}

    # Generate versions file
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gencode_orf_mapper: "1.1.0"
        python: \$(python3 --version | sed 's/Python //')
        biopython: \$(python3 -c "import Bio; print(Bio.__version__)")
        bedtools: \$(bedtools --version | sed 's/bedtools v//')
        gffread: \$(gffread --version 2>&1 | head -n1 | sed 's/gffread //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${project_id}"
    """
    touch ${prefix}.orfs.fa
    touch ${prefix}.orfs.bed
    touch ${prefix}.orfs.gtf
    touch ${prefix}.orfs.out

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gencode_orf_mapper: "1.1.0"
        python: "3.9.0"
        biopython: "1.81"
        bedtools: "2.30.0"
        gffread: "0.12.7"
    END_VERSIONS
    """
}
