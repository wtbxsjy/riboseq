process RPBP_PREPARE_GENOME {
    tag "$fasta"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        (params.rpbp_container ?: 'https://depot.galaxyproject.org/singularity/rpbp:4.0.1--py312hf731ba3_0') :
        'biocontainers/rpbp:4.0.1--py312hf731ba3_0' }"

    input:
    path fasta
    path gtf
    path ribosomal_fasta

    output:
    path "transcript-index/genome.orfs-genomic.bed.gz", emit: orfs_genomic
    path "transcript-index/genome.orfs-exons.bed.gz"  , emit: orfs_exons
    path "versions.yml"                               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    mkdir -p ribosomal_index star_index

    # Create a config file for rp-bp
    cat <<EOF > config.yaml
    genome_base_path: .
    genome_name: genome
    gtf: $gtf
    fasta: $fasta
    ribosomal_fasta: $ribosomal_fasta
    ribosomal_index: ./ribosomal_index/rRNA
    star_index: ./star_index
    EOF

    # Run prepare-rpbp-genome
    # We use --num-cpus to speed up if it does parallel processing
    prepare-rpbp-genome \\
        config.yaml \\
        --num-cpus $task.cpus \\
        --overwrite \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rpbp: \$(python3 -c "import rpbp; print(rpbp.__version__)")
    END_VERSIONS
    """
}
