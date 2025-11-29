process ORFQUANT_RUN {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    // Use RiboseQC container since ORFquant needs RiboseQC to load annotation files
    // ORFquant will be installed at runtime from GitHub
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/riboseqc:1.1--r36_1' :
        'quay.io/biocontainers/riboseqc:1.1--r36_1' }"

    input:
    tuple val(meta), path(for_orfquant)   // *_for_ORFquant file from RiboseQC
    path annotation                        // *_Rannot file from RiboseQC/ORFquant annotation
    path fasta                             // Genome fasta file

    output:
    tuple val(meta), path("*_final_ORFquant_results")  , emit: results
    tuple val(meta), path("*_Detected_ORFs.gtf")       , emit: gtf, optional: true
    tuple val(meta), path("*_Protein_sequences.fasta") , emit: proteins, optional: true
    tuple val(meta), path("*_tmp_ORFquant_results")    , emit: tmp_results, optional: true
    tuple val(meta), path("*_ORFquant_plots_RData")    , emit: plots_data, optional: true
    tuple val(meta), path("*_plots")                   , emit: plots_dir, optional: true
    path "versions.yml"                                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def n_cores = task.cpus ?: 1
    // Parse optional arguments
    def write_gtf = args.contains('write_GTF_file=FALSE') ? 'FALSE' : 'TRUE'
    def write_fasta = args.contains('write_protein_fasta=FALSE') ? 'FALSE' : 'TRUE'
    def write_tmp = args.contains('write_temp_files=FALSE') ? 'FALSE' : 'TRUE'
    def plot_results = args.contains('plot_results=TRUE') ? 'TRUE' : 'FALSE'
    """
    # Ensure fasta file is available with the expected name (if it was gzipped)
    # The annotation might refer to the uncompressed name
    if [[ "${fasta}" == *.gz ]]; then
        gunzip -c ${fasta} > \$(basename ${fasta} .gz)
    fi

    # Write R script to file
    cat <<'EOF' > run_orfquant.R
    # Install ORFquant from GitHub if not available (needed for container mode)
    if (!requireNamespace("ORFquant", quietly = TRUE)) {
        message("Installing ORFquant from GitHub...")
        if (!requireNamespace("remotes", quietly = TRUE)) {
            install.packages("remotes", repos = "https://cloud.r-project.org", quiet = TRUE)
        }
        remotes::install_github("ohlerlab/ORFquant@v1.1", quiet = FALSE, upgrade = "never")
    }

    library(ORFquant)

    # Run ORFquant
    run_ORFquant(
        for_ORFquant_file = "${for_orfquant}",
        annotation_file = "${annotation}",
        n_cores = ${n_cores},
        prefix = "${prefix}",
        write_temp_files = ${write_tmp},
        write_GTF_file = ${write_gtf},
        write_protein_fasta = ${write_fasta},
        interactive = FALSE
    )

    # Optionally generate plots
    if (${plot_results}) {
        tryCatch({
            plot_ORFquant_results(
                for_ORFquant_file = "${for_orfquant}",
                ORFquant_output_file = paste0("${prefix}", "_final_ORFquant_results"),
                annotation_file = "${annotation}",
                output_plots_path = paste0("${prefix}", "_plots"),
                prefix = "${prefix}"
            )
        }, error = function(e) {
            message("Warning: Could not generate ORFquant plots: ", e\$message)
        })
    }

    # Write versions
    writeLines(
        c(
            '"${task.process}":',
            paste0('    orfquant: "', packageVersion("ORFquant"), '"'),
            paste0('    r-base: "', R.Version()\$major, ".", R.Version()\$minor, '"')
        ),
        "versions.yml"
    )
    EOF

    # Run using the explicit Rscript path from Conda
    \$CONDA_PREFIX/bin/Rscript run_orfquant.R
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_final_ORFquant_results
    touch ${prefix}_Detected_ORFs.gtf
    touch ${prefix}_Protein_sequences.fasta
    touch ${prefix}_tmp_ORFquant_results
    mkdir -p ${prefix}_plots
    touch ${prefix}_ORFquant_plots_RData

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        orfquant: "1.02"
        r-base: "4.3"
    END_VERSIONS
    """
}
