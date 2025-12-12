process ORFQUANT_RUN {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    // Use custom container with ORFquant pre-installed
    // Build from: containers/Singularity.orfquant.def
    // Or specify via params.orfquant_container
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        (params.orfquant_container ?: 'https://depot.galaxyproject.org/singularity/riboseqc:1.1--r36_1') :
        'quay.io/biocontainers/riboseqc:1.1--r36_1' }"

    input:
    tuple val(meta), path(for_orfquant)   // *_for_ORFquant file from RiboseQC
    path annotation                        // *_Rannot file from RiboseQC/ORFquant annotation
    path fasta                             // Genome fasta file
    path orfquant_pkg                      // Pre-downloaded ORFquant R package (tar.gz) - optional

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
    def use_local_pkg = orfquant_pkg.name != 'NO_FILE'
    def local_pkg_path = orfquant_pkg.name
    """
    # Ensure fasta file is available with the expected name (if it was gzipped)
    if [[ "${fasta}" == *.gz ]]; then
        gunzip -c ${fasta} > \$(basename ${fasta} .gz)
    fi

    # Write R script - ORFquant should be pre-installed in custom container
    cat > run_orfquant.R <<RSCRIPTEOF
# Check if ORFquant is available
if (!requireNamespace("ORFquant", quietly = TRUE)) {
    stop("ORFquant is not installed. Please use a container with ORFquant pre-installed or provide --orfquant_pkg parameter.")
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
        message("Warning: Could not generate ORFquant plots: ", e\\\$message)
    })
}

# Write versions
writeLines(
    c(
        '"${task.process}":',
        paste0('    orfquant: "', packageVersion("ORFquant"), '"'),
        paste0('    r-base: "', R.Version()\\\$major, ".", R.Version()\\\$minor, '"')
    ),
    "versions.yml"
)
RSCRIPTEOF

    # Run using Rscript
    Rscript run_orfquant.R
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
