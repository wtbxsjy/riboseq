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
    def local_pkg_path = "${orfquant_pkg}"
    """
    # Ensure fasta file is available with the expected name (if it was gzipped)
    if [[ "${fasta}" == *.gz ]]; then
        gunzip -c ${fasta} > \$(basename ${fasta} .gz)
    fi

    # Install ORFquant into a task-local library if needed (works for conda runs)
    export R_LIBS_USER="${task.workDir}/Rlibs"
    mkdir -p "\$R_LIBS_USER"

    # Write R script - ORFquant should be pre-installed in custom container
    cat > run_orfquant.R <<RSCRIPTEOF
install_orfquant <- function(local_pkg_tgz = NULL, tag = "1.02") {
    work <- file.path(getwd(), "orfquant_src")
    dir.create(work, showWarnings = FALSE, recursive = TRUE)

    tgz <- local_pkg_tgz
    if (!is.null(tgz) && nzchar(tgz) && file.exists(tgz) && file.info(tgz)$size > 0) {
        message("Installing ORFquant from local tar.gz: ", tgz)
    } else {
        url <- sprintf("https://github.com/lcalviell/ORFquant/archive/refs/tags/%s.tar.gz", tag)
        tgz <- file.path(work, sprintf("ORFquant-%s.tar.gz", tag))
        message("Downloading ORFquant from GitHub: ", url)
        utils::download.file(url, tgz, mode = "wb", quiet = FALSE)
    }

    utils::untar(tgz, exdir = work, tar = "internal")
    pkg_dir <- list.dirs(work, recursive = FALSE, full.names = TRUE)
    if (length(pkg_dir) != 1) {
        stop("Unexpected ORFquant source layout in: ", work)
    }

    cmd <- sprintf("R CMD INSTALL %s", shQuote(pkg_dir[[1]]))
    message(cmd)
    status <- system(cmd)
    if (status != 0) stop("R CMD INSTALL failed with status ", status)
}

# Ensure ORFquant is available
if (!requireNamespace("ORFquant", quietly = TRUE)) {
    local_pkg <- if (${use_local_pkg ? 'TRUE' : 'FALSE'}) "${local_pkg_path}" else NULL
    tryCatch({
        install_orfquant(local_pkg_tgz = local_pkg, tag = "1.02")
    }, error = function(e) {
        stop(
            "ORFquant is not installed and automatic installation failed: ", e$message, "\n",
            "Provide a pre-downloaded tarball with --orfquant_pkg (e.g. ORFquant_1.02.0.tar.gz from lcalviell/ORFquant), ",
            "or use a container with ORFquant pre-installed (e.g. --orfquant_container)."
        )
    })

    if (!requireNamespace("ORFquant", quietly = TRUE)) {
        stop("ORFquant install completed but package is still not available on library paths.")
    }
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
