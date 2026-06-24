/*
 * Prepare for_ORFquant file with corrected P-site offsets
 * 
 * This module regenerates the for_ORFquant file using the corrected
 * read length to P-site offset mapping extracted from RiboseQC results.
 * Using the corrected offsets can improve ORFquant accuracy.
 */

process PREPARE_FOR_ORFQUANT_CORRECTED {
    tag "$meta.id"
    label 'process_medium'

    // Use same conda environment and container as ORFQUANT_RUN
    conda "../../orfquant/environment.yml"
    // Use custom container with ORFquant pre-installed (same as ORFQUANT_RUN)
    // Build from: containers/Singularity.orfquant.def
    // Or specify via params.orfquant_container
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        (params.orfquant_container ?: 'https://depot.galaxyproject.org/singularity/riboseqc:1.1--r36_1') :
        'quay.io/biocontainers/riboseqc:1.1--r36_1' }"

    input:
    tuple val(meta), path(bam), path(bai)  // BAM and index files
    path annotation                          // *_Rannot file from RIBOSEQC_PREPAREANNOTATION
    tuple val(meta2), path(rl_cutoff)        // read_length to cutoff mapping from EXTRACT_RL_CUTOFF
    path bsgenome_dir                        // BSgenome package directory from PREPAREANNOTATION (optional)

    output:
    // prepare_for_ORFquant appends "_for_ORFquant" to dest_name
    // So output is: ${prefix}_corrected_for_ORFquant
    tuple val(meta), path("*_corrected_for_ORFquant"), emit: for_orfquant
    path "versions.yml"                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def chunk_size = args.contains('chunk_size=') ? args.replaceAll(/.*chunk_size=(\d+).*/, '$1') : '5e+06'
    """
    #!/bin/bash

    # Append a task-local Rlibs directory so that packages in the container remain accessible
    _local_rlibs="${task.workDir}/Rlibs"
    mkdir -p "\$_local_rlibs"
    export R_LIBS_USER="\${_local_rlibs}\${R_LIBS_USER:+:\${R_LIBS_USER}}"

# Export BSgenome dir path to R via env var (handle optional input)
if [ -n "${bsgenome_dir:-}" ] && [ "${bsgenome_dir}" != "NO_FILE" ] && [ -d "${bsgenome_dir}" ]; then
    export BSGENOME_DIR="${bsgenome_dir}"
fi

    cat <<'RSCRIPTEOF' > script.R
# Install BSgenome from PREPAREANNOTATION output (if available)
_bsgenome_dir <- Sys.getenv("BSGENOME_DIR", unset = "")
if (nchar(_bsgenome_dir) > 0 && dir.exists(_bsgenome_dir)) {
    cat("Installing BSgenome for non-model organism...\n")
    _rlibs_local <- file.path(getwd(), "rlibs")
    dir.create(_rlibs_local, showWarnings = FALSE, recursive = TRUE)
    .libPaths(c(_rlibs_local, .libPaths()))
    Sys.setenv(R_LIBS_USER = _rlibs_local)
    install.packages(_bsgenome_dir, repos = NULL, type = "source", lib = _rlibs_local, quiet = TRUE)
    cat("BSgenome installed to", _rlibs_local, "\n")
}

# Load ORFquant library
library(ORFquant)

# Patch ORFquant::load_annotation to handle NULL genome_package
# (happens with forge_BSgenome=FALSE for non-model organisms).
# 2026-06-24: is(genome, "FaFile") returns FALSE for FaFile_Circ in container.
# Fix: check genome_package first, use genome field as fallback.
fix_load_annotation <- function() {
    ns <- asNamespace("ORFquant")
    unlockBinding("load_annotation", ns)
    patched_load_annotation <- function(path) {
        GTF_annotation <- get(load(path))
        genome_pkg <- GTF_annotation\$genome_package
        if (!is.null(genome_pkg) && nchar(genome_pkg) > 0) {
            # Traditional BSgenome package reference
            library(genome_pkg, character.only = TRUE)
            genome_sequence <- get(genome_pkg)
        } else if (is.character(GTF_annotation\$genome) && nchar(GTF_annotation\$genome) > 0) {
            # forge_BSgenome=TRUE stores package name as string in genome field
            pkg_name <- GTF_annotation\$genome
            library(pkg_name, character.only = TRUE)
            genome_sequence <- get(pkg_name)
        } else if (!is.null(GTF_annotation\$genome)) {
            # FaFile or other genome object (forge_BSgenome=FALSE)
            genome_sequence <- GTF_annotation\$genome
        } else {
            genome_sequence <- NULL
        }
            genome_sequence <- NULL
        }
        GTF_annotation <<- GTF_annotation
        genome_seq <<- genome_sequence
    }
    assign("load_annotation", patched_load_annotation, envir = ns)
    lockBinding("load_annotation", ns)
    cat("Patched ORFquant::load_annotation to handle NULL genome_package\\n")
}
fix_load_annotation()

# Read the rl_cutoff file
rl_cutoff_data <- read.table(
    "${rl_cutoff}",
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE
)

cat("Read rl_cutoff file with", nrow(rl_cutoff_data), "entries\\n")
print(rl_cutoff_data)

# Create the path to rl_cutoff file for prepare_for_ORFquant
rl_cutoff_file <- "${rl_cutoff}"

# Run prepare_for_ORFquant with the corrected offsets
prepare_for_ORFquant(
    annotation_file = "${annotation}",
    bam_file = "${bam}",
    path_to_rl_cutoff_file = rl_cutoff_file,
    chunk_size = ${chunk_size},
    dest_name = "${prefix}_corrected"
)

cat("Successfully generated corrected for_ORFquant file\\n")

# Write versions
writeLines(
    c(
        '"${task.process}":',
        paste0('    orfquant: "', packageVersion("ORFquant"), '"'),
        paste0('    r-base: "', R.version[["major"]], ".", R.version[["minor"]], '"')
    ),
    "versions.yml"
)
RSCRIPTEOF

    # Use Rscript from the Conda environment if available
    if [[ -n "\$CONDA_PREFIX" ]]; then
        "\$CONDA_PREFIX/bin/Rscript" script.R
    else
        Rscript script.R
    fi
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_corrected_for_ORFquant

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        orfquant: "1.02"
        r-base: "4.1"
    END_VERSIONS
    """
}
