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

    output:
    tuple val(meta), path("*_for_ORFquant_corrected"), emit: for_orfquant
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

    cat <<'RSCRIPTEOF' > script.R
# Load ORFquant library
library(ORFquant)

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
    touch ${prefix}_for_ORFquant_corrected

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        orfquant: "1.02"
        r-base: "4.1"
    END_VERSIONS
    """
}
