process RIBOWALTZ {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        (params.ribowaltz_container ?: 'oras://community.wave.seqera.io/library/ribowaltz_bioconductor-genomicfeatures_r-data.table:latest') :
        (params.ribowaltz_container ?: 'community.wave.seqera.io/library/ribowaltz_bioconductor-genomicfeatures_r-data.table:latest') }"

    input:
    tuple val(meta), path(bam), path(bai)
    path gtf
    path fasta

    output:
    tuple val(meta), path("*_psite_offset.tsv")       , emit: psite_offset
    tuple val(meta), path("*_psite_offset.txt")       , emit: psite_offset_txt
    tuple val(meta), path("*_cds_coverage.tsv")       , emit: cds_coverage, optional: true
    tuple val(meta), path("*_codon_usage.tsv")        , emit: codon_usage, optional: true
    tuple val(meta), path("*_frame_distribution.tsv") , emit: frame_distribution, optional: true
    tuple val(meta), path("*_region_distribution.tsv"), emit: region_distribution, optional: true
    tuple val(meta), path("*_ribowaltz_plots")        , emit: plots
    path "versions.yml"                               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def read_lengths = params.ribowaltz_read_lengths ?: [28, 29, 30]
    // Convert Groovy list to R vector string: [28, 29, 30] -> "c(28, 29, 30)"
    def read_lengths_r = "c(${read_lengths.join(', ')})"
    def r_script = file("${moduleDir}/templates/run_ribowaltz.R").text
        .replace('${prefix}', prefix)
        .replace('${bam}', bam instanceof Path ? bam.name : bam.toString())
        .replace('${gtf}', gtf instanceof Path ? gtf.name : gtf.toString())
        .replace('${read_lengths_r}', read_lengths_r)
        .replace('${task.process}', task.process.toString())
    """
    #!/bin/bash
    set -euo pipefail

    # Install dependencies if not available
    cat <<'INSTALL_SCRIPT' > install_ribowaltz.R
    # Install riboWaltz
    if (!requireNamespace("riboWaltz", quietly = TRUE)) {
        cat("Installing riboWaltz...\\n")
        if (file.exists("${workflow.projectDir}/patched_packages/riboWaltz")) {
            install.packages("${workflow.projectDir}/patched_packages/riboWaltz",
                repos = NULL, type = "source")
        } else {
            if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
            remotes::install_github("LabTranslationalArchitectomics/riboWaltz", upgrade = "never")
        }
    }
    cat("riboWaltz version:", as.character(packageVersion("riboWaltz")), "\\n")

    # Install txdbmaker for Bioc 3.20+ compatibility
    if (!requireNamespace("txdbmaker", quietly = TRUE)) {
        cat("Installing txdbmaker for Bioc 3.20+ compatibility...\\n")
        if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
        BiocManager::install("txdbmaker", update = FALSE, ask = FALSE)
        cat("txdbmaker version:", as.character(packageVersion("txdbmaker")), "\\n")
    }
INSTALL_SCRIPT

    # Run installation
    Rscript install_ribowaltz.R

    # Run main analysis
    # Note: template variables are expanded manually because template()
    # in Nextflow 26.x returns a file path instead of content.
    cat <<'RSCRIPT' > script.R
${r_script}
RSCRIPT

    Rscript script.R
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_psite_offset.tsv
    touch ${prefix}_psite_offset.txt
    touch ${prefix}_cds_coverage.tsv
    touch ${prefix}_codon_usage.tsv
    touch ${prefix}_frame_distribution.tsv
    touch ${prefix}_region_distribution.tsv
    mkdir -p ${prefix}_ribowaltz_plots

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ribowaltz: "2.0"
        r-data.table: "1.14.8"
        r-ggplot2: "3.4.4"
    END_VERSIONS
    """
}
