process RIBOSEQC_PREPAREANNOTATION {
    tag "${gtf.baseName}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    // Use ORFquant patched container if available (has RiboseQC + Bioc 3.20 fixes)
    container "${ params.orfquant_container ?:
        (workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
            'https://depot.galaxyproject.org/singularity/riboseqc:1.1--r36_1' :
            'quay.io/biocontainers/riboseqc:1.1--r36_1') }"

    input:
    path gtf
    path fasta

    output:
    path "*_Rannot"              , emit: annotation
    path "BSgenome.*"            , emit: bsgenome, optional: true
    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = gtf.baseName
    """
    #!/bin/bash

    cat <<'RSCRIPT' > script.R
    # Install RiboseQC from local source if not already available
    # (patched container has source at /opt/riboseqc/RiboseQC with Bioc 3.20 fixes)
    rlibs_local <- file.path(getwd(), "rlibs")
    dir.create(rlibs_local, showWarnings = FALSE, recursive = TRUE)
    .libPaths(c(rlibs_local, .libPaths()))
    local_src <- "/opt/riboseqc/RiboseQC"
    if (dir.exists(local_src) && !requireNamespace("RiboseQC", quietly=TRUE)) {
        cat("Installing RiboseQC from local source...\\n")
        install.packages(local_src, repos=NULL, type="source", lib=rlibs_local, quiet=TRUE)
    }
    library(RiboseQC)
    library(Biostrings)
    library(rtracklayer)

    # Build 2bit file from FASTA (required for forge_BSgenome=TRUE)
    cat("Building 2bit file from FASTA...\\n")
    twobit_path <- "genome.2bit"
    if (!file.exists(twobit_path)) {
        genome <- readDNAStringSet("${fasta}")
        genome <- replaceAmbiguities(genome, new = "N")
        export(genome, twobit_path, format = "2bit")
        cat("2bit file created:", twobit_path, "\\n")
    }

    # Prepare annotation files with BSgenome forge
    # RiboseQC 0.99.0: no genome_seq parameter, uses twobit_file directly
    cat("Preparing annotation with BSgenome forge...\\n")
    prepare_annotation_files(
        annotation_directory = ".",
        twobit_file = twobit_path,
        gtf_file = "${gtf}",
        scientific_name = "Genome.annotation",
        annotation_name = "custom",
        export_bed_tables_TxDb = FALSE,
        forge_BSgenome = TRUE,
        create_TxDb = TRUE
    )

    # Write versions
    writeLines(
        c(
            '"${task.process}":',
            paste0('    riboseqc: "', packageVersion("RiboseQC"), '"')
        ),
        "versions.yml"
    )

    cat("Annotation prepared successfully\\n")
    RSCRIPT

    Rscript script.R
    """

    stub:
    def prefix = gtf.baseName
    """
    touch ${prefix}_Rannot
    mkdir -p BSgenome.Genome.annotation.custom

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        riboseqc: "1.1"
    END_VERSIONS
    """
}
