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
    tuple val(meta), path("*_Detected_ORFs.gtf.gz")   , emit: gtf, optional: true
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

    # Append a task-local Rlibs directory instead of replacing R_LIBS_USER, so that
    # packages pre-installed in the container (e.g. ORFquant) remain accessible.
    _local_rlibs="${task.workDir}/Rlibs"
    mkdir -p "\$_local_rlibs"
    export R_LIBS_USER="\${_local_rlibs}\${R_LIBS_USER:+:\${R_LIBS_USER}}"

    # Write R script - ORFquant should be pre-installed in custom container
    cat > run_orfquant.R <<'RSCRIPTEOF'
install_orfquant <- function(local_pkg_tgz = NULL, tag = "1.02", local_src = NULL) {
    if (!is.null(local_src) && dir.exists(local_src) && file.exists(file.path(local_src, "DESCRIPTION"))) {
        message("Installing ORFquant from local source: ", local_src)
        cmd <- sprintf("R CMD INSTALL %s", shQuote(local_src))
        status <- system(cmd)
        if (status == 0) return(invisible(TRUE))
    }
    work <- file.path(getwd(), "orfquant_src")
    dir.create(work, showWarnings = FALSE, recursive = TRUE)

    tgz <- local_pkg_tgz
    if (!is.null(tgz) && nzchar(tgz) && file.exists(tgz) && file.info(tgz)[1, "size"] > 0) {
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
        install_orfquant(local_pkg_tgz = local_pkg, tag = "1.02", local_src = "/opt/ORFquant")
    }, error = function(e) {
        stop(
            "ORFquant is not installed and automatic installation failed: ", conditionMessage(e), "\n",
            "Provide a pre-downloaded tarball with --orfquant_pkg (e.g. ORFquant_1.02.0.tar.gz from lcalviell/ORFquant), ",
            "or use a container with ORFquant pre-installed (e.g. --orfquant_container)."
        )
    })

    if (!requireNamespace("ORFquant", quietly = TRUE)) {
        stop("ORFquant install completed but package is still not available on library paths.")
    }
}

# Install txdbmaker if missing (Bioc 3.20+)
if (!requireNamespace("txdbmaker", quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager", repos = "https://cloud.r-project.org", quiet = TRUE)
    BiocManager::install("txdbmaker", update = FALSE, ask = FALSE, quiet = TRUE)
}

	library(ORFquant)

			# load_annotation monkey-patch REMOVED (2026-06-28).
			# The patched ORFquant container now includes fork-safe load_annotation()
			# with FaFile->DNAStringSet conversion.



# Run ORFquant with error handling for low-quality samples
cat("Running ORFquant on sample ${prefix}...\\n")
orfquant_success <- tryCatch({
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
    TRUE
}, error = function(e) {
    error_msg <- conditionMessage(e)
    cat("\\n=== ORFquant Error ===\\n")
    cat(error_msg, "\\n")
    
    # Check for common low-signal/quality errors that should allow pipeline to continue
    is_low_signal_error <- (
        grepl("unable to find an inherited method.*summarizeOverlaps", error_msg, ignore.case = TRUE) ||
        grepl("no method.*coercing.*NULL.*GRanges", error_msg, ignore.case = TRUE) ||
        grepl("summarizeOverlaps.*GRanges.*NULL", error_msg, ignore.case = TRUE) ||
        grepl("Not enough P_sites signal", error_msg, ignore.case = TRUE) ||
        grepl("Not enough P.sites signal", error_msg, ignore.case = TRUE) ||
        grepl("insufficient.*signal", error_msg, ignore.case = TRUE) ||
        grepl("no ORFs? (were |was )?detected", error_msg, ignore.case = TRUE)
    )
    
    if (is_low_signal_error) {
        cat("\\nWARNING: ORFquant failed due to insufficient signal/ORF predictions.\\n")
        cat("This typically occurs when:\\n")
        cat("  - Sample has low ribosome profiling signal\\n")
        cat("  - Not enough P-sites signal over genomic regions\\n")
        cat("  - Very few or no ORFs meet the detection thresholds\\n")
        cat("  - P-site positioning is poor\\n")
        cat("\\nCreating empty output files to allow pipeline continuation...\\n")
        
        # Create empty output files so downstream processes can handle gracefully
        writeLines("# No ORFs detected - insufficient signal", "${prefix}_final_ORFquant_results")
        
        if (${write_gtf}) {
            writeLines("# No ORFs detected", "${prefix}_Detected_ORFs.gtf")
        }
        if (${write_fasta}) {
            writeLines("", "${prefix}_Protein_sequences.fasta")  # Empty FASTA
        }
        if (${write_tmp}) {
            writeLines("# No ORFs detected", "${prefix}_tmp_ORFquant_results")
        }
        
        return(FALSE)
    } else {
        # For other errors, re-throw
        cat("\\nUnexpected ORFquant error. Re-throwing...\\n")
        stop(e)
    }
})

if (orfquant_success) {
    cat("ORFquant completed successfully\\n")
} else {
    cat("ORFquant skipped due to insufficient data\\n")
}

# Optionally generate plots (only if ORFquant succeeded)
if (${plot_results} && orfquant_success) {
    tryCatch({
        plot_ORFquant_results(
            for_ORFquant_file = "${for_orfquant}",
            ORFquant_output_file = paste0("${prefix}", "_final_ORFquant_results"),
            annotation_file = "${annotation}",
            output_plots_path = paste0("${prefix}", "_plots"),
            prefix = "${prefix}"
        )
    }, error = function(e) {
        message("Warning: Could not generate ORFquant plots: ", conditionMessage(e))
    })
} else if (${plot_results} && !orfquant_success) {
    cat("Skipping plot generation - ORFquant did not produce results\\n")
}

# Write versions
writeLines(
    c(
        '"${task.process}":',
        paste0('    orfquant: "', packageVersion("ORFquant"), '"'),
        paste0('    r-base: "', R.Version()[["major"]], ".", R.Version()[["minor"]], '"')
    ),
    "versions.yml"
)
RSCRIPTEOF

    # Run using Rscript
    Rscript run_orfquant.R

    # Fix ORFquant GTF: when ORFs_tx is empty, run_ORFquant() takes the
    # genomic-only export path where ORFs_gen lacks metadata columns (mcols),
    # producing a GTF with empty attributes (column 9 = '.').  When ORFs_tx
    # is non-empty, the transcript-aware path preserves attributes and this
    # fix is skipped to avoid overwriting a correct GTF.
    if [ -f ${prefix}_final_ORFquant_results ] && [ -s ${prefix}_final_ORFquant_results ]; then
        NEED_FIX=false
        if [ -f ${prefix}_Detected_ORFs.gtf ] && [ -s ${prefix}_Detected_ORFs.gtf ]; then
            HAS_ATTRS=\$(grep -v '^#' ${prefix}_Detected_ORFs.gtf | head -1 | awk -F'\t' '{print \$9}')
            if [ "\$HAS_ATTRS" = "." ] || [ -z "\$HAS_ATTRS" ]; then
                NEED_FIX=true
            fi
        else
            NEED_FIX=true
        fi
        if [ "\$NEED_FIX" = "true" ]; then
            echo "[GTF fix] Original GTF has no attributes — rebuilding from final_results + FASTA"
            cat > fix_orfquant_gtf.R <<'FIXGTFEOF'
suppressMessages({
    library(ORFquant)
    library(rtracklayer)
})

args <- commandArgs(trailingOnly = TRUE)
prefix <- args[1]

cat("[GTF fix] Loading ORFquant results...\n")
load(paste0(prefix, "_final_ORFquant_results"))

if (!exists("ORFquant_results") || length(ORFquant_results$ORFs_gen) == 0) {
    cat("[GTF fix] No ORFs found, skipping GTF rewrite\n")
    quit(save = "no", status = 0)
}

g <- ORFquant_results$ORFs_gen

# Parse protein FASTA for per-ORF metadata
fa_file <- paste0(prefix, "_Protein_sequences.fasta")
if (!file.exists(fa_file)) {
    cat("[GTF fix] WARNING: protein FASTA not found, exporting with coordinate-only attributes\n")
    fa_meta <- NULL
} else {
    fa_lines <- readLines(fa_file)
    fa_hdrs <- sub("^>", "", grep("^>", fa_lines, value = TRUE))
    # Format: TRANSCRIPT_START_END|gene_biotype|gene_id|orf_type|ORF_category
    tx_info <- strsplit(fa_hdrs, "\\|")
    fa_meta <- data.frame(
        orf_name  = sapply(tx_info, `[`, 1),
        gene_biotype = sapply(tx_info, `[`, 2),
        gene_id   = sapply(tx_info, `[`, 3),
        orf_type  = sapply(tx_info, `[`, 4),
        orf_category = ifelse(lengths(tx_info) >= 5, sapply(tx_info, `[`, 5), "NA"),
        stringsAsFactors = FALSE
    )
    rownames(fa_meta) <- fa_meta$orf_name
    cat(sprintf("[GTF fix] Parsed %d FASTA entries\n", nrow(fa_meta)))
}

# Extract ORF names from GRanges (multi-exon ORFs share the same name)
orf_names <- names(g)
if (is.null(orf_names)) {
    # Fallback: use sequential IDs
    orf_names <- paste0("ORFquant_", seq_along(g))
    names(g) <- orf_names
}
unique_orfs <- unique(orf_names)
cat(sprintf("[GTF fix] %d CDS features across %d unique ORFs\n", length(g), length(unique_orfs)))

# Build metadata columns: one row per CDS feature
n <- length(g)
g\$type <- "CDS"
mcols(g)$ORF_id    <- orf_names
mcols(g)$gene_id   <- rep("NA", n)
mcols(g)$gene_biotype <- rep("NA", n)
mcols(g)$orf_type  <- rep("NA", n)
mcols(g)$orf_category <- rep("NA", n)

# Fill from FASTA metadata where available
if (!is.null(fa_meta)) {
    matched <- orf_names %in% rownames(fa_meta)
    if (sum(matched) > 0) {
        mcols(g)$gene_id[matched]      <- fa_meta[orf_names[matched], "gene_id"]
        mcols(g)$gene_biotype[matched] <- fa_meta[orf_names[matched], "gene_biotype"]
        mcols(g)$orf_type[matched]     <- fa_meta[orf_names[matched], "orf_type"]
        mcols(g)$orf_category[matched] <- fa_meta[orf_names[matched], "orf_category"]
        cat(sprintf("[GTF fix] Metadata assigned to %d / %d features\n", sum(matched), n))
    }
}

# Set source
mcols(g)$source <- "ORFquant"

# Remove old GTF and export new one
gtf_file <- paste0(prefix, "_Detected_ORFs.gtf")
cat(sprintf("[GTF fix] Writing %d features to %s\n", length(g), gtf_file))
rtracklayer::export(g, gtf_file, format = "gtf")
cat("[GTF fix] Done\n")
FIXGTFEOF
            Rscript fix_orfquant_gtf.R ${prefix} 2>&1 || {
                echo "WARNING: GTF attribute fix failed, keeping original GTF"
            }
        else
            echo "[GTF fix] GTF already has attributes — skipping rebuild"
        fi
    else
        echo "[GTF fix] No ORFquant results file found, skipping GTF fix"
    fi

    # Compress text outputs to save disk space
    for f in *_Detected_ORFs.gtf; do
        [ -f "\$f" ] && gzip -f "\$f" || true
    done
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_final_ORFquant_results
    touch ${prefix}_Detected_ORFs.gtf
    gzip -f ${prefix}_Detected_ORFs.gtf
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
