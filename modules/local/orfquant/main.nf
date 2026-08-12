process ORFQUANT_RUN {
    tag "$meta.id"
    label 'process_high'

    // Dual-mode: Nextflow picks the right backend based on active profile
    //   -profile conda      → uses pre-built conda environment
    //   -profile singularity → uses SIF container (portable, self-contained)
    conda "${params.orfquant_conda_env ?: "${moduleDir}/environment.yml"}"
    container "${params.orfquant_mirai_container}"

    input:
    tuple val(meta), path(for_orfquant)   // *_for_ORFquant file from RiboseQC
    path annotation                        // *_Rannot file from RiboseQC/ORFquant annotation
    path fasta                             // Genome fasta file

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
    def n_cores = params.orfquant_mirai_daemons ?: task.cpus ?: 32
    // Parse optional arguments
    def write_gtf = args.contains('write_GTF_file=FALSE') ? 'FALSE' : 'TRUE'
    def write_fasta = args.contains('write_protein_fasta=FALSE') ? 'FALSE' : 'TRUE'
    def write_tmp = args.contains('write_temp_files=FALSE') ? 'FALSE' : 'TRUE'
    def plot_results = args.contains('plot_results=TRUE') ? 'TRUE' : 'FALSE'
    // Parallel backend: "mirai" (socket daemon, default) or "mclapply" (fork, legacy)
    def parallel_backend = params.orfquant_parallel_backend ?: 'mirai'
    """
    # Ensure fasta file is available with the expected name (if it was gzipped)
    if [[ "${fasta}" == *.gz ]]; then
        gunzip -c ${fasta} > \$(basename ${fasta} .gz)
    fi

    # Force re-run: .safe_field() plain-list fix (2026-08-12)
    export ORFQUANT_MIRAI_RERUN=3

    # BLAS thread control — prevent each R process from spawning threads_per_core threads
    export OMP_NUM_THREADS=1
    export OPENBLAS_NUM_THREADS=1
    export MKL_NUM_THREADS=1
    export GOTO_NUM_THREADS=1

    # Write R script — ORFquant, mirai, and all dependencies are pre-installed
    cat > run_orfquant.R <<'RSCRIPTEOF'
    # Suppress BLAS threading inside R (belt + suspenders with shell exports)
    Sys.setenv(OMP_NUM_THREADS="1"); Sys.setenv(OPENBLAS_NUM_THREADS="1")
    Sys.setenv(MKL_NUM_THREADS="1"); Sys.setenv(GOTO_NUM_THREADS="1")
    options(mc.cores = 1)
    tryCatch(BiocParallel::register(BiocParallel::SerialParam()), error = function(e) NULL)

    # Load mirai-optimized parallel backend (disk-backed FaFile streaming)
    source("/opt/orfquant_mirai_optimized.R")
    library(ORFquant)

    # --- Hotfix: patch .safe_field() to handle plain-list ORFs_tx_position ---
    # select_quantify_ORFs() uses lapply() at line 2115 which degrades
    # ORFs_tx_position from GRangesList to a plain list.  The per-gene RDS
    # files contain valid GRanges elements inside a plain list, but the
    # original .safe_field() rejects any plain list as an error guard.
    # This silently discards ALL transcript-space ORFs (ORFs_tx = 0), which
    # forces run_ORFquant() down the genomic-only GTF export path (no attrs).
    #
    # The patch replaces the plain-list guard with a conversion attempt:
    # if the list is non-empty and contains GRanges objects, wrap it in
    # GRangesList() instead of returning an empty GRanges.  Empty lists and
    # unconvertible types still return GRanges() as before.
    cat("[patch] Applying .safe_field() plain-list fix\\n")
    ns <- asNamespace("ORFquant")
    unlockBinding(".safe_field", ns)
    .safe_field_patched <- function(x, field) {
        val <- x[[field]]
        if (is.null(val)) return(GRanges())
        if (is.list(val) && !is(val, "GRanges") && !is(val, "GRangesList") &&
            !is(val, "CompressedGRangesList")) {
            if (length(val) == 0) return(GRanges())
            val <- tryCatch(GRangesList(val), error = function(e) {
                cat("[patch] cannot convert", field, "to GRangesList:",
                    conditionMessage(e), "\\n")
                return(GRanges())
            })
            if (!is(val, "GRangesList") && !is(val, "CompressedGRangesList")) return(GRanges())
        }
        val <- unlist(val)
        if (is.null(val) || length(val) == 0) return(GRanges())
        val
    }
    assign(".safe_field", .safe_field_patched, envir = ns)
    lockBinding(".safe_field", ns)
    unlockBinding(".safe_nested", ns)
    .safe_nested_patched <- function(x, outer, inner) {
        out <- x[[outer]]
        if (is.null(out)) return(GRanges())
        val <- out[[inner]]
        if (is.null(val)) return(GRanges())
        if (is.list(val) && !is(val, "GRanges") && !is(val, "GRangesList") &&
            !is(val, "CompressedGRangesList")) {
            if (length(val) == 0) return(GRanges())
            val <- tryCatch(GRangesList(val), error = function(e) {
                cat("[patch] cannot convert", outer, "/", inner, "to GRangesList:",
                    conditionMessage(e), "\\n")
                return(GRanges())
            })
            if (!is(val, "GRangesList") && !is(val, "CompressedGRangesList")) return(GRanges())
        }
        val <- unlist(val)
        if (is.null(val) || length(val) == 0) return(GRanges())
        val
    }
    assign(".safe_nested", .safe_nested_patched, envir = ns)
    lockBinding(".safe_nested", ns)
    cat("[patch] .safe_field() and .safe_nested() patched\\n")

    # Run ORFquant with error handling for low-quality samples
    cat("Running ORFquant on sample ${prefix}...\\n")
    orfquant_success <- tryCatch({
        run_ORFquant(
            for_ORFquant_file = "${for_orfquant}",
            annotation_file   = "${annotation}",
            n_cores           = ${n_cores},
            prefix            = "${prefix}",
            write_temp_files  = ${write_tmp},
            write_GTF_file    = ${write_gtf},
            write_protein_fasta = ${write_fasta},
            interactive       = FALSE,
            parallel_backend  = "${parallel_backend}"
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
            paste0('    mirai: "', as.character(packageVersion("mirai")), '"'),
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
            HAS_ATTRS=\$(grep -v '^#' ${prefix}_Detected_ORFs.gtf | awk -F'\t' 'NR==1 {print \$9}')
            if [ "\$HAS_ATTRS" = "." ] || [ -z "\$HAS_ATTRS" ]; then
                NEED_FIX=true
            fi
        else
            NEED_FIX=true
        fi
        if [ "\$NEED_FIX" = "true" ]; then
            echo "[GTF fix] Original GTF has no attributes — rebuilding from final_results + FASTA"
            Rscript ${moduleDir}/templates/fix_orfquant_gtf.R ${prefix} 2>&1 || {
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
        orfquant: "1.3.3"
        mirai: "2.7.2"
        r-base: "4.4"
    END_VERSIONS
    """
}
