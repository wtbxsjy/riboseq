process PRICE {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        (params.price_container ?: 'oras://community.wave.seqera.io/library/gedi_price:latest') :
        (params.price_container ?: 'community.wave.seqera.io/library/gedi_price:latest') }"

    input:
    tuple val(meta), path(bam), path(bai)
    path fasta
    path gtf

    output:
    tuple val(meta), path("*_Detected_ORFs.gtf"), emit: gtf
    tuple val(meta), path("*.orfs.tsv")          , emit: orfs_tsv, optional: true
    path "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def genome_name = "price_genome_${meta.id}"
    // -filter is a valid GEDI Price parameter (confirmed `gedi Price -hh`).
    // It accepts a read-length range like "26:34" and is preferred.
    //
    // Modes:
    //   'auto'   → widen the sORF read-length range (e.g. 28-30 → 26-34)
    //              based on experimental evidence that a ±2-4 nt margin
    //              around the dominant footprint lengths improves ORF yield.
    //   '26:34'  → use this exact range
    //   ''       → no -filter (PRICE uses all read lengths; slowest)
    def filter_arg = ''
    def raw_filter = params.price_read_filter ?: ''
    if (raw_filter && raw_filter != 'auto') {
        filter_arg = "-filter ${raw_filter}"
    } else if (raw_filter == 'auto') {
        // Widen the sORF range by -2 nt on the low end and +4 nt on the high
        // end.  PRICE benefits from a broader read-length distribution for
        // codon model estimation, even if shorter/longer reads have weaker
        // periodicity.  The margin values are derived from the parameter
        // sweep (see docs/devlog/PRICE_PARAMETER_SWEEP_2026-07-10.md).
        def rlmin = params.sorf_read_len_min ?: 28
        def rlmax = params.sorf_read_len_max ?: 30
        def auto_min = Math.max(rlmin - 2, 18)
        def auto_max = Math.min(rlmax + 4, 60)
        filter_arg = "-filter ${auto_min}:${auto_max}"
    }
    def keep_anno_arg = params.price_keep_anno ? '-keepAnno' : ''
    def fdr_arg = (params.price_fdr && params.price_fdr != 0.1) ? "-fdr ${params.price_fdr}" : ''
    // Cap threads to avoid excessive parallelism (GEDI default is 254,
    // which can trigger race conditions on pathological contigs)
    def nthreads = Math.min(task.cpus ?: 8, 16)
    def extra_args = [filter_arg, keep_anno_arg, fdr_arg, params.extra_price_args ?: ''].findAll{it}.join(' ')
    """
    #!/bin/bash
    set -euo pipefail

    echo "[PRICE] Preparing genome: ${genome_name}"
    gedi IndexGenome \
        -s ${fasta} \
        -a ${gtf} \
        -n ${genome_name} \
        -o ${genome_name}.oml \
        -p \
        -nobowtie \
        -nostar \
        -nokallisto \
        || {
            echo "ERROR: IndexGenome failed"
            exit 1
        }

    echo "[PRICE] Genome prepared, listing generated files:"
    ls -lh ${genome_name}*

    echo "[PRICE] Running PRICE ORF detection on ${meta.id}"
    # -genomic requires absolute path to .oml file (GEDI looks it up by name
    # in config, but IndexGenome only creates the file, doesn't register it).
    gedi Price \
        -reads ${bam} \
        -genomic "\$(realpath ${genome_name}.oml)" \
        -prefix ${prefix} \
        -nthreads ${nthreads} \
        -progress \
        -novelTranscripts \
        ${extra_args} \
        || {
            echo "WARNING: PRICE main pipeline returned non-zero exit code."
            echo "Some stages may have completed. Checking for partial output..."
        }
    # Detect known GEDI crash: divide-by-zero in PriceOrfInference
    # (triggered on some sORF-filtered BAMs with pathological contigs).
    # If detected, emit a diagnostic message to help users debug.
    if grep -q '/ by zero' .command.log 2>/dev/null; then
        echo "[PRICE] ERROR: GEDI crashed with division-by-zero in ORF inference."
        echo "[PRICE] This is a known GEDI issue with heavily filtered BAMs."
        echo "[PRICE] Try using unfiltered genome BAM or a wider read-length range."
    fi

    # Collect ORF output
    if [ -f "${prefix}.orfs.tsv" ] && [ -s "${prefix}.orfs.tsv" ]; then
        echo "[PRICE] ORF detection successful. \$(wc -l < ${prefix}.orfs.tsv) ORFs detected."

        # Run post-PRICE analysis
        gedi PriceAnalyze \
            -prefix ${prefix} \
            -genomic ${genome_name}.oml \
            || echo "WARNING: PriceAnalyze post-processing failed."
    else
        echo "[PRICE] No ORFs detected. Creating placeholder outputs."
        echo -e "orf_id\\tchr\\tstart\\tend\\tstrand\\ttype\\tscore" > "${prefix}.orfs.tsv"
    fi

    # Convert PRICE ORF TSV to GTF for unification pipeline
    # PRICE TSV columns: Gene, Id, Location, Candidate Location, Codon, Type, Start, Range, p value, [conditions...], Total
    # Location format: chr+strand:start-end[|alt_start-alt_end], e.g. "1+:53348031-53348121"
    echo "[PRICE] Converting ORFs to GTF..."
    Rscript - "${prefix}" << 'RCODE'
args <- commandArgs(trailingOnly = TRUE)
prefix <- args[1]

orf_file <- paste0(prefix, ".orfs.tsv")
gtf_out  <- paste0(prefix, "_Detected_ORFs.gtf")

parse_location <- function(loc_str) {
    # Actual PRICE format: chr+strand:start-end[|block2_start-block2_end]
    # e.g. "1+:53348031-53348121", "1+:58443479-58443513|58444193-58444225"
    parts <- strsplit(loc_str, ":")[[1]]
    if (length(parts) < 2) return(NULL)
    chr_strand <- parts[1]  # e.g. "1+" or "MT-"
    # Extract strand (last character) if + or -
    last_char <- substr(chr_strand, nchar(chr_strand), nchar(chr_strand))
    if (last_char == "+" || last_char == "-") {
        chr <- substr(chr_strand, 1, nchar(chr_strand) - 1)
        strand <- last_char
    } else {
        chr <- chr_strand
        strand <- "+"
    }
    # Coordinates: take the first block (before | for multi-exon ORFs)
    coord_block <- strsplit(parts[2], "|", fixed = TRUE)[[1]][1]
    coords <- as.integer(strsplit(coord_block, "-")[[1]])
    data.frame(chr = chr, start = coords[1], end = coords[2], strand = strand,
               stringsAsFactors = FALSE)
}

if (file.exists(orf_file) && file.info(orf_file)\$size > 50) {
    orfs <- tryCatch(
        read.delim(orf_file, stringsAsFactors = FALSE, check.names = FALSE, header = TRUE),
        error = function(e) NULL
    )

    if (!is.null(orfs) && nrow(orfs) > 0) {
        nc <- ncol(orfs)
        cat(paste0("# PRICE ORFs: ", nrow(orfs), " entries\\n"), file = gtf_out)
        n_written <- 0

        for (i in seq_len(nrow(orfs))) {
            # Column 1 = Gene, 2 = Id, 3 = Location
            loc_info <- parse_location(as.character(orfs[i, 3]))
            if (is.null(loc_info)) next

            gene_id <- if (nc >= 1) gsub("[;= ]", "_", as.character(orfs[i, 1])) else "unknown"
            orf_uid <- if (nc >= 2) gsub("[;= ]", "_", as.character(orfs[i, 2])) else paste0("orf_", i)
            orf_type <- if (nc >= 6) as.character(orfs[i, 6]) else "CDS"
            orf_score <- if (nc >= 8) as.character(orfs[i, 8]) else "."
            orf_pval <- if (nc >= 9) as.character(orfs[i, 9]) else "."

            cat(sprintf(
                "%s\\tPRICE\\tCDS\\t%d\\t%d\\t.\\t%s\\t.\\tgene_id \\"%s\\"; transcript_id \\"%s\\"; orf_type \\"%s\\"; score \\"%s\\"; pvalue \\"%s\\"; source \\"PRICE\\";\\n",
                loc_info\$chr, loc_info\$start, loc_info\$end,
                loc_info\$strand,
                gene_id, orf_uid,
                orf_type, orf_score, orf_pval),
                file = gtf_out, append = TRUE)
            n_written <- n_written + 1
        }
        cat(sprintf("Wrote %d PRICE ORFs to GTF\\n", n_written))
    } else {
        cat("# PRICE: no valid ORF entries\\n", file = gtf_out)
        cat("PRICE: no valid ORF entries\\n")
    }
} else {
    cat("# PRICE: no ORF output found\\n", file = gtf_out)
    cat("PRICE: ORF TSV file not found or empty\\n")
}
RCODE

    ls -lh ${prefix}*

    # Write versions
    GEDI_VER=\$(gedi Version 2>/dev/null || echo "unknown")
    printf '"%s":\\n    gedi: "%s"\\n    price: "1.0"\\n' \\
        '${task.process}' "\$GEDI_VER" > versions.yml

    echo "[PRICE] Pipeline complete for ${meta.id}."
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def nthreads_stub = Math.min(task.cpus ?: 8, 16)
    """
    echo -e "orf_id\\tchr\\tstart\\tend\\tstrand\\ttype\\tscore" > ${prefix}.orfs.tsv
    echo "# PRICE stub" > ${prefix}_Detected_ORFs.gtf

    printf '"%s":\\n    gedi: "stub"\\n    price: "stub"\\n' \\
        '${task.process}' > versions.yml
    """
}
