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
    // PRICE read-length filter: by default use 28:30 (canonical RPFs).
    // Use != null check (not Elvis ?:) because Groovy treats 0 as falsy,
    // so `params.sorf_read_len_min ?: 28` would silently convert an explicit
    // `--sorf_read_len_min 0` (user wants "no filter") into the default 28.
    // When both min/max are 0, pass an empty filter so PRICE uses all reads.
    def length_min = params.sorf_read_len_min != null ? params.sorf_read_len_min : 28
    def length_max = params.sorf_read_len_max != null ? params.sorf_read_len_max : 30
    def price_filter = (length_min == 0 && length_max == 0) ? '' : "${length_min}:${length_max}"
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
    gedi Price \
        -reads ${bam} \
        -genomic ${genome_name}.oml \
        -prefix ${prefix} \
        ${price_filter ? "-filter ${price_filter}" : ''} \
        -nthreads ${task.cpus} \
        -progress \
        || {
            echo "WARNING: PRICE main pipeline returned non-zero exit code."
            echo "Some stages may have completed. Checking for partial output..."
        }

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
    """
    echo -e "orf_id\\tchr\\tstart\\tend\\tstrand\\ttype\\tscore" > ${prefix}.orfs.tsv
    echo "# PRICE stub" > ${prefix}_Detected_ORFs.gtf

    printf '"%s":\\n    gedi: "stub"\\n    price: "stub"\\n' \\
        '${task.process}' > versions.yml
    """
}
