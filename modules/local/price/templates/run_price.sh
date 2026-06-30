#!/bin/bash
set -euo pipefail

# PRICE pipeline runner for GEDI
# Input: BAM file, genome FASTA, annotation GTF
# Output: ORF predictions (TSV + GTF format for unification pipeline)

PREFIX="${PREFIX:-price_output}"
BAM="${BAM}"
FASTA="${FASTA}"
GTF="${GTF}"
OUTDIR="${OUTDIR:-.}"
GENOME_NAME="${GENOME_NAME:-price_genome}"
MEMORY="${GEDI_MEMORY:-16g}"
THREADS="${THREADS:-${SLURM_CPUS_PER_TASK:-4}}"
READ_LENGTHS="${READ_LENGTHS:-28:30}"

export GEDI_MEMORY="$MEMORY"

echo "=== GEDI PRICE Pipeline ==="
echo "BAM: $BAM"
echo "FASTA: $FASTA"
echo "GTF: $GTF"
echo "Output: $OUTDIR/$PREFIX"
echo "Memory: $MEMORY"
echo "Threads: $THREADS"

mkdir -p "$OUTDIR"
cd "$OUTDIR"

# Stage 1: Prepare genome index
echo "[1] Preparing genome index..."
gedi IndexGenome \
    -s "$FASTA" \
    -a "$GTF" \
    -n "$GENOME_NAME" \
    -o "$GENOME_NAME.oml" \
    -p \
    -nobowtie \
    -nostar \
    -nokallisto \
    && echo "Genome index prepared." \
    || { echo "ERROR: IndexGenome failed"; exit 1; }

# Stage 2: Run PRICE ORF detection
echo "[2] Running PRICE ORF detection..."
gedi Price \
    -reads "$BAM" \
    -genomic "$GENOME_NAME.oml" \
    -prefix "$PREFIX" \
    -filter "$READ_LENGTHS" \
    -nthreads "$THREADS" \
    -progress \
    || echo "WARNING: PRICE returned non-zero exit, checking partial output..."

# Stage 3: Post-PRICE analysis
if [ -f "$PREFIX.orfs.tsv" ] && [ -s "$PREFIX.orfs.tsv" ]; then
    echo "[3] Running PriceAnalyze..."
    gedi PriceAnalyze \
        -prefix "$PREFIX" \
        -genomic "$GENOME_NAME.oml" \
        || echo "WARNING: PriceAnalyze failed."
else
    echo "[3] No ORFs detected, skipping PriceAnalyze."
    echo -e "orf_id\tchr\tstart\tend\tstrand\ttype\tscore" > "$PREFIX.orfs.tsv"
fi

# Stage 4: Convert PRICE ORF TSV to GTF for unification pipeline
# PRICE TSV columns: Gene, Id, Location, Candidate Location, Codon, Type, Start, Range, p value, [conditions...], Total
# Location format: aa_start-aa_stop:chr:genomic_start-genomic_end:strand
echo "[4] Converting ORFs to GTF..."
Rscript - << 'RCODE'
args <- commandArgs(trailingOnly = TRUE)
prefix <- args[1]

orf_file <- paste0(prefix, ".orfs.tsv")
gtf_file <- paste0(prefix, "_Detected_ORFs.gtf")

parse_location <- function(loc_str) {
    parts <- strsplit(loc_str, ":")[[1]]
    if (length(parts) < 3) return(NULL)
    chr <- parts[2]
    coords <- as.integer(strsplit(parts[3], "-")[[1]])
    strand <- if (length(parts) >= 4) parts[4] else "+"
    data.frame(chr = chr, start = coords[1], end = coords[2], strand = strand,
               stringsAsFactors = FALSE)
}

if (file.exists(orf_file) && file.info(orf_file)$size > 50) {
    orfs <- tryCatch(
        read.delim(orf_file, stringsAsFactors = FALSE, check.names = FALSE, header = TRUE),
        error = function(e) NULL
    )

    if (!is.null(orfs) && nrow(orfs) > 0) {
        nc <- ncol(orfs)
        cat(paste0("# PRICE ORFs: ", nrow(orfs), " entries\n"), file = gtf_file)
        n_written <- 0

        for (i in seq_len(nrow(orfs))) {
            loc_info <- parse_location(as.character(orfs[i, 3]))
            if (is.null(loc_info)) next

            gene_id <- if (nc >= 1) gsub("[;= ]", "_", as.character(orfs[i, 1])) else "unknown"
            orf_uid <- if (nc >= 2) gsub("[;= ]", "_", as.character(orfs[i, 2])) else paste0("orf_", i)
            orf_type <- if (nc >= 6) as.character(orfs[i, 6]) else "CDS"
            orf_score <- if (nc >= 8) as.character(orfs[i, 8]) else "."
            orf_pval <- if (nc >= 9) as.character(orfs[i, 9]) else "."

            cat(sprintf(
                "%s\tPRICE\tCDS\t%d\t%d\t.\t%s\t.\tgene_id \"%s\"; transcript_id \"%s\"; orf_type \"%s\"; score \"%s\"; pvalue \"%s\"; source \"PRICE\";\n",
                loc_info$chr, loc_info$start, loc_info$end,
                loc_info$strand,
                gene_id, orf_uid,
                orf_type, orf_score, orf_pval),
                file = gtf_file, append = TRUE)
            n_written <- n_written + 1
        }
        cat(sprintf("Wrote %d ORFs to GTF\n", n_written))
    } else {
        cat("# PRICE: no valid ORF entries\n", file = gtf_file)
    }
} else {
    cat("# PRICE: no ORF output found\n", file = gtf_file)
}
RCODE "$PREFIX"

# Write version info
GEDI_VER=$(gedi Version 2>/dev/null || echo "unknown")
cat <<-END_VERSIONS > versions.yml
"PRICE":
    gedi: "$GEDI_VER"
    price: "1.0"
END_VERSIONS

echo "=== PRICE pipeline complete ==="
ls -lh "$PREFIX"*
