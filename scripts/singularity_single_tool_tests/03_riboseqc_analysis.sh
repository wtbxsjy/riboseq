#!/usr/bin/env bash
set -euo pipefail

# Mirrors: modules/local/riboseqc/analysis/main.nf
IMG_URL="https://depot.galaxyproject.org/singularity/riboseqc:1.1--r36_1"
SAMTOOLS_IMG_URL="https://depot.galaxyproject.org/singularity/samtools:1.21--h50ea8bc_0"

usage() {
  cat <<'EOF'
Usage:
  03_riboseqc_analysis.sh --sample ID --bam in.bam --annotation X_Rannot --fasta genome.fa [--outdir DIR]

Required:
  --sample      sample ID (prefix)
  --bam         input BAM (sorted + indexed recommended)
  --annotation  RiboseQC annotation file (*_Rannot)
  --fasta       genome FASTA

Options:
  --outdir    output directory (default: ./out_riboseqc_analysis)
  --fast-mode TRUE|FALSE (default: TRUE)
  --cpus      threads for samtools index (default: 4)

Env:
  BIND_EXTRA extra singularity binds, comma-separated (e.g. /mnt:/mnt)
EOF
}

SAMPLE=""
BAM=""
ANNOT=""
FASTA=""
OUTDIR="./out_riboseqc_analysis"
FAST_MODE="TRUE"
CPUS=4

while [[ $# -gt 0 ]]; do
  case "$1" in
    --sample) SAMPLE="$2"; shift 2;;
    --bam) BAM="$2"; shift 2;;
    --annotation) ANNOT="$2"; shift 2;;
    --fasta) FASTA="$2"; shift 2;;
    --outdir) OUTDIR="$2"; shift 2;;
    --fast-mode) FAST_MODE="$2"; shift 2;;
    --cpus) CPUS="$2"; shift 2;;
    -h|--help) usage; exit 0;;
    *) echo "Unknown arg: $1"; usage; exit 2;;
  esac
done

if [[ -z "$SAMPLE" || -z "$BAM" || -z "$ANNOT" || -z "$FASTA" ]]; then
  usage
  exit 2
fi

mkdir -p "$OUTDIR" "./containers"
OUTDIR="$(cd "$OUTDIR" && pwd)"
WORKDIR="$(pwd)"

pull_img() {
  local url="$1"
  local base
  base="$(basename "$url")"
  base="${base//:/_}"
  local sif="$(pwd)/containers/${base}.sif"
  if [[ ! -f "$sif" ]]; then
    singularity pull --disable-cache --force "$sif" "$url"
  fi
  echo "$sif"
}

ensure_bai() {
  local bam="$1"
  local img="$2"
  if [[ -f "${bam}.bai" || -f "${bam%.bam}.bai" || -f "${bam}.csi" ]]; then
    return 0
  fi
  echo "[INFO] Missing BAM index; creating with samtools index"
  singularity exec \
    --bind "$WORKDIR:$WORKDIR${BIND_EXTRA:+,$BIND_EXTRA}" \
    --pwd "$WORKDIR" \
    "$img" \
    samtools index -@ "$CPUS" "$bam"
}

IMG="$(pull_img "$IMG_URL")"
SAMTOOLS_IMG="$(pull_img "$SAMTOOLS_IMG_URL")"

ensure_bai "$BAM" "$SAMTOOLS_IMG"

singularity exec \
  --bind "$WORKDIR:$WORKDIR${BIND_EXTRA:+,$BIND_EXTRA}" \
  --pwd "$WORKDIR" \
  "$IMG" \
  bash -lc "
set -euo pipefail
cd '$OUTDIR'

cat <<'RSCRIPT' > script.R
library(RiboseQC)
library(Rsamtools)

RiboseQC_analysis(
  annotation_file = \"$ANNOT\",
  bam_files = \"$BAM\",
  genome_seq = \"$FASTA\",
  dest_names = \"$SAMPLE\",
  sample_names = \"$SAMPLE\",
  fast_mode = ${FAST_MODE},
  create_report = FALSE,
  write_tmp_files = TRUE
)

writeLines(
  c(
    'RIBOSEQC_ANALYSIS:',
    paste0('    riboseqc: \"', packageVersion('RiboseQC'), '\"')
  ),
  'versions.yml'
)
RSCRIPT

Rscript script.R

# Convert P-sites bedgraphs to ggRibo TSV (same as module)
if [ -f \"${SAMPLE}_P_sites_plus.bedgraph\" ]; then
  awk -v OFS='\t' '{print $4, $1, $3, "+"}' \"${SAMPLE}_P_sites_plus.bedgraph\" > \"${SAMPLE}_ggribo.tsv\"
fi
if [ -f \"${SAMPLE}_P_sites_minus.bedgraph\" ]; then
  awk -v OFS='\t' '{print $4, $1, $3, "-"}' \"${SAMPLE}_P_sites_minus.bedgraph\" >> \"${SAMPLE}_ggribo.tsv\"
fi
"

echo "[OK] RiboseQC outputs in: $OUTDIR (look for ${SAMPLE}_for_ORFquant)"
