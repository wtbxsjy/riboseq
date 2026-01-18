#!/usr/bin/env bash
set -euo pipefail

# Avoid noisy locale warnings inside the container (host locales may not exist there)
export SINGULARITYENV_LANG=${SINGULARITYENV_LANG:-C}
export SINGULARITYENV_LC_ALL=${SINGULARITYENV_LC_ALL:-C}
export APPTAINERENV_LANG=${APPTAINERENV_LANG:-C}
export APPTAINERENV_LC_ALL=${APPTAINERENV_LC_ALL:-C}

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

# Auto-detect and bind mount input file directories
auto_bind_paths() {
  local bam_path="$1"
  local annot_path="$2"
  local fasta_path="$3"
  local bind_paths=""

  # Convert to absolute paths
  bam_abs="$(cd "$(dirname "$bam_path")" && pwd)/$(basename "$bam_path")"
  annot_abs="$(cd "$(dirname "$annot_path")" && pwd)/$(basename "$annot_path")"
  fasta_abs="$(cd "$(dirname "$fasta_path")" && pwd)/$(basename "$fasta_path")"

  # Extract parent directories
  bam_dir="$(dirname "$bam_abs")"
  annot_dir="$(dirname "$annot_abs")"
  fasta_dir="$(dirname "$fasta_abs")"

  # Add unique directories to bind list
  for dir in "$bam_dir" "$annot_dir" "$fasta_dir" "$OUTDIR"; do
    if [[ ":$bind_paths:" != *":$dir:"* ]]; then
      bind_paths="${bind_paths:+$bind_paths,}$dir:$dir"
    fi
  done

  echo "$bind_paths"
}

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
  local binds="$3"
  if [[ -f "${bam}.bai" || -f "${bam%.bam}.bai" || -f "${bam}.csi" ]]; then
    return 0
  fi
  echo "[INFO] Missing BAM index; creating with samtools index"
  singularity exec \
    --bind "$binds" \
    --pwd "$WORKDIR" \
    "$img" \
    samtools index -@ "$CPUS" "$bam"
}

# Auto-detect bind mounts
AUTO_BINDS="$(auto_bind_paths "$BAM" "$ANNOT" "$FASTA")"
echo "[INFO] Auto-detected bind mounts: $AUTO_BINDS"

# Convert inputs to absolute paths for container
BAM="$(cd "$(dirname "$BAM")" && pwd)/$(basename "$BAM")"
ANNOT="$(cd "$(dirname "$ANNOT")" && pwd)/$(basename "$ANNOT")"
FASTA="$(cd "$(dirname "$FASTA")" && pwd)/$(basename "$FASTA")"

IMG="$(pull_img "$IMG_URL")"
SAMTOOLS_IMG="$(pull_img "$SAMTOOLS_IMG_URL")"

# Combine auto-detected binds with BIND_EXTRA
ALL_BINDS="$WORKDIR:$WORKDIR,$AUTO_BINDS${BIND_EXTRA:+,$BIND_EXTRA}"

echo "[INFO] Final bind mounts: $ALL_BINDS"

ensure_bai "$BAM" "$SAMTOOLS_IMG" "$ALL_BINDS"

singularity exec \
  --bind "$ALL_BINDS" \
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
  awk -v OFS='\t' '{print \$4, \$1, \$3, \"+\"}' \"${SAMPLE}_P_sites_plus.bedgraph\" > \"${SAMPLE}_ggribo.tsv\"
fi
if [ -f \"${SAMPLE}_P_sites_minus.bedgraph\" ]; then
  awk -v OFS='\t' '{print \$4, \$1, \$3, \"-\"}' \"${SAMPLE}_P_sites_minus.bedgraph\" >> \"${SAMPLE}_ggribo.tsv\"
fi
"

echo "[OK] RiboseQC outputs in: $OUTDIR (look for ${SAMPLE}_for_ORFquant)"
