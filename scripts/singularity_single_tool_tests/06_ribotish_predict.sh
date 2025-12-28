#!/usr/bin/env bash
set -euo pipefail

# Mirrors: modules/nf-core/ribotish/predict/main.nf (per-sample use)
IMG_URL="https://depot.galaxyproject.org/singularity/ribotish:0.2.7--pyhdfd78af_0"
SAMTOOLS_IMG_URL="https://depot.galaxyproject.org/singularity/samtools:1.21--h50ea8bc_0"

usage() {
  cat <<'EOF'
Usage:
  06_ribotish_predict.sh --sample ID --bam in.bam --gtf annot.gtf --fasta genome.fa --ribopara S1.para.py [--cpus N] [--outdir DIR]

Required:
  --sample   sample ID (prefix)
  --bam      ribo BAM (filtered BAM recommended)
  --gtf      GTF
  --fasta    genome FASTA
  --ribopara offset para file from ribotish quality (*.para.py)

Options:
  --cpus     threads (default: 4)
  --outdir   output directory (default: ./out_ribotish_predict)
  --args     extra args passed to ribotish predict (quoted string)

Env:
  BIND_EXTRA extra singularity binds, comma-separated (e.g. /mnt:/mnt)
EOF
}

SAMPLE=""
BAM=""
GTF=""
FASTA=""
RIBOPARA=""
CPUS=4
OUTDIR="./out_ribotish_predict"
EXTRA_ARGS=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --sample) SAMPLE="$2"; shift 2;;
    --bam) BAM="$2"; shift 2;;
    --gtf) GTF="$2"; shift 2;;
    --fasta) FASTA="$2"; shift 2;;
    --ribopara) RIBOPARA="$2"; shift 2;;
    --cpus) CPUS="$2"; shift 2;;
    --outdir) OUTDIR="$2"; shift 2;;
    --args) EXTRA_ARGS="$2"; shift 2;;
    -h|--help) usage; exit 0;;
    *) echo "Unknown arg: $1"; usage; exit 2;;
  esac
done

if [[ -z "$SAMPLE" || -z "$BAM" || -z "$GTF" || -z "$FASTA" || -z "$RIBOPARA" ]]; then
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

ribotish predict \
  -b '$BAM' \
  --ribopara '$RIBOPARA' \
  -f '$FASTA' \
  -g '$GTF' \
  -o '${SAMPLE}_pred.txt' \
  --allresult '${SAMPLE}_all.txt' \
  --transprofile '${SAMPLE}_transprofile.py' \
  -p '$CPUS' \
  $EXTRA_ARGS

ribotish --version > versions.ribotish.txt
"

echo "[OK] Outputs: $OUTDIR/${SAMPLE}_pred.txt"
