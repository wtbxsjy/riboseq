#!/usr/bin/env bash
set -euo pipefail

# Mirrors: modules/nf-core/ribotricer/detectorfs/main.nf
IMG_URL="https://depot.galaxyproject.org/singularity/ribotricer:1.3.3--pyhdfd78af_0"
SAMTOOLS_IMG_URL="https://depot.galaxyproject.org/singularity/samtools:1.21--h50ea8bc_0"

usage() {
  cat <<'EOF'
Usage:
  08_ribotricer_detectorfs.sh --sample ID --bam in.bam --index PREFIX_candidate_orfs.tsv [--stranded forward|reverse|unstranded] [--outdir DIR]

Required:
  --sample sample ID (prefix)
  --bam    input BAM (filtered BAM recommended)
  --index  ribotricer index file (*_candidate_orfs.tsv)

Options:
  --stranded forward|reverse|unstranded (default: forward)
  --outdir   output directory (default: ./out_ribotricer_detect)
  --args     extra args passed to ribotricer detect-orfs (quoted string)

Env:
  BIND_EXTRA extra singularity binds, comma-separated (e.g. /mnt:/mnt)
EOF
}

SAMPLE=""
BAM=""
INDEX=""
STRANDED="forward"
OUTDIR="./out_ribotricer_detect"
EXTRA_ARGS=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --sample) SAMPLE="$2"; shift 2;;
    --bam) BAM="$2"; shift 2;;
    --index) INDEX="$2"; shift 2;;
    --stranded) STRANDED="$2"; shift 2;;
    --outdir) OUTDIR="$2"; shift 2;;
    --args) EXTRA_ARGS="$2"; shift 2;;
    -h|--help) usage; exit 0;;
    *) echo "Unknown arg: $1"; usage; exit 2;;
  esac
done

if [[ -z "$SAMPLE" || -z "$BAM" || -z "$INDEX" ]]; then
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
    samtools index -@ 2 "$bam"
}

IMG="$(pull_img "$IMG_URL")"
SAMTOOLS_IMG="$(pull_img "$SAMTOOLS_IMG_URL")"

ensure_bai "$BAM" "$SAMTOOLS_IMG"

STR_CMD=""
case "$STRANDED" in
  forward) STR_CMD="--stranded yes";;
  reverse) STR_CMD="--stranded reverse";;
  unstranded) STR_CMD="";;
  *) echo "Invalid --stranded: $STRANDED"; exit 2;;
esac

singularity exec \
  --bind "$WORKDIR:$WORKDIR${BIND_EXTRA:+,$BIND_EXTRA}" \
  --pwd "$WORKDIR" \
  "$IMG" \
  bash -lc "
set -euo pipefail
cd '$OUTDIR'

ribotricer detect-orfs \
  --bam '$BAM' \
  --ribotricer_index '$INDEX' \
  --prefix '$SAMPLE' \
  $STR_CMD \
  $EXTRA_ARGS

ribotricer --version 2>&1 | head -n 1 > versions.ribotricer.txt
"

echo "[OK] Output: $OUTDIR/${SAMPLE}_translating_ORFs.tsv"
