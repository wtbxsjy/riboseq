#!/usr/bin/env bash
set -euo pipefail

# Mirrors: modules/nf-core/ribotricer/detectorfs/main.nf
IMG_URL="https://depot.galaxyproject.org/singularity/ribotricer:1.3.3--pyhdfd78af_0"
SAMTOOLS_IMG_URL="https://depot.galaxyproject.org/singularity/samtools:1.21--h50ea8bc_0"

usage() {
  cat <<'EOF'
Usage:
  08_ribotricer_detectorfs.sh --sample ID --bam in.bam --index PREFIX_candidate_orfs.tsv [--stranded forward|reverse|unstranded] [--cpus N] [--outdir DIR]

Required:
  --sample sample ID (prefix)
  --bam    input BAM (filtered BAM recommended)
  --index  ribotricer index file (*_candidate_orfs.tsv)

Options:
  --stranded forward|reverse|unstranded (default: forward)
  --cpus     threads for samtools index (default: 4)
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
CPUS=4

while [[ $# -gt 0 ]]; do
  case "$1" in
    --sample) SAMPLE="$2"; shift 2;;
    --bam) BAM="$2"; shift 2;;
    --index) INDEX="$2"; shift 2;;
    --stranded) STRANDED="$2"; shift 2;;
    --cpus) CPUS="$2"; shift 2;;
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

# Auto-detect and bind mount input file directories
auto_bind_paths() {
  local bam_path="$1"
  local index_path="$2"
  local bind_paths=""

  # Convert to absolute paths
  bam_abs="$(cd "$(dirname "$bam_path")" && pwd)/$(basename "$bam_path")"
  index_abs="$(cd "$(dirname "$index_path")" && pwd)/$(basename "$index_path")"

  # Extract parent directories
  bam_dir="$(dirname "$bam_abs")"
  index_dir="$(dirname "$index_abs")"

  # Add unique directories to bind list
  for dir in "$bam_dir" "$index_dir" "$OUTDIR"; do
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
AUTO_BINDS="$(auto_bind_paths "$BAM" "$INDEX")"
echo "[INFO] Auto-detected bind mounts: $AUTO_BINDS"

# Convert BAM and INDEX to absolute paths for container
BAM="$(cd "$(dirname "$BAM")" && pwd)/$(basename "$BAM")"
INDEX="$(cd "$(dirname "$INDEX")" && pwd)/$(basename "$INDEX")"

IMG="$(pull_img "$IMG_URL")"
SAMTOOLS_IMG="$(pull_img "$SAMTOOLS_IMG_URL")"

# Combine auto-detected binds with BIND_EXTRA
ALL_BINDS="$WORKDIR:$WORKDIR,$AUTO_BINDS${BIND_EXTRA:+,$BIND_EXTRA}"

echo "[INFO] Final bind mounts: $ALL_BINDS"

ensure_bai "$BAM" "$SAMTOOLS_IMG" "$ALL_BINDS"

STR_CMD=""
case "$STRANDED" in
  forward) STR_CMD="--stranded yes";;
  reverse) STR_CMD="--stranded reverse";;
  unstranded) STR_CMD="";;
  *) echo "Invalid --stranded: $STRANDED"; exit 2;;
esac

singularity exec \
  --bind "$ALL_BINDS" \
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
