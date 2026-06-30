#!/usr/bin/env bash
set -euo pipefail

# Mirrors: modules/nf-core/ribotricer/prepareorfs/main.nf
IMG_URL="https://depot.galaxyproject.org/singularity/ribotricer:1.3.3--pyhdfd78af_0"

usage() {
  cat <<'EOF'
Usage:
  07_ribotricer_prepareorfs.sh --prefix NAME --gtf annot.gtf --fasta genome.fa [--outdir DIR]

Required:
  --prefix output prefix
  --gtf    GTF
  --fasta  genome FASTA

Options:
  --outdir output directory (default: ./out_ribotricer_index)

Env:
  BIND_EXTRA extra singularity binds, comma-separated (e.g. /mnt:/mnt)
EOF
}

PREFIX=""
GTF=""
FASTA=""
OUTDIR="./out_ribotricer_index"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --prefix) PREFIX="$2"; shift 2;;
    --gtf) GTF="$2"; shift 2;;
    --fasta) FASTA="$2"; shift 2;;
    --outdir) OUTDIR="$2"; shift 2;;
    -h|--help) usage; exit 0;;
    *) echo "Unknown arg: $1"; usage; exit 2;;
  esac
done

if [[ -z "$PREFIX" || -z "$GTF" || -z "$FASTA" ]]; then
  usage
  exit 2
fi

mkdir -p "$OUTDIR" "./containers"
OUTDIR="$(cd "$OUTDIR" && pwd)"
WORKDIR="$(pwd)"

# Auto-detect and bind mount input file directories
auto_bind_paths() {
  local gtf_path="$1"
  local fasta_path="$2"
  local bind_paths=""

  # Convert to absolute paths
  gtf_abs="$(cd "$(dirname "$gtf_path")" && pwd)/$(basename "$gtf_path")"
  fasta_abs="$(cd "$(dirname "$fasta_path")" && pwd)/$(basename "$fasta_path")"

  # Extract parent directories
  gtf_dir="$(dirname "$gtf_abs")"
  fasta_dir="$(dirname "$fasta_abs")"

  # Add unique directories to bind list
  for dir in "$gtf_dir" "$fasta_dir" "$OUTDIR"; do
    if [[ ":$bind_paths:" != *":$dir:"* ]]; then
      bind_paths="${bind_paths:+$bind_paths,}$dir:$dir"
    fi
  done

  echo "$bind_paths"
}

# Auto-detect bind mounts
AUTO_BINDS="$(auto_bind_paths "$GTF" "$FASTA")"
echo "[INFO] Auto-detected bind mounts: $AUTO_BINDS"

# Convert GTF and FASTA to absolute paths for container
GTF="$(cd "$(dirname "$GTF")" && pwd)/$(basename "$GTF")"
FASTA="$(cd "$(dirname "$FASTA")" && pwd)/$(basename "$FASTA")"

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

IMG="$(pull_img "$IMG_URL")"

# Combine auto-detected binds with BIND_EXTRA
ALL_BINDS="$WORKDIR:$WORKDIR,$AUTO_BINDS${BIND_EXTRA:+,$BIND_EXTRA}"

echo "[INFO] Final bind mounts: $ALL_BINDS"

singularity exec \
  --bind "$ALL_BINDS" \
  --pwd "$WORKDIR" \
  "$IMG" \
  bash -lc "
set -euo pipefail
cd '$OUTDIR'

ribotricer prepare-orfs \
  --gtf '$GTF' \
  --fasta '$FASTA' \
  --prefix '$PREFIX'

ribotricer --version 2>&1 | head -n 1 > versions.ribotricer.txt
"

echo "[OK] Output: $OUTDIR/${PREFIX}_candidate_orfs.tsv"
