#!/usr/bin/env bash
set -euo pipefail

# Mirrors: modules/local/convert_ribotish_to_gencode/main.nf
# Converts Ribo-TISH predict output to gencode-riboseqORFs compatible format

IMG_URL="https://depot.galaxyproject.org/singularity/biopython:1.81"

usage() {
  cat <<'EOF'
Usage:
  14_ribotish_to_gencode.sh --sample ID --predict ribotish_pred.txt --fasta genome.fa [--outdir DIR]

Required:
  --sample   sample ID (used as study_id in output)
  --predict  Ribo-TISH predict output file (*_pred.txt)
  --fasta    genome FASTA file

Options:
  --outdir      output directory (default: ./out_ribotish_to_gencode)
  --min-length  minimum ORF length in amino acids (default: 16)
  --args        extra args passed to ribotish_to_gencode.py (quoted string)

Env:
  BIND_EXTRA extra singularity binds, comma-separated (e.g. /mnt:/mnt)

Output:
  ${SAMPLE}.gencode.fa   - FASTA in gencode-riboseqORFs format
  ${SAMPLE}.gencode.bed  - BED (1-based) in gencode-riboseqORFs format
EOF
}

SAMPLE=""
PREDICT=""
FASTA=""
OUTDIR="./out_ribotish_to_gencode"
MIN_LENGTH=16
EXTRA_ARGS=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --sample) SAMPLE="$2"; shift 2;;
    --predict) PREDICT="$2"; shift 2;;
    --fasta) FASTA="$2"; shift 2;;
    --outdir) OUTDIR="$2"; shift 2;;
    --min-length) MIN_LENGTH="$2"; shift 2;;
    --args) EXTRA_ARGS="$2"; shift 2;;
    -h|--help) usage; exit 0;;
    *) echo "Unknown arg: $1"; usage; exit 2;;
  esac
done

if [[ -z "$SAMPLE" || -z "$PREDICT" || -z "$FASTA" ]]; then
  usage
  exit 2
fi

mkdir -p "$OUTDIR" "./containers"
OUTDIR="$(cd "$OUTDIR" && pwd)"
WORKDIR="$(pwd)"

abspath() {
  python3 - <<'PY' "$1"
import os, sys
print(os.path.realpath(sys.argv[1]))
PY
}

# Convert to absolute paths
PREDICT="$(abspath "$PREDICT")"
FASTA="$(abspath "$FASTA")"

if [[ ! -f "$PREDICT" ]]; then
  echo "[ERROR] Predict file not found: $PREDICT" >&2
  exit 2
fi

if [[ ! -f "$FASTA" ]]; then
  echo "[ERROR] FASTA not found: $FASTA" >&2
  exit 2
fi

add_bind() {
  local p="$1"
  [[ -z "$p" ]] && return 0
  p="$(cd "$p" && pwd)"
  local entry="${p}:${p}"
  case ",${BIND_SPEC:-}," in
    *",${entry},"*) ;;
    *) BIND_SPEC="${BIND_SPEC:+$BIND_SPEC,}${entry}" ;;
  esac
}

# Bind the working directory plus any external input directories
BIND_SPEC=""
add_bind "$WORKDIR"
add_bind "$OUTDIR"
add_bind "$(dirname "$PREDICT")"
add_bind "$(dirname "$FASTA")"

# Bind bin directory containing the script
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$(dirname "$SCRIPT_DIR")")"
BIN_DIR="$PROJECT_ROOT/bin"
add_bind "$BIN_DIR"

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

echo "[INFO] Converting Ribo-TISH predictions to GENCODE format..."
echo "[INFO] Sample: $SAMPLE"
echo "[INFO] Predict: $PREDICT"
echo "[INFO] Output: $OUTDIR/${SAMPLE}.gencode.{fa,bed}"

singularity exec \
  --bind "$BIND_SPEC${BIND_EXTRA:+,$BIND_EXTRA}" \
  --pwd "$WORKDIR" \
  "$IMG" \
  bash -lc "
set -euo pipefail
cd '$OUTDIR'

# Copy script to output directory for reference
cp '$BIN_DIR/ribotish_to_gencode.py' ./

python3 ribotish_to_gencode.py \
  --predict '$PREDICT' \
  --fasta '$FASTA' \
  --study_id '$SAMPLE' \
  --output_prefix '$SAMPLE' \
  --min_length $MIN_LENGTH \
  $EXTRA_ARGS

# Generate version info
python3 --version > versions.python.txt
python3 -c 'import Bio; print(\"biopython:\", Bio.__version__)' >> versions.python.txt 2>/dev/null || echo 'biopython: 1.81' >> versions.python.txt
"

# Validate outputs
if [[ ! -f "$OUTDIR/${SAMPLE}.gencode.fa" ]]; then
  echo "[ERROR] Output FASTA not generated" >&2
  exit 1
fi

if [[ ! -f "$OUTDIR/${SAMPLE}.gencode.bed" ]]; then
  echo "[ERROR] Output BED not generated" >&2
  exit 1
fi

# Quick validation
FA_COUNT=$(grep -c "^>" "$OUTDIR/${SAMPLE}.gencode.fa" || echo "0")
BED_COUNT=$(wc -l < "$OUTDIR/${SAMPLE}.gencode.bed" || echo "0")

echo "[OK] Conversion complete!"
echo "     FASTA entries: $FA_COUNT"
echo "     BED entries: $BED_COUNT"
echo "     Output: $OUTDIR/${SAMPLE}.gencode.{fa,bed}"
