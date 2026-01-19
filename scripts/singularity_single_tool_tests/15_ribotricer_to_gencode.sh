#!/usr/bin/env bash
set -euo pipefail

# Mirrors: modules/local/convert_ribotricer_to_gencode/main.nf
# Converts Ribotricer translating_ORFs.tsv to gencode-riboseqORFs compatible format

IMG_URL="https://depot.galaxyproject.org/singularity/biopython:1.81"

usage() {
  cat <<'EOF'
Usage:
  15_ribotricer_to_gencode.sh --sample ID --tsv translating_ORFs.tsv --fasta genome.fa [OPTIONS]

Required:
  --sample   sample ID (used as study_id in output)
  --tsv      Ribotricer translating_ORFs.tsv output file
  --fasta    genome FASTA file

Options:
  --outdir           output directory (default: ./out_ribotricer_to_gencode)
  --image            path to biopython Singularity image (.sif file)
                     if not specified, will auto-pull to ./containers/biopython_1.81.sif
  --min-length       minimum ORF length in amino acids (default: 16)
  --min-phase-score  minimum phase score for filtering (default: 0.5)
  --args             extra args passed to ribotricer_to_gencode.py (quoted string)

Env:
  BIND_EXTRA extra singularity binds, comma-separated (e.g. /mnt:/mnt)

Output:
  ${SAMPLE}.gencode.fa   - FASTA in gencode-riboseqORFs format
  ${SAMPLE}.gencode.bed  - BED (1-based) in gencode-riboseqORFs format

Examples:
  # Use pre-built image (recommended for HPC)
  15_ribotricer_to_gencode.sh --sample S1 --tsv orfs.tsv --fasta genome.fa \
    --image /path/to/containers/biopython_1.81.sif

  # Auto-download image (first run only)
  15_ribotricer_to_gencode.sh --sample S1 --tsv orfs.tsv --fasta genome.fa
EOF
}

SAMPLE=""
TSV=""
FASTA=""
OUTDIR="./out_ribotricer_to_gencode"
IMAGE=""
MIN_LENGTH=16
MIN_PHASE_SCORE=0.5
EXTRA_ARGS=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --sample) SAMPLE="$2"; shift 2;;
    --tsv) TSV="$2"; shift 2;;
    --fasta) FASTA="$2"; shift 2;;
    --outdir) OUTDIR="$2"; shift 2;;
    --image) IMAGE="$2"; shift 2;;
    --min-length) MIN_LENGTH="$2"; shift 2;;
    --min-phase-score) MIN_PHASE_SCORE="$2"; shift 2;;
    --args) EXTRA_ARGS="$2"; shift 2;;
    -h|--help) usage; exit 0;;
    *) echo "Unknown arg: $1"; usage; exit 2;;
  esac
done

if [[ -z "$SAMPLE" || -z "$TSV" || -z "$FASTA" ]]; then
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
TSV="$(abspath "$TSV")"
FASTA="$(abspath "$FASTA")"

if [[ ! -f "$TSV" ]]; then
  echo "[ERROR] TSV file not found: $TSV" >&2
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
add_bind "$(dirname "$TSV")"
add_bind "$(dirname "$FASTA")"

# Bind scripts directory containing the converter script
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$(dirname "$SCRIPT_DIR")")"
CONVERTER_DIR="$PROJECT_ROOT/scripts/gencode_converters"
add_bind "$CONVERTER_DIR"

pull_img() {
  local url="$1"
  local base
  base="$(basename "$url")"
  base="${base//:/_}"
  local sif="$(pwd)/containers/${base}.sif"
  if [[ ! -f "$sif" ]]; then
    echo "[INFO] Image not found at $sif, downloading..."
    mkdir -p "$(dirname "$sif")"
    singularity pull --disable-cache --force "$sif" "$url"
  fi
  echo "$sif"
}

# Determine which image to use
if [[ -n "$IMAGE" ]]; then
  # User specified image path
  if [[ ! -f "$IMAGE" ]]; then
    echo "[ERROR] Specified image not found: $IMAGE" >&2
    exit 2
  fi
  IMG="$(abspath "$IMAGE")"
  echo "[INFO] Using user-specified image: $IMG"
else
  # Auto-pull image
  echo "[INFO] No --image specified, will auto-pull if needed"
  IMG="$(pull_img "$IMG_URL")"
  echo "[INFO] Using image: $IMG"
fi

# Bind the image directory
add_bind "$(dirname "$IMG")"

echo "[INFO] Converting Ribotricer predictions to GENCODE format..."
echo "[INFO] Sample: $SAMPLE"
echo "[INFO] TSV: $TSV"
echo "[INFO] Min length: $MIN_LENGTH aa"
echo "[INFO] Min phase score: $MIN_PHASE_SCORE"
echo "[INFO] Output: $OUTDIR/${SAMPLE}.gencode.{fa,bed}"

singularity exec \
  --bind "$BIND_SPEC${BIND_EXTRA:+,$BIND_EXTRA}" \
  --pwd "$WORKDIR" \
  "$IMG" \
  bash -lc "
set -euo pipefail
cd '$OUTDIR'

# Copy script to output directory for reference
cp '$CONVERTER_DIR/ribotricer_to_gencode.py' ./

python3 ribotricer_to_gencode.py \
  --tsv '$TSV' \
  --fasta '$FASTA' \
  --study_id '$SAMPLE' \
  --output_prefix '$SAMPLE' \
  --min_length $MIN_LENGTH \
  --min_phase_score $MIN_PHASE_SCORE \
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
