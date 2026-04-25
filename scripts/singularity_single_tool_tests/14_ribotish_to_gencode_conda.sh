#!/usr/bin/env bash
set -euo pipefail

# Legacy standalone converter; no longer mirrored by an active Nextflow module.
# Converts Ribo-TISH predict output to gencode-riboseqORFs compatible format
# Conda version - simpler and faster than Singularity

CONDA_ENV="ribotish_to_gencode"

usage() {
  cat <<'EOF'
Usage:
  14_ribotish_to_gencode_conda.sh --sample ID --predict ribotish_pred.txt --fasta genome.fa [OPTIONS]

Required:
  --sample   sample ID (used as study_id in output)
  --predict  Ribo-TISH predict output file (*_pred.txt)
  --fasta    genome FASTA file

Options:
  --outdir      output directory (default: ./out_ribotish_to_gencode)
  --min-length  minimum ORF length in amino acids (default: 16)
  --args        extra args passed to ribotish_to_gencode.py (quoted string)
  --conda-env   conda environment name (default: ribotish_to_gencode)

Output:
  ${SAMPLE}.gencode.fa   - FASTA in gencode-riboseqORFs format
  ${SAMPLE}.gencode.bed  - BED (1-based) in gencode-riboseqORFs format

Examples:
  # Using conda environment
  14_ribotish_to_gencode_conda.sh --sample S1 --predict pred.txt --fasta genome.fa

  # With custom conda env name
  14_ribotish_to_gencode_conda.sh --sample S1 --predict pred.txt --fasta genome.fa \
    --conda-env my_custom_env

Prerequisites:
  1. Create conda environment first:
     conda env create -f scripts/envs/ribotish_to_gencode.yml

  2. Or manually install dependencies:
     conda create -n ribotish_to_gencode python=3.9 biopython=1.81 pyfaidx bedtools samtools -c bioconda -c conda-forge
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
    --conda-env) CONDA_ENV="$2"; shift 2;;
    -h|--help) usage; exit 0;;
    *) echo "Unknown arg: $1"; usage; exit 2;;
  esac
done

if [[ -z "$SAMPLE" || -z "$PREDICT" || -z "$FASTA" ]]; then
  usage
  exit 2
fi

mkdir -p "$OUTDIR"
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

# Find converter script
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$(dirname "$SCRIPT_DIR")")"
CONVERTER_SCRIPT="$PROJECT_ROOT/scripts/gencode_converters/ribotish_to_gencode.py"

if [[ ! -f "$CONVERTER_SCRIPT" ]]; then
  echo "[ERROR] Converter script not found: $CONVERTER_SCRIPT" >&2
  exit 2
fi

# Check if conda environment exists
if ! conda env list | grep -q "^${CONDA_ENV} "; then
  echo "[ERROR] Conda environment '${CONDA_ENV}' not found" >&2
  echo "[HINT] Create it with: conda env create -f scripts/envs/ribotish_to_gencode.yml" >&2
  exit 2
fi

echo "[INFO] Converting Ribo-TISH predictions to GENCODE format..."
echo "[INFO] Sample: $SAMPLE"
echo "[INFO] Predict: $PREDICT"
echo "[INFO] FASTA: $FASTA"
echo "[INFO] Conda env: $CONDA_ENV"
echo "[INFO] Output: $OUTDIR/${SAMPLE}.gencode.{fa,bed}"

# Activate conda environment and run
eval "$(conda shell.bash hook)"
conda activate "$CONDA_ENV"

cd "$OUTDIR"

# Copy script to output directory for reference
cp "$CONVERTER_SCRIPT" ./

# Run the converter
python3 ribotish_to_gencode.py \
  --predict "$PREDICT" \
  --fasta "$FASTA" \
  --study_id "$SAMPLE" \
  --output_prefix "$SAMPLE" \
  --min_length $MIN_LENGTH \
  $EXTRA_ARGS

# Generate version info
python3 --version > versions.python.txt
python3 -c 'import Bio; print("biopython:", Bio.__version__)' >> versions.python.txt 2>/dev/null || echo 'biopython: 1.81' >> versions.python.txt
python3 -c 'import pyfaidx; print("pyfaidx:", pyfaidx.__version__)' >> versions.python.txt 2>/dev/null || echo 'pyfaidx: not installed' >> versions.python.txt
bedtools --version >> versions.python.txt 2>/dev/null || echo 'bedtools: not available' >> versions.python.txt

conda deactivate

cd "$WORKDIR"

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
