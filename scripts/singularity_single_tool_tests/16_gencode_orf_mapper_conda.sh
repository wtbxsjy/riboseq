#!/usr/bin/env bash
set -euo pipefail

# Mirrors: modules/local/gencode_orf_mapper/main.nf
# Runs gencode-riboseqORFs ORF_mapper_to_GENCODE.py to unify ORF annotations
# Conda version - simpler and faster than Singularity

CONDA_ENV="gencode_orf_mapper"

usage() {
  cat <<'EOF'
Usage:
  16_gencode_orf_mapper_conda.sh --project ID --fasta merged.fa --bed merged.bed --ensembl-dir DIR [OPTIONS]

Required:
  --project       project/study ID (output prefix)
  --fasta         merged ORF sequences in GENCODE format (*.gencode.fa)
  --bed           merged ORF coordinates in GENCODE format (*.gencode.bed, 1-based)
  --ensembl-dir   directory containing Ensembl annotation files
                  (PROTEOME_FASTA, TRANSCRIPTOME_FASTA, SORTED_TRANSCRIPTOME_GTF, etc.)

Options:
  --outdir               output directory (default: ./out_gencode_orf_mapper)
  --min-length           minimum ORF length for mapping (default: 16)
  --collapse-threshold   similarity threshold for ORF merging (default: 0.9)
  --collapse-method      method for ORF deduplication: longest_string or psite_overlap (default: longest_string)
  --args                 extra args passed to ORF_mapper_to_GENCODE (quoted string)
  --conda-env            conda environment name (default: gencode_orf_mapper)

Output:
  ${PROJECT}.orfs.fa       - Unified ORF sequences
  ${PROJECT}.orfs.bed      - Unified ORF coordinates (1-based)
  ${PROJECT}.orfs.gtf      - GENCODE-format GTF annotation
  ${PROJECT}.orfs.out      - Detailed ORF features table
  ${PROJECT}.altmapped     - Alternative mappings
  ${PROJECT}.unmapped      - Unmapped ORFs

Examples:
  # Using conda environment
  16_gencode_orf_mapper_conda.sh --project PROJ1 --fasta merged.fa --bed merged.bed \
    --ensembl-dir /path/to/Ens110_GRCh38

Prerequisites:
  1. Create conda environment first:
     conda env create -f scripts/envs/gencode_orf_mapper.yml

  2. Or manually install dependencies:
     conda create -n gencode_orf_mapper python=3.9 biopython=1.81 pandas bedtools samtools -c bioconda -c conda-forge
EOF
}

PROJECT=""
FASTA=""
BED=""
ENSEMBL_DIR=""
OUTDIR="./out_gencode_orf_mapper"
MIN_LENGTH=16
COLLAPSE_THRESHOLD=0.9
COLLAPSE_METHOD="longest_string"
EXTRA_ARGS=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --project) PROJECT="$2"; shift 2;;
    --fasta) FASTA="$2"; shift 2;;
    --bed) BED="$2"; shift 2;;
    --ensembl-dir) ENSEMBL_DIR="$2"; shift 2;;
    --outdir) OUTDIR="$2"; shift 2;;
    --min-length) MIN_LENGTH="$2"; shift 2;;
    --collapse-threshold) COLLAPSE_THRESHOLD="$2"; shift 2;;
    --collapse-method) COLLAPSE_METHOD="$2"; shift 2;;
    --args) EXTRA_ARGS="$2"; shift 2;;
    --conda-env) CONDA_ENV="$2"; shift 2;;
    -h|--help) usage; exit 0;;
    *) echo "Unknown arg: $1"; usage; exit 2;;
  esac
done

if [[ -z "$PROJECT" || -z "$FASTA" || -z "$BED" || -z "$ENSEMBL_DIR" ]]; then
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
FASTA="$(abspath "$FASTA")"
BED="$(abspath "$BED")"
ENSEMBL_DIR="$(abspath "$ENSEMBL_DIR")"

if [[ ! -f "$FASTA" ]]; then
  echo "[ERROR] FASTA not found: $FASTA" >&2
  exit 2
fi

if [[ ! -f "$BED" ]]; then
  echo "[ERROR] BED not found: $BED" >&2
  exit 2
fi

if [[ ! -d "$ENSEMBL_DIR" ]]; then
  echo "[ERROR] Ensembl directory not found: $ENSEMBL_DIR" >&2
  exit 2
fi

# Validate Ensembl directory contents
REQUIRED_FILES=(
  "PROTEOME_FASTA"
  "TRANSCRIPTOME_FASTA"
  "SORTED_TRANSCRIPTOME_GTF"
  "TRANSCRIPT_SUPPORT"
  "PSITES_BED"
)

for file in "${REQUIRED_FILES[@]}"; do
  if [[ ! -f "$ENSEMBL_DIR/$file" ]]; then
    echo "[ERROR] Missing required Ensembl file: $ENSEMBL_DIR/$file" >&2
    echo "[HINT] Run prepare_ensembl_annotation to download Ensembl files" >&2
    exit 2
  fi
done

# Find gencode-riboseqORFs scripts
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$(dirname "$SCRIPT_DIR")")"
GENCODE_SCRIPTS_DIR="$PROJECT_ROOT/scripts/gencode-riboseqORFs"

if [[ ! -f "$GENCODE_SCRIPTS_DIR/ORF_mapper_to_GENCODE_v1.1.py" ]]; then
  echo "[ERROR] ORF mapper script not found: $GENCODE_SCRIPTS_DIR/ORF_mapper_to_GENCODE_v1.1.py" >&2
  exit 2
fi

# Check if conda environment exists
if ! conda env list | grep -q "^${CONDA_ENV} "; then
  echo "[ERROR] Conda environment '${CONDA_ENV}' not found" >&2
  echo "[HINT] Create it with: conda env create -f scripts/envs/gencode_orf_mapper.yml" >&2
  exit 2
fi

echo "[INFO] Running GENCODE ORF mapper..."
echo "[INFO] Project: $PROJECT"
echo "[INFO] Input FASTA: $FASTA"
echo "[INFO] Input BED: $BED"
echo "[INFO] Ensembl dir: $ENSEMBL_DIR"
echo "[INFO] Min length: $MIN_LENGTH aa"
echo "[INFO] Collapse threshold: $COLLAPSE_THRESHOLD"
echo "[INFO] Collapse method: $COLLAPSE_METHOD"
echo "[INFO] Conda env: $CONDA_ENV"
echo "[INFO] Output: $OUTDIR/${PROJECT}.orfs.*"

# Activate conda environment and run
eval "$(conda shell.bash hook)"
conda activate "$CONDA_ENV"

cd "$OUTDIR"

# Copy gencode-riboseqORFs scripts to output directory for reference
cp "$GENCODE_SCRIPTS_DIR/ORF_mapper_to_GENCODE_v1.1.py" ./
cp "$GENCODE_SCRIPTS_DIR/functions.py" ./

# Run the ORF mapper script from the local copy
python3 ORF_mapper_to_GENCODE_v1.1.py \
  -d "$ENSEMBL_DIR" \
  -f "$FASTA" \
  -b "$BED" \
  -o "$PROJECT" \
  -l $MIN_LENGTH \
  -c $COLLAPSE_THRESHOLD \
  -m $COLLAPSE_METHOD \
  $EXTRA_ARGS

# Generate version info
python3 --version > versions.python.txt
python3 -c 'import pandas; print("pandas:", pandas.__version__)' >> versions.python.txt 2>/dev/null || echo 'pandas: 1.3.5' >> versions.python.txt
python3 -c 'import Bio; print("biopython:", Bio.__version__)' >> versions.python.txt 2>/dev/null || echo 'biopython: 1.81' >> versions.python.txt
bedtools --version >> versions.python.txt 2>/dev/null || echo 'bedtools: 2.30.0' >> versions.python.txt

conda deactivate

cd "$WORKDIR"

# Validate outputs
EXPECTED_OUTPUTS=(
  "${PROJECT}.orfs.fa"
  "${PROJECT}.orfs.bed"
  "${PROJECT}.orfs.gtf"
  "${PROJECT}.orfs.out"
)

MISSING=0
for out in "${EXPECTED_OUTPUTS[@]}"; do
  if [[ ! -f "$OUTDIR/$out" ]]; then
    echo "[WARNING] Expected output not found: $out" >&2
    MISSING=$((MISSING + 1))
  fi
done

if [[ $MISSING -gt 0 ]]; then
  echo "[ERROR] $MISSING expected output files missing" >&2
  exit 1
fi

# Quick statistics
if [[ -f "$OUTDIR/${PROJECT}.orfs.out" ]]; then
  ORF_COUNT=$(tail -n +2 "$OUTDIR/${PROJECT}.orfs.out" | wc -l || echo "0")
  echo "[OK] GENCODE ORF mapping complete!"
  echo "     Unified ORFs: $ORF_COUNT"
  echo "     Output: $OUTDIR/${PROJECT}.orfs.*"
else
  echo "[OK] GENCODE ORF mapping complete!"
  echo "     Output: $OUTDIR/${PROJECT}.orfs.*"
fi
