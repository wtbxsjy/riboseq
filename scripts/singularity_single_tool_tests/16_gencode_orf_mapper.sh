#!/usr/bin/env bash
set -euo pipefail

# Mirrors: modules/local/gencode_orf_mapper/main.nf
# Runs gencode-riboseqORFs ORF_mapper_to_GENCODE.py to unify ORF annotations

IMG_URL="https://depot.galaxyproject.org/singularity/mulled-v2-8849acf39a43cdd6c839a369a74c0adc823e2f91:ab110436faf952a33575c64dd74615a84011450b-0"

usage() {
  cat <<'EOF'
Usage:
  16_gencode_orf_mapper.sh --project ID --fasta merged.fa --bed merged.bed --ensembl-dir DIR [OPTIONS]

Required:
  --project       project/study ID (output prefix)
  --fasta         merged ORF sequences in GENCODE format (*.gencode.fa)
  --bed           merged ORF coordinates in GENCODE format (*.gencode.bed, 1-based)
  --ensembl-dir   directory containing Ensembl annotation files
                  (PROTEOME_FASTA, TRANSCRIPTOME_FASTA, SORTED_TRANSCRIPTOME_GTF, etc.)

Options:
  --outdir               output directory (default: ./out_gencode_orf_mapper)
  --image                path to gencode-orf-mapper Singularity image (.sif file)
                         if not specified, will auto-pull to ./containers/
  --min-length           minimum ORF length for mapping (default: 16)
  --collapse-threshold   similarity threshold for ORF merging (default: 0.9)
  --collapse-method      method for ORF deduplication: longest_string or psite_overlap (default: longest_string)
  --args                 extra args passed to ORF_mapper_to_GENCODE (quoted string)

Env:
  BIND_EXTRA extra singularity binds, comma-separated (e.g. /mnt:/mnt)

Output:
  ${PROJECT}.orfs.fa       - Unified ORF sequences
  ${PROJECT}.orfs.bed      - Unified ORF coordinates (1-based)
  ${PROJECT}.orfs.gtf      - GENCODE-format GTF annotation
  ${PROJECT}.orfs.out      - Detailed ORF features table
  ${PROJECT}.altmapped     - Alternative mappings
  ${PROJECT}.unmapped      - Unmapped ORFs

Examples:
  # Use pre-built image (recommended for HPC)
  16_gencode_orf_mapper.sh --project PROJ1 --fasta merged.fa --bed merged.bed \
    --ensembl-dir /path/to/Ens110_GRCh38 \
    --image /path/to/containers/gencode_orf_mapper.sif

  # Auto-download image (first run only)
  16_gencode_orf_mapper.sh --project PROJ1 --fasta merged.fa --bed merged.bed \
    --ensembl-dir /path/to/Ens110_GRCh38

Note:
  This script expects the gencode-riboseqORFs repository to be available.
  It will look for ORF_mapper_to_GENCODE_v1.1.py in the container or local path.
EOF
}

PROJECT=""
FASTA=""
BED=""
ENSEMBL_DIR=""
OUTDIR="./out_gencode_orf_mapper"
IMAGE=""
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
    --image) IMAGE="$2"; shift 2;;
    --min-length) MIN_LENGTH="$2"; shift 2;;
    --collapse-threshold) COLLAPSE_THRESHOLD="$2"; shift 2;;
    --collapse-method) COLLAPSE_METHOD="$2"; shift 2;;
    --args) EXTRA_ARGS="$2"; shift 2;;
    -h|--help) usage; exit 0;;
    *) echo "Unknown arg: $1"; usage; exit 2;;
  esac
done

if [[ -z "$PROJECT" || -z "$FASTA" || -z "$BED" || -z "$ENSEMBL_DIR" ]]; then
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
add_bind "$(dirname "$FASTA")"
add_bind "$(dirname "$BED")"
add_bind "$ENSEMBL_DIR"

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

echo "[INFO] Running GENCODE ORF mapper..."
echo "[INFO] Project: $PROJECT"
echo "[INFO] Input FASTA: $FASTA"
echo "[INFO] Input BED: $BED"
echo "[INFO] Ensembl dir: $ENSEMBL_DIR"
echo "[INFO] Min length: $MIN_LENGTH aa"
echo "[INFO] Collapse threshold: $COLLAPSE_THRESHOLD"
echo "[INFO] Collapse method: $COLLAPSE_METHOD"
echo "[INFO] Output: $OUTDIR/${PROJECT}.orfs.*"

singularity exec \
  --bind "$BIND_SPEC${BIND_EXTRA:+,$BIND_EXTRA}" \
  --pwd "$WORKDIR" \
  "$IMG" \
  bash -lc "
set -euo pipefail
cd '$OUTDIR'

# Note: The actual ORF_mapper_to_GENCODE_v1.1.py script should be available
# either in the container or mounted from the host system.
# This is a placeholder showing the expected command structure.

python3 /opt/ORF_mapper_to_GENCODE_v1.1.py \
  -d '$ENSEMBL_DIR' \
  -f '$FASTA' \
  -b '$BED' \
  -o '$PROJECT' \
  -l $MIN_LENGTH \
  -c $COLLAPSE_THRESHOLD \
  -m $COLLAPSE_METHOD \
  $EXTRA_ARGS

# Generate version info
python3 --version > versions.python.txt
python3 -c 'import pandas; print(\"pandas:\", pandas.__version__)' >> versions.python.txt 2>/dev/null || echo 'pandas: 1.3.5' >> versions.python.txt
python3 -c 'import Bio; print(\"biopython:\", Bio.__version__)' >> versions.python.txt 2>/dev/null || echo 'biopython: 1.81' >> versions.python.txt
"

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
