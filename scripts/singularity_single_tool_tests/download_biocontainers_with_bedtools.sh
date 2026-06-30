#!/usr/bin/env bash
set -euo pipefail

################################################################################
# Script: download_biocontainers_with_bedtools.sh
# Purpose: Download biocontainers image with Python, Biopython, and bedtools
################################################################################

usage() {
  cat <<'EOF'
Usage:
  download_biocontainers_with_bedtools.sh [--outdir DIR]

Options:
  --outdir DIR    Output directory (default: ./containers)

Description:
  Downloads a biocontainers image that includes:
  - Python 3
  - Biopython
  - bedtools
  - Other bioinformatics tools

EOF
}

OUTDIR="./containers"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --outdir) OUTDIR="$2"; shift 2;;
    -h|--help) usage; exit 0;;
    *) echo "Unknown arg: $1"; usage; exit 2;;
  esac
done

mkdir -p "$OUTDIR"
OUTDIR="$(cd "$OUTDIR" && pwd)"

echo "========================================================================"
echo "  Downloading biocontainers image with bedtools"
echo "========================================================================"
echo ""

# Try multiple image sources
IMAGES=(
  "docker://quay.io/biocontainers/bedtools:2.31.1--hf5e1c6e_1"
  "docker://quay.io/biocontainers/mulled-v2-8186960447c5cb2faa697666dc1e6d919ad23f3e:3127fcae6b6bdaf8181e21a26ae61231030a9fcb-0"
)

OUTPUT_NAME="biotools_bedtools_python.sif"

for img in "${IMAGES[@]}"; do
  echo "------------------------------------------------------------------------"
  echo "Trying: $img"
  echo "Output: $OUTDIR/$OUTPUT_NAME"
  echo "------------------------------------------------------------------------"

  if singularity pull --force "$OUTDIR/$OUTPUT_NAME" "$img"; then
    echo ""
    echo "[OK] Successfully downloaded image!"
    echo "Size: $(du -h "$OUTDIR/$OUTPUT_NAME" | cut -f1)"

    # Verify contents
    echo ""
    echo "Verifying image contents..."
    echo "=== Python ==="
    singularity exec "$OUTDIR/$OUTPUT_NAME" python3 --version 2>&1 || echo "Python not found, trying python..."
    singularity exec "$OUTDIR/$OUTPUT_NAME" python --version 2>&1 || echo "No python found"

    echo ""
    echo "=== Bedtools ==="
    singularity exec "$OUTDIR/$OUTPUT_NAME" bedtools --version 2>&1 || \
    singularity exec "$OUTDIR/$OUTPUT_NAME" intersectBed 2>&1 | head -3 || \
    echo "Bedtools not found in expected location"

    echo ""
    echo "=== Biopython ==="
    singularity exec "$OUTDIR/$OUTPUT_NAME" python3 -c "import Bio; print('Biopython:', Bio.__version__)" 2>&1 || \
    singularity exec "$OUTDIR/$OUTPUT_NAME" python -c "import Bio; print('Biopython:', Bio.__version__)" 2>&1 || \
    echo "Biopython not found"

    echo ""
    echo "========================================================================"
    echo "  ✓ Download complete!"
    echo "========================================================================"
    echo ""
    echo "Use with:"
    echo "  bash 16_gencode_orf_mapper.sh --image $OUTDIR/$OUTPUT_NAME ..."
    echo ""
    exit 0
  else
    echo "[WARNING] Failed to download from this source"
    echo ""
  fi
done

echo ""
echo "[ERROR] All image sources failed"
echo ""
echo "Alternative: Install bedtools on the host system:"
echo "  sudo apt-get install bedtools"
echo "  # or"
echo "  conda install -c bioconda bedtools"
echo ""
exit 1
