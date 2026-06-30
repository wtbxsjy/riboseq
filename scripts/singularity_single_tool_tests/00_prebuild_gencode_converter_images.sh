#!/usr/bin/env bash
set -euo pipefail

################################################################################
# Script: 00_prebuild_gencode_converter_images.sh
# Purpose: Pre-build Singularity images for gencode converter scripts (14 & 15)
#
# This script downloads the required Singularity image for:
#   - 14_ribotish_to_gencode.sh
#   - 15_ribotricer_to_gencode.sh
#
# Both scripts use the same biopython image.
################################################################################

usage() {
  cat <<'EOF'
Usage:
  00_prebuild_gencode_converter_images.sh [--outdir DIR] [--cache-dir DIR]

Options:
  --outdir DIR       Output directory for .sif files (default: ./containers)
  --cache-dir DIR    Singularity cache directory (default: system default)
  --force            Force rebuild even if image exists
  -h, --help         Show this help message

Description:
  Downloads and builds the Singularity image required for gencode converter
  scripts (14_ribotish_to_gencode.sh and 15_ribotricer_to_gencode.sh).

  After building, you can transfer the entire containers/ directory to your
  HPC cluster and use it directly with the converter scripts.

Examples:
  # Build to default location (./containers/)
  bash 00_prebuild_gencode_converter_images.sh

  # Build to custom location
  bash 00_prebuild_gencode_converter_images.sh --outdir /path/to/containers

  # Force rebuild
  bash 00_prebuild_gencode_converter_images.sh --force

  # Use custom cache directory (for systems with limited /tmp space)
  bash 00_prebuild_gencode_converter_images.sh --cache-dir /scratch/singularity_cache

Transfer to server:
  # On local machine after building:
  rsync -avP containers/ user@server:/path/to/riboseq/containers/

  # Or use scp:
  scp containers/biopython_1.81.sif user@server:/path/to/riboseq/containers/

EOF
}

OUTDIR="./containers"
CACHE_DIR=""
FORCE=false

while [[ $# -gt 0 ]]; do
  case "$1" in
    --outdir) OUTDIR="$2"; shift 2;;
    --cache-dir) CACHE_DIR="$2"; shift 2;;
    --force) FORCE=true; shift;;
    -h|--help) usage; exit 0;;
    *) echo "[ERROR] Unknown argument: $1"; usage; exit 2;;
  esac
done

# Create output directory
mkdir -p "$OUTDIR"
OUTDIR="$(cd "$OUTDIR" && pwd)"

# Set cache directory if specified
if [[ -n "$CACHE_DIR" ]]; then
  mkdir -p "$CACHE_DIR"
  export SINGULARITY_CACHEDIR="$(cd "$CACHE_DIR" && pwd)"
  echo "[INFO] Using custom Singularity cache directory: $SINGULARITY_CACHEDIR"
fi

echo "========================================================================"
echo "  Pre-building Singularity images for GENCODE converters"
echo "========================================================================"
echo ""
echo "Output directory: $OUTDIR"
echo ""

# Image 1: Biopython for scripts 14 and 15
IMG_URL_BIOPYTHON="https://depot.galaxyproject.org/singularity/biopython:1.81"
IMG_NAME_BIOPYTHON="biopython_1.81.sif"
IMG_PATH_BIOPYTHON="$OUTDIR/$IMG_NAME_BIOPYTHON"

# Image 2: Mulled container for script 16 (includes biopython, pandas, bedtools, gffread)
IMG_URL_GENCODE="https://depot.galaxyproject.org/singularity/mulled-v2-8849acf39a43cdd6c839a369a74c0adc823e2f91:ab110436faf952a33575c64dd74615a84011450b-0"
IMG_NAME_GENCODE="gencode_orf_mapper_mulled.sif"
IMG_PATH_GENCODE="$OUTDIR/$IMG_NAME_GENCODE"

# Build Biopython image
echo "------------------------------------------------------------------------"
echo "Image 1: Biopython 1.81"
echo "------------------------------------------------------------------------"
echo "URL: $IMG_URL_BIOPYTHON"
echo "Output: $IMG_PATH_BIOPYTHON"
echo "Used by: 14_ribotish_to_gencode.sh, 15_ribotricer_to_gencode.sh"
echo ""

if [[ -f "$IMG_PATH_BIOPYTHON" ]] && [[ "$FORCE" == "false" ]]; then
  echo "[INFO] Image already exists: $IMG_PATH_BIOPYTHON"
  echo "[INFO] File size: $(du -h "$IMG_PATH_BIOPYTHON" | cut -f1)"
  echo "[INFO] Use --force to rebuild"
  echo ""
else
  echo "[INFO] Pulling image from Galaxy Depot..."
  echo "[INFO] This may take several minutes depending on network speed..."
  echo ""

  singularity pull --force --disable-cache "$IMG_PATH_BIOPYTHON" "$IMG_URL_BIOPYTHON"

  if [[ $? -eq 0 ]]; then
    echo ""
    echo "[OK] Image built successfully!"
    echo "     Path: $IMG_PATH_BIOPYTHON"
    echo "     Size: $(du -h "$IMG_PATH_BIOPYTHON" | cut -f1)"
    echo ""
  else
    echo ""
    echo "[ERROR] Failed to build image from $IMG_URL_BIOPYTHON"
    exit 1
  fi
fi

# Build GENCODE ORF mapper image
echo "------------------------------------------------------------------------"
echo "Image 2: GENCODE ORF Mapper (mulled container)"
echo "------------------------------------------------------------------------"
echo "URL: $IMG_URL_GENCODE"
echo "Output: $IMG_PATH_GENCODE"
echo "Used by: 16_gencode_orf_mapper.sh"
echo "Contains: biopython, pandas, bedtools, gffread"
echo ""

if [[ -f "$IMG_PATH_GENCODE" ]] && [[ "$FORCE" == "false" ]]; then
  echo "[INFO] Image already exists: $IMG_PATH_GENCODE"
  echo "[INFO] File size: $(du -h "$IMG_PATH_GENCODE" | cut -f1)"
  echo "[INFO] Use --force to rebuild"
  echo ""
else
  echo "[INFO] Pulling image from Galaxy Depot..."
  echo "[INFO] This may take several minutes depending on network speed..."
  echo ""

  singularity pull --force --disable-cache "$IMG_PATH_GENCODE" "$IMG_URL_GENCODE"

  if [[ $? -eq 0 ]]; then
    echo ""
    echo "[OK] Image built successfully!"
    echo "     Path: $IMG_PATH_GENCODE"
    echo "     Size: $(du -h "$IMG_PATH_GENCODE" | cut -f1)"
    echo ""
  else
    echo ""
    echo "[ERROR] Failed to build image from $IMG_URL_GENCODE"
    exit 1
  fi
fi

# Verify image
echo "------------------------------------------------------------------------"
echo "Verifying images..."
echo "------------------------------------------------------------------------"

echo ""
echo "=== Image 1: Biopython ==="
if ! singularity exec "$IMG_PATH_BIOPYTHON" python3 --version; then
  echo "[ERROR] Image verification failed: cannot run python3"
  exit 1
fi

if ! singularity exec "$IMG_PATH_BIOPYTHON" python3 -c "import Bio; print('Biopython version:', Bio.__version__)"; then
  echo "[ERROR] Image verification failed: cannot import Biopython"
  exit 1
fi

echo ""
echo "=== Image 2: GENCODE ORF Mapper ==="
if ! singularity exec "$IMG_PATH_GENCODE" python3 --version; then
  echo "[ERROR] Image verification failed: cannot run python3"
  exit 1
fi

if ! singularity exec "$IMG_PATH_GENCODE" python3 -c "import Bio; print('Biopython version:', Bio.__version__)"; then
  echo "[ERROR] Image verification failed: cannot import Biopython"
  exit 1
fi

if ! singularity exec "$IMG_PATH_GENCODE" python3 -c "import pandas; print('pandas version:', pandas.__version__)"; then
  echo "[ERROR] Image verification failed: cannot import pandas"
  exit 1
fi

if ! singularity exec "$IMG_PATH_GENCODE" bedtools --version; then
  echo "[ERROR] Image verification failed: cannot run bedtools"
  exit 1
fi

echo ""
echo "========================================================================"
echo "  ✓ Pre-build complete!"
echo "========================================================================"
echo ""
echo "Built images:"
ls -lh "$OUTDIR"/*.sif
echo ""
echo "------------------------------------------------------------------------"
echo "Next steps:"
echo "------------------------------------------------------------------------"
echo ""
echo "1. Transfer to your HPC cluster:"
echo "   rsync -avP $OUTDIR/ user@server:/path/to/riboseq/containers/"
echo ""
echo "2. Or compress and transfer:"
echo "   tar -czf gencode_converter_images.tar.gz -C \"$OUTDIR\" ."
echo "   scp gencode_converter_images.tar.gz user@server:/path/to/"
echo "   # On server:"
echo "   tar -xzf gencode_converter_images.tar.gz -C /path/to/riboseq/containers/"
echo ""
echo "3. Verify on server:"
echo "   singularity exec biopython_1.81.sif python3 --version"
echo "   singularity exec gencode_orf_mapper_mulled.sif python3 -c 'import pandas'"
echo ""
echo "4. Run the converter scripts:"
echo "   # Scripts 14 and 15 use biopython_1.81.sif"
echo "   bash 14_ribotish_to_gencode.sh --sample test --predict data.txt --fasta genome.fa \\"
echo "     --image /path/to/biopython_1.81.sif"
echo ""
echo "   bash 15_ribotricer_to_gencode.sh --sample test --tsv orfs.tsv --fasta genome.fa \\"
echo "     --image /path/to/biopython_1.81.sif"
echo ""
echo "   # Script 16 uses gencode_orf_mapper_mulled.sif"
echo "   bash 16_gencode_orf_mapper.sh --project test --fasta merged.fa --bed merged.bed \\"
echo "     --ensembl-dir /path/to/ensembl --image /path/to/gencode_orf_mapper_mulled.sif"
echo ""
echo "========================================================================"
