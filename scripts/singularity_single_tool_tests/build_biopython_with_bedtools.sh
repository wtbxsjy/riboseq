#!/usr/bin/env bash
set -euo pipefail

################################################################################
# Script: build_biopython_with_bedtools.sh
# Purpose: Build Singularity image with Python, Biopython, and bedtools
################################################################################

usage() {
  cat <<'EOF'
Usage:
  build_biopython_with_bedtools.sh [OPTIONS]

Options:
  --outdir DIR       Output directory for .sif file (default: ./containers)
  --name NAME        Output image name (default: biopython_1.81.sif)
  --use-sudo         Use sudo for building (required on most systems)
  --use-fakeroot     Use --fakeroot instead of sudo (if available)
  -h, --help         Show this help message

Description:
  Builds a Singularity image containing:
  - Python 3.12
  - Biopython 1.81
  - bedtools 2.30
  - pandas, numpy

  This image can be used for scripts 14, 15, and 16.

Examples:
  # Build with sudo (most common)
  bash build_biopython_with_bedtools.sh --use-sudo

  # Build with fakeroot (if available on your system)
  bash build_biopython_with_bedtools.sh --use-fakeroot

  # Build to custom location
  bash build_biopython_with_bedtools.sh --use-sudo --outdir /data/containers

EOF
}

OUTDIR="./containers"
IMAGE_NAME="biopython_1.81.sif"
USE_SUDO=false
USE_FAKEROOT=false

while [[ $# -gt 0 ]]; do
  case "$1" in
    --outdir) OUTDIR="$2"; shift 2;;
    --name) IMAGE_NAME="$2"; shift 2;;
    --use-sudo) USE_SUDO=true; shift;;
    --use-fakeroot) USE_FAKEROOT=true; shift;;
    -h|--help) usage; exit 0;;
    *) echo "[ERROR] Unknown argument: $1"; usage; exit 2;;
  esac
done

if [[ "$USE_SUDO" == "true" ]] && [[ "$USE_FAKEROOT" == "true" ]]; then
  echo "[ERROR] Cannot use both --use-sudo and --use-fakeroot"
  exit 2
fi

if [[ "$USE_SUDO" == "false" ]] && [[ "$USE_FAKEROOT" == "false" ]]; then
  echo "[ERROR] Must specify either --use-sudo or --use-fakeroot"
  echo ""
  usage
  exit 2
fi

# Create output directory
mkdir -p "$OUTDIR"
OUTDIR="$(cd "$OUTDIR" && pwd)"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DEF_FILE="$SCRIPT_DIR/Singularity.biopython_with_bedtools.def"

if [[ ! -f "$DEF_FILE" ]]; then
  echo "[ERROR] Definition file not found: $DEF_FILE"
  exit 2
fi

OUTPUT_PATH="$OUTDIR/$IMAGE_NAME"

echo "========================================================================"
echo "  Building Biopython + bedtools Singularity Image"
echo "========================================================================"
echo ""
echo "Definition file: $DEF_FILE"
echo "Output: $OUTPUT_PATH"
echo ""

# Build command
BUILD_CMD="singularity build"
if [[ "$USE_SUDO" == "true" ]]; then
  BUILD_CMD="sudo $BUILD_CMD"
  echo "Build method: sudo"
elif [[ "$USE_FAKEROOT" == "true" ]]; then
  BUILD_CMD="$BUILD_CMD --fakeroot"
  echo "Build method: fakeroot"
fi
echo ""

echo "------------------------------------------------------------------------"
echo "Starting build..."
echo "------------------------------------------------------------------------"
echo "This may take 5-10 minutes depending on network speed and system resources"
echo ""

# Execute build
if $BUILD_CMD --force "$OUTPUT_PATH" "$DEF_FILE"; then
  echo ""
  echo "========================================================================"
  echo "  ✓ Build successful!"
  echo "========================================================================"
  echo ""
  echo "Image: $OUTPUT_PATH"
  echo "Size: $(du -h "$OUTPUT_PATH" | cut -f1)"
  echo ""

  # Verify the image
  echo "------------------------------------------------------------------------"
  echo "Verifying image contents..."
  echo "------------------------------------------------------------------------"
  echo ""

  echo "=== Python ==="
  singularity exec "$OUTPUT_PATH" python3 --version

  echo ""
  echo "=== Biopython ==="
  singularity exec "$OUTPUT_PATH" python3 -c "import Bio; print('Biopython version:', Bio.__version__)"

  echo ""
  echo "=== bedtools ==="
  singularity exec "$OUTPUT_PATH" bedtools --version

  echo ""
  echo "=== pandas (optional) ==="
  singularity exec "$OUTPUT_PATH" python3 -c "import pandas; print('pandas version:', pandas.__version__)" || echo "pandas not found (non-critical)"

  echo ""
  echo "========================================================================"
  echo "  ✓ All dependencies verified!"
  echo "========================================================================"
  echo ""
  echo "You can now use this image with scripts 14, 15, and 16:"
  echo ""
  echo "  bash 14_ribotish_to_gencode.sh --image $OUTPUT_PATH ..."
  echo "  bash 15_ribotricer_to_gencode.sh --image $OUTPUT_PATH ..."
  echo "  bash 16_gencode_orf_mapper.sh --image $OUTPUT_PATH ..."
  echo ""
else
  echo ""
  echo "========================================================================"
  echo "  ✗ Build failed!"
  echo "========================================================================"
  echo ""
  echo "Common issues:"
  echo "  1. Insufficient permissions - try with --use-sudo"
  echo "  2. Network connectivity - check internet connection"
  echo "  3. Disk space - ensure enough space in $OUTDIR"
  echo ""
  exit 1
fi
