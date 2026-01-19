#!/usr/bin/env bash
set -euo pipefail

################################################################################
# Script: diagnose_gencode_image.sh
# Purpose: Diagnose and fix GENCODE ORF mapper image issues
################################################################################

usage() {
  cat <<'EOF'
Usage:
  diagnose_gencode_image.sh [--image PATH]

Options:
  --image PATH    Path to the gencode_orf_mapper_mulled.sif image to diagnose

Description:
  This script diagnoses issues with the GENCODE ORF mapper Singularity image
  and suggests solutions.

EOF
}

IMAGE=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --image) IMAGE="$2"; shift 2;;
    -h|--help) usage; exit 0;;
    *) echo "Unknown arg: $1"; usage; exit 2;;
  esac
done

if [[ -z "$IMAGE" ]]; then
  echo "[ERROR] Please specify --image path"
  usage
  exit 2
fi

if [[ ! -f "$IMAGE" ]]; then
  echo "[ERROR] Image not found: $IMAGE"
  exit 2
fi

echo "========================================================================"
echo "  Diagnosing GENCODE ORF Mapper Image"
echo "========================================================================"
echo ""
echo "Image: $IMAGE"
echo "Size: $(du -h "$IMAGE" | cut -f1)"
echo ""

echo "------------------------------------------------------------------------"
echo "Test 1: Check if image is valid"
echo "------------------------------------------------------------------------"
if ! singularity inspect "$IMAGE" &>/dev/null; then
  echo "[ERROR] Image is corrupted or invalid"
  echo "[SOLUTION] Re-download the image"
  exit 1
fi
echo "[OK] Image is valid"
echo ""

echo "------------------------------------------------------------------------"
echo "Test 2: List available executables in PATH"
echo "------------------------------------------------------------------------"
echo "Searching for python/python3..."
singularity exec "$IMAGE" sh -c 'find /usr /opt /bin 2>/dev/null | grep -E "bin/(python|python3)$" | head -20' || echo "No python found in standard locations"
echo ""

echo "------------------------------------------------------------------------"
echo "Test 3: Check available commands"
echo "------------------------------------------------------------------------"
echo "Checking for python alternatives..."
singularity exec "$IMAGE" sh -c 'command -v python || command -v python3 || command -v python2 || echo "No python found"'
echo ""

echo "------------------------------------------------------------------------"
echo "Test 4: List all files in /usr/bin"
echo "------------------------------------------------------------------------"
singularity exec "$IMAGE" ls -la /usr/bin/ 2>/dev/null | grep python || echo "No python in /usr/bin"
echo ""

echo "------------------------------------------------------------------------"
echo "Test 5: Check image metadata"
echo "------------------------------------------------------------------------"
singularity inspect --deffile "$IMAGE" 2>/dev/null || echo "No definition file available"
echo ""

echo "------------------------------------------------------------------------"
echo "Test 6: Interactive shell test"
echo "------------------------------------------------------------------------"
echo "Starting interactive shell... (type 'exit' to quit)"
echo "Try commands: which python3, find / -name python3 2>/dev/null | head"
echo ""
singularity shell "$IMAGE"
