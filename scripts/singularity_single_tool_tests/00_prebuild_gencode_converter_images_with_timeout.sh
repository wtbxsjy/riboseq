#!/usr/bin/env bash
set -euo pipefail

################################################################################
# Script: 00_prebuild_gencode_converter_images_with_timeout.sh
# Purpose: Pre-build Singularity images with extended timeout and retry logic
#
# Enhanced features:
# - Extended timeout settings
# - Automatic retry on failure
# - Support for wget/curl fallback download
# - Resume capability
################################################################################

usage() {
  cat <<'EOF'
Usage:
  00_prebuild_gencode_converter_images_with_timeout.sh [OPTIONS]

Options:
  --outdir DIR       Output directory for .sif files (default: ./containers)
  --cache-dir DIR    Singularity cache directory (default: system default)
  --timeout SECONDS  HTTP timeout in seconds (default: 3600)
  --retry COUNT      Number of retry attempts (default: 3)
  --method METHOD    Download method: singularity, wget, or curl (default: singularity)
  --force            Force rebuild even if image exists
  -h, --help         Show this help message

Environment Variables:
  SINGULARITY_CACHEDIR     Custom cache directory
  SINGULARITY_TMPDIR       Custom temp directory
  http_proxy               HTTP proxy server
  https_proxy              HTTPS proxy server

Examples:
  # Basic usage with extended timeout
  bash 00_prebuild_gencode_converter_images_with_timeout.sh --timeout 7200

  # Use wget for more reliable download
  bash 00_prebuild_gencode_converter_images_with_timeout.sh --method wget

  # With retry and custom output
  bash 00_prebuild_gencode_converter_images_with_timeout.sh \
    --outdir /data/containers --retry 5 --timeout 3600

  # Use proxy
  http_proxy=http://proxy:8080 bash 00_prebuild_gencode_converter_images_with_timeout.sh

EOF
}

OUTDIR="./containers"
CACHE_DIR=""
TIMEOUT=3600  # 1 hour default
RETRY_COUNT=3
METHOD="singularity"
FORCE=false

while [[ $# -gt 0 ]]; do
  case "$1" in
    --outdir) OUTDIR="$2"; shift 2;;
    --cache-dir) CACHE_DIR="$2"; shift 2;;
    --timeout) TIMEOUT="$2"; shift 2;;
    --retry) RETRY_COUNT="$2"; shift 2;;
    --method) METHOD="$2"; shift 2;;
    --force) FORCE=true; shift;;
    -h|--help) usage; exit 0;;
    *) echo "[ERROR] Unknown argument: $1"; usage; exit 2;;
  esac
done

# Validate method
if [[ ! "$METHOD" =~ ^(singularity|wget|curl)$ ]]; then
  echo "[ERROR] Invalid method: $METHOD (must be singularity, wget, or curl)"
  exit 2
fi

# Create output directory
mkdir -p "$OUTDIR"
OUTDIR="$(cd "$OUTDIR" && pwd)"

# Set cache directory if specified
if [[ -n "$CACHE_DIR" ]]; then
  mkdir -p "$CACHE_DIR"
  export SINGULARITY_CACHEDIR="$(cd "$CACHE_DIR" && pwd)"
  export SINGULARITY_TMPDIR="$SINGULARITY_CACHEDIR/tmp"
  mkdir -p "$SINGULARITY_TMPDIR"
  echo "[INFO] Using custom Singularity cache: $SINGULARITY_CACHEDIR"
fi

echo "========================================================================"
echo "  Pre-building Singularity images with extended timeout"
echo "========================================================================"
echo ""
echo "Output directory: $OUTDIR"
echo "Timeout: ${TIMEOUT}s"
echo "Retry count: $RETRY_COUNT"
echo "Download method: $METHOD"
echo ""

# Function: Download with singularity pull
download_with_singularity() {
  local url="$1"
  local output="$2"
  local attempt=1

  while [[ $attempt -le $RETRY_COUNT ]]; do
    echo "[INFO] Attempt $attempt/$RETRY_COUNT: Pulling with singularity..."

    # Set timeout environment variables
    export SINGULARITY_DOWNLOAD_TIMEOUT="$TIMEOUT"

    if timeout "${TIMEOUT}s" singularity pull --force --disable-cache "$output" "$url"; then
      echo "[OK] Download successful!"
      return 0
    else
      echo "[WARNING] Attempt $attempt failed"
      if [[ $attempt -lt $RETRY_COUNT ]]; then
        echo "[INFO] Retrying in 10 seconds..."
        sleep 10
      fi
      attempt=$((attempt + 1))
    fi
  done

  echo "[ERROR] All attempts failed"
  return 1
}

# Function: Download with wget (supports resume)
download_with_wget() {
  local url="$1"
  local output="$2"
  local temp_file="${output}.tmp"
  local attempt=1

  # Convert docker:// URL to HTTP URL for Galaxy Depot
  if [[ "$url" == docker://* ]]; then
    echo "[WARNING] wget cannot handle docker:// URLs directly"
    echo "[INFO] Falling back to singularity method"
    download_with_singularity "$url" "$output"
    return $?
  fi

  while [[ $attempt -le $RETRY_COUNT ]]; do
    echo "[INFO] Attempt $attempt/$RETRY_COUNT: Downloading with wget..."

    # wget supports resume with -c flag
    if wget --timeout="$TIMEOUT" \
           --tries=1 \
           --continue \
           --show-progress \
           -O "$temp_file" \
           "$url"; then
      mv "$temp_file" "$output"
      echo "[OK] Download successful!"
      return 0
    else
      echo "[WARNING] Attempt $attempt failed"
      if [[ $attempt -lt $RETRY_COUNT ]]; then
        echo "[INFO] Retrying in 10 seconds..."
        sleep 10
      fi
      attempt=$((attempt + 1))
    fi
  done

  rm -f "$temp_file"
  echo "[ERROR] All attempts failed"
  return 1
}

# Function: Download with curl (supports resume)
download_with_curl() {
  local url="$1"
  local output="$2"
  local temp_file="${output}.tmp"
  local attempt=1

  # Convert docker:// URL to HTTP URL for Galaxy Depot
  if [[ "$url" == docker://* ]]; then
    echo "[WARNING] curl cannot handle docker:// URLs directly"
    echo "[INFO] Falling back to singularity method"
    download_with_singularity "$url" "$output"
    return $?
  fi

  while [[ $attempt -le $RETRY_COUNT ]]; do
    echo "[INFO] Attempt $attempt/$RETRY_COUNT: Downloading with curl..."

    # curl supports resume with -C - flag
    if curl --max-time "$TIMEOUT" \
           --connect-timeout 60 \
           --retry 0 \
           --continue-at - \
           --progress-bar \
           -L -o "$temp_file" \
           "$url"; then
      mv "$temp_file" "$output"
      echo "[OK] Download successful!"
      return 0
    else
      echo "[WARNING] Attempt $attempt failed"
      if [[ $attempt -lt $RETRY_COUNT ]]; then
        echo "[INFO] Retrying in 10 seconds..."
        sleep 10
      fi
      attempt=$((attempt + 1))
    fi
  done

  rm -f "$temp_file"
  echo "[ERROR] All attempts failed"
  return 1
}

# Function: Smart download dispatcher
download_image() {
  local url="$1"
  local output="$2"

  case "$METHOD" in
    singularity)
      download_with_singularity "$url" "$output"
      ;;
    wget)
      download_with_wget "$url" "$output"
      ;;
    curl)
      download_with_curl "$url" "$output"
      ;;
  esac
}

# Image 1: Biopython for scripts 14 and 15
IMG_URL_BIOPYTHON="https://depot.galaxyproject.org/singularity/biopython:1.81"
IMG_NAME_BIOPYTHON="biopython_1.81.sif"
IMG_PATH_BIOPYTHON="$OUTDIR/$IMG_NAME_BIOPYTHON"

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
  if ! download_image "$IMG_URL_BIOPYTHON" "$IMG_PATH_BIOPYTHON"; then
    echo "[ERROR] Failed to download Biopython image"
    exit 1
  fi
  echo "[INFO] Size: $(du -h "$IMG_PATH_BIOPYTHON" | cut -f1)"
  echo ""
fi

# Image 2: Mulled container for script 16
IMG_URL_GENCODE="https://depot.galaxyproject.org/singularity/mulled-v2-8849acf39a43cdd6c839a369a74c0adc823e2f91:ab110436faf952a33575c64dd74615a84011450b-0"
IMG_NAME_GENCODE="gencode_orf_mapper_mulled.sif"
IMG_PATH_GENCODE="$OUTDIR/$IMG_NAME_GENCODE"

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
  if ! download_image "$IMG_URL_GENCODE" "$IMG_PATH_GENCODE"; then
    echo "[ERROR] Failed to download GENCODE ORF mapper image"
    exit 1
  fi
  echo "[INFO] Size: $(du -h "$IMG_PATH_GENCODE" | cut -f1)"
  echo ""
fi

# Verify images
echo "------------------------------------------------------------------------"
echo "Verifying images..."
echo "------------------------------------------------------------------------"
echo ""

echo "=== Image 1: Biopython ==="
if ! singularity exec "$IMG_PATH_BIOPYTHON" python3 --version; then
  echo "[ERROR] Verification failed: cannot run python3"
  exit 1
fi

if ! singularity exec "$IMG_PATH_BIOPYTHON" python3 -c "import Bio; print('Biopython version:', Bio.__version__)"; then
  echo "[ERROR] Verification failed: cannot import Biopython"
  exit 1
fi

echo ""
echo "=== Image 2: GENCODE ORF Mapper ==="
if ! singularity exec "$IMG_PATH_GENCODE" python3 --version; then
  echo "[ERROR] Verification failed: cannot run python3"
  exit 1
fi

if ! singularity exec "$IMG_PATH_GENCODE" python3 -c "import Bio; print('Biopython:', Bio.__version__)"; then
  echo "[ERROR] Verification failed: cannot import Biopython"
  exit 1
fi

if ! singularity exec "$IMG_PATH_GENCODE" python3 -c "import pandas; print('pandas:', pandas.__version__)"; then
  echo "[ERROR] Verification failed: cannot import pandas"
  exit 1
fi

if ! singularity exec "$IMG_PATH_GENCODE" bedtools --version 2>&1 | head -1; then
  echo "[WARNING] bedtools not found in container (may not be critical)"
fi

echo ""
echo "========================================================================"
echo "  ✓ Pre-build complete!"
echo "========================================================================"
echo ""
echo "Built images:"
ls -lh "$OUTDIR"/*.sif 2>/dev/null || echo "No images found"
echo ""
echo "------------------------------------------------------------------------"
echo "Next steps:"
echo "------------------------------------------------------------------------"
echo ""
echo "1. Transfer to your HPC cluster:"
echo "   rsync -avP $OUTDIR/ user@server:/path/to/riboseq/containers/"
echo ""
echo "2. Or compress and transfer:"
echo "   tar -czf gencode_images.tar.gz -C \"$OUTDIR\" ."
echo "   scp gencode_images.tar.gz user@server:/path/to/"
echo ""
echo "3. Run the converter scripts:"
echo "   bash 14_ribotish_to_gencode.sh --image biopython_1.81.sif ..."
echo "   bash 15_ribotricer_to_gencode.sh --image biopython_1.81.sif ..."
echo "   bash 16_gencode_orf_mapper.sh --image gencode_orf_mapper_mulled.sif ..."
echo ""
echo "========================================================================"
