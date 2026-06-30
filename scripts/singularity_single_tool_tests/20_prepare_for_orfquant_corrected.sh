#!/usr/bin/env bash
set -euo pipefail

# Mirrors: modules/local/prepare_for_orfquant/main.nf
# Regenerates for_ORFquant file with corrected P-site offsets
# Uses ORFquant's prepare_for_ORFquant() function with path_to_rl_cutoff_file

# Uses the same container as ORFquant
IMG_URL="docker://ghcr.io/orfquant-patched:latest"

usage() {
  cat <<'EOF'
Usage:
  20_prepare_for_orfquant_corrected.sh --sample ID --for-orfquant FILE --rl-cutoff FILE --annot FILE [OPTIONS]

Required:
  --sample         Sample ID (prefix)
  --for-orfquant   Original for_ORFquant file from RiboseQC
  --rl-cutoff      P-site offset file from 19_extract_rl_cutoff.sh
  --annot          RiboseQC annotation file (*_Rannot)

Options:
  --outdir    Output directory (default: ./out_prepare_for_orfquant_corrected)
  --image     Path to custom ORFquant container (default: ./containers/orfquant_patched.sif)

Env:
  BIND_EXTRA  Extra singularity binds (e.g. /mnt:/mnt)

Output:
  {sample}_corrected_for_ORFquant  - Regenerated for_ORFquant file with corrected offsets
EOF
}

SAMPLE=""
FOR_ORFQUANT=""
RL_CUTOFF=""
ANNOT=""
OUTDIR="./out_prepare_for_orfquant_corrected"
IMAGE=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --sample) SAMPLE="$2"; shift 2;;
    --for-orfquant) FOR_ORFQUANT="$2"; shift 2;;
    --rl-cutoff) RL_CUTOFF="$2"; shift 2;;
    --annot) ANNOT="$2"; shift 2;;
    --outdir) OUTDIR="$2"; shift 2;;
    --image) IMAGE="$2"; shift 2;;
    -h|--help) usage; exit 0;;
    *) echo "Unknown arg: $1"; usage; exit 2;;
  esac
done

if [[ -z "$SAMPLE" || -z "$FOR_ORFQUANT" || -z "$RL_CUTOFF" || -z "$ANNOT" ]]; then
  usage
  exit 2
fi

# Create output and container directories
mkdir -p "$OUTDIR" "./containers"
OUTDIR="$(cd "$OUTDIR" && pwd)"
WORKDIR="$(pwd)"

# Get absolute path
abspath() {
  python3 -c "import os,sys; print(os.path.abspath(sys.argv[1]))" "$1"
}

FOR_ORFQUANT="$(abspath "$FOR_ORFQUANT")"
RL_CUTOFF="$(abspath "$RL_CUTOFF")"
ANNOT="$(abspath "$ANNOT")"

# Build bind mounts
BIND_SPEC="$WORKDIR:$WORKDIR,$OUTDIR:$OUTDIR"
BIND_SPEC="$BIND_SPEC,$(dirname "$FOR_ORFQUANT"):$(dirname "$FOR_ORFQUANT")"
BIND_SPEC="$BIND_SPEC,$(dirname "$RL_CUTOFF"):$(dirname "$RL_CUTOFF")"
BIND_SPEC="$BIND_SPEC,$(dirname "$ANNOT"):$(dirname "$ANNOT")"

# Find ORFquant container
if [[ -n "$IMAGE" ]]; then
  IMG="$(abspath "$IMAGE")"
elif [[ -f "./containers/orfquant_patched.sif" ]]; then
  IMG="$(abspath ./containers/orfquant_patched.sif)"
elif [[ -f "$(dirname "${BASH_SOURCE[0]}")/../../orfquant_patched.sif" ]]; then
  IMG="$(abspath "$(dirname "${BASH_SOURCE[0]}")/../../orfquant_patched.sif")"
else
  echo "[ERROR] ORFquant container not found. Please specify --image or place orfquant_patched.sif in containers/"
  exit 1
fi

echo "[INFO] Using ORFquant container: $IMG"
echo "[INFO] Input for_ORFquant: $FOR_ORFQUANT"
echo "[INFO] RL cutoff file: $RL_CUTOFF"
echo "[INFO] Annotation: $ANNOT"
echo "[INFO] Output: $OUTDIR/${SAMPLE}_corrected_for_ORFquant"

singularity exec \
  --bind "$BIND_SPEC${BIND_EXTRA:+,$BIND_EXTRA}" \
  --pwd "$WORKDIR" \
  "$IMG" \
  bash -lc "
set -euo pipefail
cd '$OUTDIR'

cat <<'RSCRIPT' > script.R
library(ORFquant)

# Prepare for ORFquant with corrected P-site offsets
# This regenerates the for_ORFquant file using the optimal offsets

# Load original for_ORFquant to extract the RiboseQC_result path
original_data <- readRDS('$FOR_ORFQUANT')

# Extract the RiboseQC results file path
# The for_ORFquant file contains references to the results
riboseqc_result <- original_data\$psites_all_per_read_length

cat('Regenerating for_ORFquant with corrected P-site offsets...\\n')
cat('Using RL cutoff file: $RL_CUTOFF\\n')

# Call prepare_for_ORFquant with the rl_cutoff file
prepare_for_ORFquant(
    annotation_file = '$ANNOT',
    for_ORFquant_file = '$FOR_ORFQUANT',
    dest_name = '${SAMPLE}_corrected',
    path_to_rl_cutoff_file = '$RL_CUTOFF'
)

cat('Successfully generated ${SAMPLE}_corrected_for_ORFquant\\n')
RSCRIPT

Rscript script.R
"

echo "[OK] Corrected for_ORFquant generated: $OUTDIR/${SAMPLE}_corrected_for_ORFquant"
ls -la "$OUTDIR/${SAMPLE}_corrected_for_ORFquant"
