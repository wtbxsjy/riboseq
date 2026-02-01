#!/usr/bin/env bash
set -euo pipefail

# Wraps: scripts/unify_orf_predictions.py
# Unifies ORF predictions from Ribo-TISH, Ribotricer, and ORFquant

# We use a standard python image and install dependencies (biopython, pyfaidx) at runtime
# to avoid needing a specific large container.
IMG_URL="https://depot.galaxyproject.org/singularity/python:3.9"

usage() {
  cat <<'EOF'
Usage:
  17_unify_predictions.sh --gtf ref.gtf --fasta genome.fa --output prefix [OPTIONS]

Required:
  --gtf        Reference GTF file
  --fasta      Genome FASTA file
  --output     Output prefix (e.g. ./unified_results/my_study)

Input Options (at least one required):
  --ribotish   "file1 file2 ..."  (Space-separated list of Ribo-TISH *_pred.txt)
  --ribotricer "file1 file2 ..."  (Space-separated list of Ribotricer *.tsv)
  --orfquant   "file1 file2 ..."  (Space-separated list of ORFquant *.gtf)

P-site Statistics Options (for unified P-site stats from RiboseQC):
  --bedgraph-dir  Directory containing RiboseQC P-site bedgraph files
  --sample-list   Comma-separated list of sample names (matching bedgraph file prefixes)

Other Options:
  --min-len    Minimum AA length (default: 10)
  --outdir     Output directory for logs/scripts (default: ./out_unify)
  --image      Path to Singularity image (default: auto-pull python:3.9)
  --cpus       Number of threads (default: 1)

Env:
  BIND_EXTRA   Extra singularity binds (e.g. /mnt:/mnt)

Examples:
  # Basic unification
  17_unify_predictions.sh \
    --gtf gencode.gtf --fasta genome.fa \
    --ribotish "s1_pred.txt s2_pred.txt" \
    --ribotricer "s1_orfs.tsv" \
    --output ./results/unified

  # With P-site statistics from RiboseQC
  17_unify_predictions.sh \
    --gtf gencode.gtf --fasta genome.fa \
    --ribotish "s1_pred.txt" --orfquant "s1_orfquant.gtf" \
    --bedgraph-dir ./out_riboseqc_analysis \
    --sample-list "s1,s2" \
    --output ./results/unified_with_psite_stats
EOF
}

GTF=""
FASTA=""
OUTPUT_PREFIX=""
RIBOTISH_FILES=""
RIBOTRICER_FILES=""
ORFQUANT_FILES=""
BEDGRAPH_DIR=""
SAMPLE_LIST=""
MIN_LEN=10
OUTDIR="./out_unify"
IMAGE=""
CPUS=1

while [[ $# -gt 0 ]]; do
  case "$1" in
    --gtf) GTF="$2"; shift 2;;
    --fasta) FASTA="$2"; shift 2;;
    --output) OUTPUT_PREFIX="$2"; shift 2;;
    --ribotish) RIBOTISH_FILES="$2"; shift 2;;
    --ribotricer) RIBOTRICER_FILES="$2"; shift 2;;
    --orfquant) ORFQUANT_FILES="$2"; shift 2;;
    --bedgraph-dir) BEDGRAPH_DIR="$2"; shift 2;;
    --sample-list) SAMPLE_LIST="$2"; shift 2;;
    --min-len) MIN_LEN="$2"; shift 2;;
    --outdir) OUTDIR="$2"; shift 2;;
    --image) IMAGE="$2"; shift 2;;
    --cpus) CPUS="$2"; shift 2;;
    -h|--help) usage; exit 0;;
    *) echo "Unknown arg: $1"; usage; exit 2;;
  esac
done

if [[ -z "$GTF" || -z "$FASTA" || -z "$OUTPUT_PREFIX" ]]; then
  echo "[ERROR] Missing required arguments."
  usage
  exit 2
fi

# Ensure output directory exists
OUTPUT_DIR_PATH="$(dirname "$OUTPUT_PREFIX")"
mkdir -p "$OUTPUT_DIR_PATH" "$OUTDIR" "./containers"

# Absolute paths
abspath() { python3 -c "import os,sys; print(os.path.abspath(sys.argv[1]))" "$1"; }

GTF="$(abspath "$GTF")"
FASTA="$(abspath "$FASTA")"
# Handle output prefix relative to current dir, convert to absolute
OUTPUT_PREFIX_ABS="$(abspath "$(dirname "$OUTPUT_PREFIX")")/$(basename "$OUTPUT_PREFIX")"

OUTDIR="$(abspath "$OUTDIR")"
WORKDIR="$(pwd)"

# Binds
BIND_SPEC=""
add_bind() {
  local p="$1"
  [[ -z "$p" ]] && return 0
  if [[ -e "$p" ]]; then
    p="$(abspath "$p")"
    if [[ -d "$p" ]]; then
        local entry="${p}:${p}"
    else
        local d
        d="$(dirname "$p")"
        local entry="${d}:${d}"
    fi
    # Simple deduplication check could be added here
    BIND_SPEC="${BIND_SPEC:+$BIND_SPEC,}${entry}"
  fi
}

add_bind "$WORKDIR"
add_bind "$OUTDIR"
add_bind "$GTF"
add_bind "$FASTA"
add_bind "$(dirname "$OUTPUT_PREFIX_ABS")"

# Bind bedgraph directory if provided
if [[ -n "$BEDGRAPH_DIR" ]]; then
    BEDGRAPH_DIR_ABS="$(abspath "$BEDGRAPH_DIR")"
    add_bind "$BEDGRAPH_DIR_ABS"
fi

# Bind inputs
for f in $RIBOTISH_FILES $RIBOTRICER_FILES $ORFQUANT_FILES; do
    add_bind "$f"
done

# Bind script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$(dirname "$SCRIPT_DIR")")"
SCRIPTS_ROOT="$PROJECT_ROOT/scripts"
add_bind "$SCRIPTS_ROOT"

# Container
pull_img() {
  local url="$1"
  local base="$(basename "$url" | sed 's/:/_/g')"
  local sif="$(pwd)/containers/${base}.sif"
  if [[ ! -f "$sif" ]]; then
    echo "[INFO] Pulling $url..."
    singularity pull --disable-cache --force "$sif" "$url"
  fi
  echo "$sif"
}

if [[ -n "$IMAGE" ]]; then
  IMG="$(abspath "$IMAGE")"
else
  IMG="$(pull_img "$IMG_URL")"
fi

echo "[INFO] Running Unified Prediction..."
echo "[INFO] Output: $OUTPUT_PREFIX_ABS"

# Construct command args
CMD_ARGS="--gtf '$GTF' --fasta '$FASTA' --output '$OUTPUT_PREFIX_ABS' --min_len $MIN_LEN --threads $CPUS"

if [[ -n "$RIBOTISH_FILES" ]]; then
    CMD_ARGS="$CMD_ARGS --ribotish $RIBOTISH_FILES"
fi
if [[ -n "$RIBOTRICER_FILES" ]]; then
    CMD_ARGS="$CMD_ARGS --ribotricer $RIBOTRICER_FILES"
fi
if [[ -n "$ORFQUANT_FILES" ]]; then
    CMD_ARGS="$CMD_ARGS --orfquant $ORFQUANT_FILES"
fi

# Add bedgraph options for P-site statistics
if [[ -n "$BEDGRAPH_DIR" ]]; then
    CMD_ARGS="$CMD_ARGS --bedgraph-dir '$BEDGRAPH_DIR_ABS'"
fi
if [[ -n "$SAMPLE_LIST" ]]; then
    CMD_ARGS="$CMD_ARGS --sample-list '$SAMPLE_LIST'"
fi

singularity exec \
  --bind "$BIND_SPEC${BIND_EXTRA:+,$BIND_EXTRA}" \
  --pwd "$WORKDIR" \
  "$IMG" \
  bash -c "
    set -e
    # Install deps if missing (assuming writable tmp or just trying)
    # Since singularity images are read-only, we install to a user location or venv
    # But usually pip install --user works if HOME is bound or we set PYTHONUSERBASE
    
    export PYTHONUSERBASE=$OUTDIR/pylibs
    export PATH=\$PYTHONUSERBASE/bin:\$PATH
    mkdir -p \$PYTHONUSERBASE
    
    echo '[INFO] Installing dependencies (biopython, pyfaidx)...'
    pip install --user --quiet --no-warn-script-location biopython pyfaidx
    
    echo '[INFO] Running script...'
    python3 '$SCRIPTS_ROOT/unify_orf_predictions.py' $CMD_ARGS
    
    echo '[INFO] Versions:'
    python3 --version
    pip show biopython | grep Version || true
    pip show pyfaidx | grep Version || true
"

echo "[OK] Unified prediction complete."
ls -l "${OUTPUT_PREFIX_ABS}".*
