#!/usr/bin/env bash
set -euo pipefail

# Wraps: scripts/classify_orfs_wrapper.py
# Runs ORF classification using one of three modes: gencode, orfquant, orf_type

usage() {
  cat <<'EOF'
Usage:
  18_classify_orfs.sh --mode MODE --input INPUT --output-dir DIR [OPTIONS]

Required:
  --mode        Classification mode: 'gencode', 'orfquant', or 'orf_type'
  --input       Input prefix or file (output from 17_unify_predictions.sh)
                e.g. ./results/unified_study (will look for .bed, .gtf, .metadata.tsv)
  --output-dir  Output directory for results

Ref Options (Required based on mode):
  --gtf         Reference GTF (Required for 'orfquant', 'orf_type')
  --fasta       Reference FASTA (Required for 'gencode' if input is prefix)
  --ensembl-dir Ensembl directory (Required for 'gencode')
                Contains PROTEOME_FASTA, TRANSCRIPTOME_GTF etc.

Other Options:
  --image       Override Singularity image path
  --cpus        CPUs (default: 1)
  --outdir      Temp/Log directory (default: ./out_classify)

Examples:
  # Gencode Mode
  18_classify_orfs.sh --mode gencode --input ./unified/res \
    --ensembl-dir ./refs/ensembl --output-dir ./class_res

  # ORFquant Mode
  18_classify_orfs.sh --mode orfquant --input ./unified/res \
    --gtf ./refs/gen.gtf --output-dir ./class_res

  # ORFtype Mode
  18_classify_orfs.sh --mode orf_type --input ./unified/res \
    --gtf ./refs/gen.gtf --output-dir ./class_res
EOF
}

MODE=""
INPUT=""
OUTPUT_DIR=""
GTF=""
FASTA=""
ENSEMBL_DIR=""
IMAGE=""
CPUS=1
OUTDIR="./out_classify"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --mode) MODE="$2"; shift 2;;
    --input) INPUT="$2"; shift 2;;
    --output-dir) OUTPUT_DIR="$2"; shift 2;;
    --gtf) GTF="$2"; shift 2;;
    --fasta) FASTA="$2"; shift 2;;
    --ensembl-dir) ENSEMBL_DIR="$2"; shift 2;;
    --image) IMAGE="$2"; shift 2;;
    --cpus) CPUS="$2"; shift 2;;
    --outdir) OUTDIR="$2"; shift 2;;
    -h|--help) usage; exit 0;;
    *) echo "Unknown arg: $1"; usage; exit 2;;
  esac
done

if [[ -z "$MODE" || -z "$INPUT" || -z "$OUTPUT_DIR" ]]; then
  echo "[ERROR] Missing required arguments."
  usage
  exit 2
fi

mkdir -p "$OUTPUT_DIR" "$OUTDIR" "./containers"

# Absolute paths
abspath() { python3 -c "import os,sys; print(os.path.abspath(sys.argv[1]))" "$1"; }

INPUT="$(abspath "$INPUT")"
OUTPUT_DIR="$(abspath "$OUTPUT_DIR")"
OUTDIR="$(abspath "$OUTDIR")"
WORKDIR="$(pwd)"

if [[ -n "$GTF" ]]; then GTF="$(abspath "$GTF")"; fi
if [[ -n "$FASTA" ]]; then FASTA="$(abspath "$FASTA")"; fi
if [[ -n "$ENSEMBL_DIR" ]]; then ENSEMBL_DIR="$(abspath "$ENSEMBL_DIR")"; fi

# Determine Image based on Mode
IMG_URL=""
case "$MODE" in
    gencode)
        # Image with python, pandas, bedtools, biopython
        # Same as 16_gencode_orf_mapper.sh
        IMG_URL="https://depot.galaxyproject.org/singularity/mulled-v2-8849acf39a43cdd6c839a369a74c0adc823e2f91:ab110436faf952a33575c64dd74615a84011450b-0"
        if [[ -z "$ENSEMBL_DIR" ]]; then echo "[ERROR] --ensembl-dir required for gencode mode"; exit 1; fi
        ;;
    orfquant)
        # Image with R, ORFquant
        # Same as 04_orfquant_run.sh
        IMG_URL="https://depot.galaxyproject.org/singularity/orfquant:1.1.0--r40_1"
        if [[ -z "$GTF" ]]; then echo "[ERROR] --gtf required for orfquant mode"; exit 1; fi
        ;;
    orf_type)
        # Image with python standard lib (and maybe multiprocessing)
        # python:3.9 is sufficient
        IMG_URL="https://depot.galaxyproject.org/singularity/python:3.9"
        if [[ -z "$GTF" ]]; then echo "[ERROR] --gtf required for orf_type mode"; exit 1; fi
        ;;
    *)
        echo "[ERROR] Unknown mode: $MODE"
        exit 1
        ;;
esac

# Container pull
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
    BIND_SPEC="${BIND_SPEC:+$BIND_SPEC,}${entry}"
  fi
}

add_bind "$WORKDIR"
add_bind "$OUTDIR"
add_bind "$OUTPUT_DIR"
# Bind input directory (assuming input is prefix or file)
if [[ -d "$INPUT" ]]; then
    add_bind "$INPUT"
else
    add_bind "$(dirname "$INPUT")"
fi

add_bind "$GTF"
add_bind "$FASTA"
add_bind "$ENSEMBL_DIR"

# Bind scripts
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$(dirname "$SCRIPT_DIR")")"
SCRIPTS_ROOT="$PROJECT_ROOT/scripts"
add_bind "$SCRIPTS_ROOT"

echo "[INFO] Running Classification ($MODE)..."
echo "[INFO] Image: $IMG"

# Construct command
CMD_ARGS="--mode '$MODE' --input '$INPUT' --output_dir '$OUTPUT_DIR' --cpus '$CPUS'"
if [[ -n "$GTF" ]]; then CMD_ARGS="$CMD_ARGS --gtf '$GTF'"; fi
if [[ -n "$FASTA" ]]; then CMD_ARGS="$CMD_ARGS --fasta '$FASTA'"; fi
if [[ -n "$ENSEMBL_DIR" ]]; then CMD_ARGS="$CMD_ARGS --ensembl_dir '$ENSEMBL_DIR'"; fi

# Execute
# Note: For orf_type mode, we might need to pip install nothing as it uses stdlib.
# For gencode mode, image has deps.
# For orfquant mode, image has R deps.

singularity exec \
  --bind "$BIND_SPEC${BIND_EXTRA:+,$BIND_EXTRA}" \
  --pwd "$WORKDIR" \
  "$IMG" \
  bash -c "
    set -e
    # Run the wrapper
    # The wrapper will call python3 or Rscript depending on mode.
    # We assume 'python3' and 'Rscript' are in PATH of the container.
    
    python3 '$SCRIPTS_ROOT/classify_orfs_wrapper.py' $CMD_ARGS
"

echo "[OK] Classification complete."
ls -l "$OUTPUT_DIR"
