#!/usr/bin/env bash
set -euo pipefail

# ORFquant-specific annotation preparation
# Based on the provided Nextflow example pipeline:
#   - Create .2bit from genome FASTA (faToTwoBit)
#   - Run ORFquant::prepare_annotation_files()
# NOTE: ORFquant Rannot is NOT compatible with RiboseQC Rannot.

# Avoid noisy locale warnings inside containers
export SINGULARITYENV_LANG=${SINGULARITYENV_LANG:-C}
export SINGULARITYENV_LC_ALL=${SINGULARITYENV_LC_ALL:-C}
export APPTAINERENV_LANG=${APPTAINERENV_LANG:-C}
export APPTAINERENV_LC_ALL=${APPTAINERENV_LC_ALL:-C}

# ORFquant container from the example (R 4.0)
ORFQUANT_IMG_URL="https://depot.galaxyproject.org/singularity/orfquant:1.1.0--r40_1"
# faToTwoBit container from the example (ORAS). Override with --fatotwobit-container if needed.
FATOTWOBIT_ORAS="oras://community.wave.seqera.io/library/ucsc-fatotwobit:482--1d5005b012bd3271"

usage() {
  cat <<'EOF'
Usage:
  13_orfquant_prepareannotation.sh --gtf annot.gtf[.gz] --fasta genome.fa[.gz] --outdir DIR [options]

Required:
  --gtf      Annotation GTF (can be .gz)
  --fasta    Genome FASTA (can be .gz)
  --outdir   Output directory (will contain *_Rannot)

Options:
  --twobit               Existing .2bit file (skip faToTwoBit)
  --species              Scientific name (e.g. Mus.musculus). Default: inferred from GTF filename when possible, else 'Homo.sapiens'
  --annotation-name       Annotation name. Default: inferred from GTF filename base
  --container             ORFquant SIF path (optional). If not set, pulls/uses orfquant:1.1.0--r40_1
  --fatotwobit-container  Container reference or SIF for faToTwoBit (default: ORAS reference from example)

Env:
  BIND_EXTRA  Extra container binds, comma-separated (e.g. /mnt:/mnt)

Outputs:
  - *_Rannot in --outdir
EOF
}

GTF=""
FASTA=""
OUTDIR=""
TWOBIT=""
SPECIES=""
ANN_NAME=""
CONTAINER_SIF=""
FATOTWOBIT_CONTAINER=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --gtf) GTF="$2"; shift 2;;
    --fasta) FASTA="$2"; shift 2;;
    --outdir) OUTDIR="$2"; shift 2;;
    --twobit) TWOBIT="$2"; shift 2;;
    --species) SPECIES="$2"; shift 2;;
    --annotation-name) ANN_NAME="$2"; shift 2;;
    --container) CONTAINER_SIF="$2"; shift 2;;
    --fatotwobit-container) FATOTWOBIT_CONTAINER="$2"; shift 2;;
    -h|--help) usage; exit 0;;
    *) echo "Unknown arg: $1"; usage; exit 2;;
  esac
done

if [[ -z "$GTF" || -z "$FASTA" || -z "$OUTDIR" ]]; then
  usage
  exit 2
fi

if [[ ! -f "$GTF" ]]; then
  echo "[ERROR] --gtf not found: $GTF" >&2
  exit 2
fi
if [[ ! -f "$FASTA" ]]; then
  echo "[ERROR] --fasta not found: $FASTA" >&2
  exit 2
fi

if [[ -n "$TWOBIT" && ! -f "$TWOBIT" ]]; then
  echo "[ERROR] --twobit not found: $TWOBIT" >&2
  exit 2
fi

if [[ -n "$CONTAINER_SIF" && ! -f "$CONTAINER_SIF" ]]; then
  echo "[ERROR] --container not found: $CONTAINER_SIF" >&2
  exit 2
fi

if [[ -z "$FATOTWOBIT_CONTAINER" ]]; then
  FATOTWOBIT_CONTAINER="$FATOTWOBIT_ORAS"
fi

detect_runtime() {
  if command -v apptainer >/dev/null 2>&1; then
    echo "apptainer"
    return 0
  fi
  if command -v singularity >/dev/null 2>&1; then
    echo "singularity"
    return 0
  fi
  echo "[ERROR] Neither 'apptainer' nor 'singularity' is available on PATH." >&2
  exit 2
}

RUNTIME="$(detect_runtime)"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORKDIR="$(pwd)"
mkdir -p "$OUTDIR" "./containers"
OUTDIR="$(cd "$OUTDIR" && pwd)"

pull_img() {
  local url="$1"
  local base
  base="$(basename "$url")"
  base="${base//:/_}"
  local sif="$WORKDIR/containers/${base}.sif"
  if [[ ! -f "$sif" ]]; then
    "$RUNTIME" pull --disable-cache --force "$sif" "$url"
  fi
  echo "$sif"
}

if [[ -n "$CONTAINER_SIF" ]]; then
  ORFQ_IMG="$(python3 -c 'import os,sys; print(os.path.abspath(sys.argv[1]))' "$CONTAINER_SIF")"
else
  ORFQ_IMG="$(pull_img "$ORFQUANT_IMG_URL")"
fi

# Infer species/annotation_name from GTF filename if not provided (same as example)
gtf_base="$(basename "$GTF")"
gtf_base="${gtf_base%.gz}"
if [[ -z "$ANN_NAME" ]]; then
  ANN_NAME="${gtf_base%.gtf}"
fi
if [[ -z "$SPECIES" ]]; then
  # Example logic: species <- sub("_",".",sub("\\..+","",gtf_file_name))
  # Here: take substring before first '.' and replace '_' with '.'
  species_guess="${gtf_base%%.*}"
  species_guess="${species_guess//_/.}"
  if [[ -n "$species_guess" ]]; then
    SPECIES="$species_guess"
  else
    SPECIES="Homo.sapiens"
  fi
fi

# Build binds (input dirs + outdir + workdir)
collect_binds() {
  local items=()
  add() {
    local p="$1"
    [[ -z "$p" ]] && return 0
    if [[ -e "$p" ]]; then
      local d
      if [[ -d "$p" ]]; then d="$p"; else d="$(cd "$(dirname "$p")" && pwd)"; fi
      items+=("$d:$d")
    fi
  }
  add "$WORKDIR"
  add "$OUTDIR"
  add "$GTF"
  add "$FASTA"
  add "$TWOBIT"

  if [[ -n "${BIND_EXTRA:-}" ]]; then
    IFS=',' read -r -a extra <<<"$BIND_EXTRA"
    for b in "${extra[@]}"; do
      [[ -n "$b" ]] && items+=("$b")
    done
  fi

  python3 - "$@" <<'PY'
import sys
seen=set(); out=[]
for b in sys.argv[1:]:
    if b and b not in seen:
        seen.add(b); out.append(b)
print(",".join(out))
PY
}

BIND_SPEC="$(collect_binds "${WORKDIR}:${WORKDIR}" "${OUTDIR}:${OUTDIR}" )"

# Prepare twobit if needed
TWOBIT_ABS=""
if [[ -n "$TWOBIT" ]]; then
  TWOBIT_ABS="$(python3 -c 'import os,sys; print(os.path.abspath(sys.argv[1]))' "$TWOBIT")"
else
  name="$(basename "$FASTA")"
  name="${name%.gz}"
  name="${name%.fa}"
  name="${name%.fasta}"
  TWOBIT_ABS="$OUTDIR/${name}.2bit"

  if [[ -f "$TWOBIT_ABS" ]]; then
    echo "[INFO] faToTwoBit skipped (exists): $TWOBIT_ABS"
  elif command -v faToTwoBit >/dev/null 2>&1; then
    echo "[INFO] Creating .2bit with host faToTwoBit"
    if [[ "$FASTA" == *.gz ]]; then
      gunzip -c "$FASTA" | faToTwoBit stdin "$TWOBIT_ABS"
    else
      faToTwoBit "$FASTA" "$TWOBIT_ABS"
    fi
  else
    echo "[INFO] Creating .2bit with container faToTwoBit: $FATOTWOBIT_CONTAINER"
    # Run faToTwoBit in its container; write output into OUTDIR
    "$RUNTIME" exec \
      --bind "$BIND_SPEC" \
      --pwd "$WORKDIR" \
      "$FATOTWOBIT_CONTAINER" \
      bash -lc 'set -euo pipefail; if [[ "'$FASTA'" == *.gz ]]; then gunzip -c "'$FASTA'" > "'$OUTDIR'/__genome.fa"; faToTwoBit "'$OUTDIR'/__genome.fa" "'$TWOBIT_ABS'"; rm -f "'$OUTDIR'/__genome.fa"; else faToTwoBit "'$FASTA'" "'$TWOBIT_ABS'"; fi'
  fi
fi

if [[ ! -f "$TWOBIT_ABS" ]]; then
  echo "[ERROR] Failed to create/find .2bit: $TWOBIT_ABS" >&2
  echo "        Provide an existing one with --twobit" >&2
  exit 2
fi

# Run ORFquant prepare_annotation_files inside ORFquant container
cd "$OUTDIR"
cat > run_orfquant_prepareannotation.R <<'RS'
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) {
  stop("Expected 5 args: gtf fasta twobit species ann_name")
}
gtf <- args[[1]]
fasta <- args[[2]]
twobit <- args[[3]]
species <- args[[4]]
ann_name <- args[[5]]

suppressPackageStartupMessages(library(ORFquant))
suppressPackageStartupMessages(library(magrittr))

prepare_annotation_files(
  annotation_directory = ".",
  twobit_file = twobit,
  gtf_file = gtf,
  scientific_name = species,
  annotation_name = ann_name,
  export_bed_tables_TxDb = TRUE,
  forge_BSgenome = FALSE,
  genome_seq = fasta
)
RS

echo "[INFO] ORFquant annotation: species=$SPECIES ann_name=$ANN_NAME"
"$RUNTIME" exec \
  --bind "$BIND_SPEC" \
  --pwd "$WORKDIR" \
  "$ORFQ_IMG" \
  Rscript --vanilla "$OUTDIR/run_orfquant_prepareannotation.R" \
    "$(python3 -c 'import os,sys; print(os.path.abspath(sys.argv[1]))' "$GTF")" \
    "$(python3 -c 'import os,sys; print(os.path.abspath(sys.argv[1]))' "$FASTA")" \
    "$TWOBIT_ABS" \
    "$SPECIES" \
    "$ANN_NAME"

if ! ls -1 *_Rannot >/dev/null 2>&1; then
  echo "[ERROR] ORFquant annotation did not produce *_Rannot in: $OUTDIR" >&2
  ls -la >&2 || true
  exit 2
fi

echo "[OK] ORFquant annotation ready in: $OUTDIR"
ls -la *_Rannot || true
