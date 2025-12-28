#!/usr/bin/env bash
set -euo pipefail

# Pipeline wrapper (single-tool chain)
# Start from a filtered BAM, run RiboseQC (prepareannotation + analysis), then ORFquant.
# Reuses the existing single-tool scripts in this directory:
#   02_riboseqc_prepareannotation.sh
#   03_riboseqc_analysis.sh
#   04_orfquant_run.sh

usage() {
  cat <<'EOF'
Usage:
  12_riboseqc_orfquant_from_filtered_bam.sh --sample ID --bam filtered.bam --gtf annot.gtf --fasta genome.fa [options]

Required:
  --sample    Sample ID (prefix)
  --bam       Filtered BAM (sorted + indexed recommended)
  --gtf       Annotation GTF (can be .gz)
  --fasta     Genome FASTA (can be .gz)

Options:
  --outdir        Output directory (default: ./out_riboseqc_orfquant)
  --cpus          Threads (default: 4)
  --fast-mode     TRUE|FALSE for RiboseQC_analysis (default: TRUE)
  --annotation    Existing RiboseQC annotation file (*_Rannot). If provided, skip prepareannotation.
  --orfquant-pkg  Local ORFquant source tar.gz (optional; avoids GitHub download)

Env:
  BIND_EXTRA  Extra singularity binds, comma-separated (e.g. /mnt:/mnt)

Outputs (inside --outdir):
  01_riboseqc_annot/      -> RiboseQC annotation (contains *_Rannot)
  02_riboseqc_analysis/   -> RiboseQC analysis (contains ${sample}_for_ORFquant)
  03_orfquant/            -> ORFquant results
EOF
}

SAMPLE=""
BAM=""
GTF=""
FASTA=""
OUTDIR="./out_riboseqc_orfquant"
CPUS=4
FAST_MODE="TRUE"
RANNOT=""
ORFQUANT_PKG=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --sample) SAMPLE="$2"; shift 2;;
    --bam) BAM="$2"; shift 2;;
    --gtf) GTF="$2"; shift 2;;
    --fasta) FASTA="$2"; shift 2;;
    --outdir) OUTDIR="$2"; shift 2;;
    --cpus) CPUS="$2"; shift 2;;
    --fast-mode) FAST_MODE="$2"; shift 2;;
    --annotation) RANNOT="$2"; shift 2;;
    --orfquant-pkg) ORFQUANT_PKG="$2"; shift 2;;
    -h|--help) usage; exit 0;;
    *) echo "Unknown arg: $1"; usage; exit 2;;
  esac
done

if [[ -z "$SAMPLE" || -z "$BAM" || -z "$GTF" || -z "$FASTA" ]]; then
  usage
  exit 2
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

mkdir -p "$OUTDIR"
OUTDIR="$(cd "$OUTDIR" && pwd)"

OUT_ANNOT="$OUTDIR/01_riboseqc_annot"
OUT_ANALYSIS="$OUTDIR/02_riboseqc_analysis"
OUT_ORFQUANT="$OUTDIR/03_orfquant"

mkdir -p "$OUT_ANNOT" "$OUT_ANALYSIS" "$OUT_ORFQUANT"

resolve_rannot() {
  local search_dir="$1"
  local gtf_path="$2"
  local base

  base="$(basename "$gtf_path")"
  base="${base%.gz}"
  base="${base%.gtf}"

  # Prefer a *_Rannot matching the GTF base name, if present.
  local preferred
  preferred="$(ls -1 "$search_dir" 2>/dev/null | awk -v b="$base" '$0 ~ b && $0 ~ /_Rannot$/ {print; exit}')" || true
  if [[ -n "$preferred" && -f "$search_dir/$preferred" ]]; then
    echo "$search_dir/$preferred"
    return 0
  fi

  # Otherwise, if only one *_Rannot exists, use it.
  local found
  mapfile -t found < <(find "$search_dir" -maxdepth 1 -type f -name "*_Rannot" | sort)
  if [[ ${#found[@]} -eq 1 ]]; then
    echo "${found[0]}"
    return 0
  fi

  # Ambiguous or not found.
  return 1
}

if [[ -z "$RANNOT" ]]; then
  echo "[INFO] Step 1/3: RiboseQC prepareannotation"
  bash "$SCRIPT_DIR/02_riboseqc_prepareannotation.sh" \
    --gtf "$GTF" \
    --fasta "$FASTA" \
    --outdir "$OUT_ANNOT"

  if ! RANNOT="$(resolve_rannot "$OUT_ANNOT" "$GTF")"; then
    echo "[ERROR] Cannot uniquely determine *_Rannot in: $OUT_ANNOT" >&2
    echo "        Please provide it explicitly with --annotation /path/to/*_Rannot" >&2
    echo "        Files found:" >&2
    find "$OUT_ANNOT" -maxdepth 1 -type f -name "*_Rannot" -print >&2 || true
    exit 2
  fi
else
  echo "[INFO] Using provided annotation: $RANNOT"
fi

echo "[INFO] Step 2/3: RiboseQC analysis"
bash "$SCRIPT_DIR/03_riboseqc_analysis.sh" \
  --sample "$SAMPLE" \
  --bam "$BAM" \
  --annotation "$RANNOT" \
  --fasta "$FASTA" \
  --outdir "$OUT_ANALYSIS" \
  --cpus "$CPUS" \
  --fast-mode "$FAST_MODE"

FOR_ORFQUANT="$OUT_ANALYSIS/${SAMPLE}_for_ORFquant"
if [[ ! -f "$FOR_ORFQUANT" ]]; then
  echo "[ERROR] Missing RiboseQC output for ORFquant: $FOR_ORFQUANT" >&2
  echo "        Check RiboseQC outputs in: $OUT_ANALYSIS" >&2
  exit 2
fi

echo "[INFO] Step 3/3: ORFquant"
ORFQUANT_ARGS=()
if [[ -n "$ORFQUANT_PKG" ]]; then
  ORFQUANT_ARGS+=(--orfquant-pkg "$ORFQUANT_PKG")
fi

bash "$SCRIPT_DIR/04_orfquant_run.sh" \
  --sample "$SAMPLE" \
  --for-orfquant "$FOR_ORFQUANT" \
  --annotation "$RANNOT" \
  --fasta "$FASTA" \
  --cpus "$CPUS" \
  --outdir "$OUT_ORFQUANT" \
  "${ORFQUANT_ARGS[@]}"

echo "[OK] Done. Outputs:"
echo "  - RiboseQC annotation: $OUT_ANNOT"
echo "  - RiboseQC analysis:   $OUT_ANALYSIS"
echo "  - ORFquant:            $OUT_ORFQUANT"
