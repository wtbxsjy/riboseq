#!/usr/bin/env bash
set -euo pipefail

# Pipeline wrapper (single-tool chain)
# Start from an input BAM, run sORF-style BAM filtering, then RiboseQC (prepareannotation + analysis), then ORFquant.
# Reuses the existing single-tool scripts in this directory:
#   01_sorf_bam_filter.sh
#   02_riboseqc_prepareannotation.sh
#   03_riboseqc_analysis.sh
#   04_orfquant_run.sh

usage() {
  cat <<'EOF'
Usage:
  12_riboseqc_orfquant_from_filtered_bam.sh --sample ID --bam in.bam --fai genome.fa.fai --gtf annot.gtf --fasta genome.fa [options]

Required:
  --sample    Sample ID (prefix)
  --bam       Input BAM (sorted + indexed recommended)
  --fai       FASTA index (.fai) for contig list (used by filtering step)
  --gtf       Annotation GTF (can be .gz)
  --fasta     Genome FASTA (can be .gz)

Options:
  --outdir        Output directory (default: ./out_riboseqc_orfquant)
  --cpus          Threads (default: 4)
  --fast-mode     TRUE|FALSE for RiboseQC_analysis (default: TRUE)
  --resume        TRUE|FALSE (default: TRUE). If TRUE, skip steps with existing expected outputs.
  --skip-filter   TRUE|FALSE (default: FALSE). If TRUE, skip filtering and treat --bam as already filtered.
  --skip-orfquant TRUE|FALSE (default: FALSE). If TRUE, skip ORFquant and only produce RiboseQC outputs.
  --unique-mode   auto|nh|mapq for filtering (default: auto)
  --mapq          MAPQ threshold for filtering (default: 60)
  --len-min       read length min for filtering (default: 28)
  --len-max       read length max for filtering (default: 30)
  --exclude-regex contig regex to EXCLUDE for filtering (default: pipeline-like contig excludes)
  --annotation    Existing RiboseQC annotation file (*_Rannot). If provided, skip prepareannotation.
  --orfquant-pkg  Local ORFquant source tar.gz (optional; avoids GitHub download)
  --orfquant-sif  ORFquant container SIF path (optional; if set, passed to 04_orfquant_run.sh --container)
  --orfquant-annotation Existing ORFquant annotation (*_Rannot). If not set, will be generated for ORFquant.

Env:
  BIND_EXTRA  Extra singularity binds, comma-separated (e.g. /mnt:/mnt)

Outputs (inside --outdir):
  00_sorf_filter/         -> filtered BAM (${sample}.sorf.filtered.bam)
  01_riboseqc_annot/      -> RiboseQC annotation (contains *_Rannot)
  02_riboseqc_analysis/   -> RiboseQC analysis (contains ${sample}_for_ORFquant)
  03_orfquant/            -> ORFquant results
EOF
}

SAMPLE=""
BAM=""
FAI=""
GTF=""
FASTA=""
OUTDIR="./out_riboseqc_orfquant"
CPUS=4
FAST_MODE="TRUE"
RANNOT=""
ORFQUANT_PKG=""
ORFQUANT_SIF=""
ORFQUANT_RANNOT=""

SKIP_FILTER="FALSE"
SKIP_ORFQUANT="FALSE"
RESUME="TRUE"
UNIQUE_MODE="auto"
MAPQ=60
LEN_MIN=28
LEN_MAX=30
EXCLUDE_REGEX='^(chr)?(M|MT|Mt|chrM|chrMT|chrMt|ChrM|ChrMT|ChrMt)$|^(chr)?(C|CP|Pt|chrC|chrCP|chrPt|ChrC|ChrCP|ChrPt)$|^chrUn_.*|.*_random$|.*_alt$|.*_fix$'

die_missing_value() {
  local opt="$1"
  echo "[ERROR] Missing value for: $opt" >&2
  usage >&2
  exit 2
}

parse_bool_or_flag_true() {
  # Usage: parse_bool_or_flag_true --opt_name "$@"; echoes value and returns shift count via global PARSE_SHIFT
  # Behavior:
  #   --opt TRUE|FALSE  -> value=TRUE/FALSE, shift=2
  #   --opt             -> value=TRUE, shift=1
  local opt="$1"; shift
  if [[ $# -ge 2 && "${2:-}" != --* ]]; then
    echo "$2"
    PARSE_SHIFT=2
  else
    echo "TRUE"
    PARSE_SHIFT=1
  fi
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --sample) [[ $# -ge 2 ]] || die_missing_value "$1"; SAMPLE="$2"; shift 2;;
    --bam) [[ $# -ge 2 ]] || die_missing_value "$1"; BAM="$2"; shift 2;;
    --fai) [[ $# -ge 2 ]] || die_missing_value "$1"; FAI="$2"; shift 2;;
    --gtf) [[ $# -ge 2 ]] || die_missing_value "$1"; GTF="$2"; shift 2;;
    --fasta) [[ $# -ge 2 ]] || die_missing_value "$1"; FASTA="$2"; shift 2;;
    --outdir) [[ $# -ge 2 ]] || die_missing_value "$1"; OUTDIR="$2"; shift 2;;
    --cpus) [[ $# -ge 2 ]] || die_missing_value "$1"; CPUS="$2"; shift 2;;
    --fast-mode) [[ $# -ge 2 ]] || die_missing_value "$1"; FAST_MODE="$2"; shift 2;;
    --resume) [[ $# -ge 2 ]] || die_missing_value "$1"; RESUME="$2"; shift 2;;
    --skip-filter)
      SKIP_FILTER="$(parse_bool_or_flag_true "$1" "$@")"
      shift "$PARSE_SHIFT";;
    --skip-orfquant)
      SKIP_ORFQUANT="$(parse_bool_or_flag_true "$1" "$@")"
      shift "$PARSE_SHIFT";;
    --unique-mode) [[ $# -ge 2 ]] || die_missing_value "$1"; UNIQUE_MODE="$2"; shift 2;;
    --mapq) [[ $# -ge 2 ]] || die_missing_value "$1"; MAPQ="$2"; shift 2;;
    --len-min) [[ $# -ge 2 ]] || die_missing_value "$1"; LEN_MIN="$2"; shift 2;;
    --len-max) [[ $# -ge 2 ]] || die_missing_value "$1"; LEN_MAX="$2"; shift 2;;
    --exclude-regex) [[ $# -ge 2 ]] || die_missing_value "$1"; EXCLUDE_REGEX="$2"; shift 2;;
    --annotation) [[ $# -ge 2 ]] || die_missing_value "$1"; RANNOT="$2"; shift 2;;
    --orfquant-pkg) [[ $# -ge 2 ]] || die_missing_value "$1"; ORFQUANT_PKG="$2"; shift 2;;
    --orfquant-sif) [[ $# -ge 2 ]] || die_missing_value "$1"; ORFQUANT_SIF="$2"; shift 2;;
    --orfquant-annotation) [[ $# -ge 2 ]] || die_missing_value "$1"; ORFQUANT_RANNOT="$2"; shift 2;;
    -h|--help) usage; exit 0;;
    *) echo "Unknown arg: $1"; usage; exit 2;;
  esac
done

if [[ -z "$SAMPLE" || -z "$BAM" || -z "$GTF" || -z "$FASTA" ]]; then
  usage
  exit 2
fi

if [[ "$SKIP_FILTER" != "TRUE" && "$SKIP_FILTER" != "FALSE" ]]; then
  echo "[ERROR] --skip-filter must be TRUE or FALSE" >&2
  exit 2
fi

if [[ "$SKIP_ORFQUANT" != "TRUE" && "$SKIP_ORFQUANT" != "FALSE" ]]; then
  echo "[ERROR] --skip-orfquant must be TRUE or FALSE" >&2
  exit 2
fi

if [[ "$RESUME" != "TRUE" && "$RESUME" != "FALSE" ]]; then
  echo "[ERROR] --resume must be TRUE or FALSE" >&2
  exit 2
fi

if [[ "$SKIP_FILTER" == "FALSE" && -z "$FAI" ]]; then
  echo "[ERROR] --fai is required unless --skip-filter TRUE" >&2
  usage
  exit 2
fi

if [[ -n "$ORFQUANT_PKG" ]]; then
  if [[ ! -f "$ORFQUANT_PKG" ]]; then
    echo "[ERROR] --orfquant-pkg not found: $ORFQUANT_PKG" >&2
    exit 2
  fi
  ORFQUANT_PKG="$(python3 -c 'import os,sys; print(os.path.abspath(sys.argv[1]))' "$ORFQUANT_PKG")"
fi

if [[ -n "$ORFQUANT_SIF" ]]; then
  if [[ ! -f "$ORFQUANT_SIF" ]]; then
    echo "[ERROR] --orfquant-sif not found: $ORFQUANT_SIF" >&2
    exit 2
  fi
  ORFQUANT_SIF="$(python3 -c 'import os,sys; print(os.path.abspath(sys.argv[1]))' "$ORFQUANT_SIF")"
fi

if [[ -n "$ORFQUANT_RANNOT" ]]; then
  if [[ ! -f "$ORFQUANT_RANNOT" ]]; then
    echo "[ERROR] --orfquant-annotation not found: $ORFQUANT_RANNOT" >&2
    exit 2
  fi
  ORFQUANT_RANNOT="$(python3 -c 'import os,sys; print(os.path.abspath(sys.argv[1]))' "$ORFQUANT_RANNOT")"
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

mkdir -p "$OUTDIR"
OUTDIR="$(cd "$OUTDIR" && pwd)"

OUT_FILTER="$OUTDIR/00_sorf_filter"
OUT_ANNOT="$OUTDIR/01_riboseqc_annot"
OUT_ANALYSIS="$OUTDIR/02_riboseqc_analysis"
OUT_ORFQUANT="$OUTDIR/03_orfquant"

mkdir -p "$OUT_FILTER" "$OUT_ANNOT" "$OUT_ANALYSIS" "$OUT_ORFQUANT"

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

if [[ -n "$RANNOT" ]]; then
  echo "[INFO] Using provided annotation: $RANNOT"
fi

BAM_FOR_DOWNSTREAM="$BAM"
if [[ "$SKIP_FILTER" == "FALSE" ]]; then
  BAM_FOR_DOWNSTREAM="$OUT_FILTER/${SAMPLE}.sorf.filtered.bam"
  if [[ "$RESUME" == "TRUE" && -f "$BAM_FOR_DOWNSTREAM" ]]; then
    echo "[INFO] Step 1/4: BAM filtering skipped (resume): $BAM_FOR_DOWNSTREAM"
  else
    echo "[INFO] Step 1/4: BAM filtering (sorf_bam_filter)"
    bash "$SCRIPT_DIR/01_sorf_bam_filter.sh" \
      --sample "$SAMPLE" \
      --bam "$BAM" \
      --fai "$FAI" \
      --unique-mode "$UNIQUE_MODE" \
      --mapq "$MAPQ" \
      --len-min "$LEN_MIN" \
      --len-max "$LEN_MAX" \
      --exclude-regex "$EXCLUDE_REGEX" \
      --cpus "$CPUS" \
      --outdir "$OUT_FILTER"

    if [[ ! -f "$BAM_FOR_DOWNSTREAM" ]]; then
      echo "[ERROR] Filtering step did not produce: $BAM_FOR_DOWNSTREAM" >&2
      exit 2
    fi
  fi
else
  echo "[INFO] Step 1/4: Skipping filtering; using provided BAM as filtered input"
fi

if [[ -z "$RANNOT" ]]; then
  if [[ "$RESUME" == "TRUE" ]] && RANNOT="$(resolve_rannot "$OUT_ANNOT" "$GTF")"; then
    echo "[INFO] Step 2/4: RiboseQC prepareannotation skipped (resume): $RANNOT"
  else
    echo "[INFO] Step 2/4: RiboseQC prepareannotation"
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
  fi
fi

FOR_ORFQUANT="$OUT_ANALYSIS/${SAMPLE}_for_ORFquant"
if [[ "$RESUME" == "TRUE" && -f "$FOR_ORFQUANT" ]]; then
  echo "[INFO] Step 3/4: RiboseQC analysis skipped (resume): $FOR_ORFQUANT"
else
  echo "[INFO] Step 3/4: RiboseQC analysis"
  bash "$SCRIPT_DIR/03_riboseqc_analysis.sh" \
    --sample "$SAMPLE" \
    --bam "$BAM_FOR_DOWNSTREAM" \
    --annotation "$RANNOT" \
    --fasta "$FASTA" \
    --outdir "$OUT_ANALYSIS" \
    --cpus "$CPUS" \
    --fast-mode "$FAST_MODE"
fi

if [[ ! -f "$FOR_ORFQUANT" ]]; then
  echo "[ERROR] Missing RiboseQC output for ORFquant: $FOR_ORFQUANT" >&2
  echo "        Check RiboseQC outputs in: $OUT_ANALYSIS" >&2
  exit 2
fi

if [[ "$SKIP_ORFQUANT" == "TRUE" ]]; thenscm-history-item:%5Cmnt%5Cc%5CUsers%5Crenzhe%5COneDrive%20-%20BGI%20Hong%20Kong%20Tech%20Co.%2C%20Limited%5CPolyu%5C2025s1%5CRiboSeq%5Cnextflow%5Criboseq?%7B%22repositoryId%22%3A%22scm0%22%2C%22historyItemId%22%3A%2250007a63e64914a284c8f22a640b6ef2af8ad312%22%2C%22historyItemParentId%22%3A%22c947a0ca3bbffab99730519e4825c819e6d0a1e9%22%2C%22historyItemDisplayId%22%3A%2250007a6%22%7D
  echo "[INFO] Step 4/4: Skipping ORFquant (--skip-orfquant TRUE)"
  echo "[OK] Done. Outputs:"
  echo "  - Filtered BAM:         $BAM_FOR_DOWNSTREAM"
  echo "  - RiboseQC annotation: $OUT_ANNOT"
  echo "  - RiboseQC analysis:   $OUT_ANALYSIS"
  exit 0
fi

echo "[INFO] Step 4/4: ORFquant"
ORFQ_ANNOT_DIR="$OUT_ORFQUANT/orfquant_annot"
mkdir -p "$ORFQ_ANNOT_DIR"

# IMPORTANT (from the referenced Nextflow example): ORFquant requires its own Rannot.
# RiboseQC *_Rannot is not guaranteed to work for ORFquant.
if [[ -z "$ORFQUANT_RANNOT" ]]; then
  if [[ "$RESUME" == "TRUE" ]] && ORFQUANT_RANNOT="$(resolve_rannot "$ORFQ_ANNOT_DIR" "$GTF")"; then
    echo "[INFO] Step 4/4: ORFquant annotation skipped (resume): $ORFQUANT_RANNOT"
  else
    echo "[INFO] Step 4/4: Preparing ORFquant annotation (ORFquant-specific *_Rannot)"
    ORFQ_ANN_ARGS=(--gtf "$GTF" --fasta "$FASTA" --outdir "$ORFQ_ANNOT_DIR")
    if [[ -n "$ORFQUANT_SIF" ]]; then
      ORFQ_ANN_ARGS+=(--container "$ORFQUANT_SIF")
    fi
    bash "$SCRIPT_DIR/13_orfquant_prepareannotation.sh" "${ORFQ_ANN_ARGS[@]}"
    if ! ORFQUANT_RANNOT="$(resolve_rannot "$ORFQ_ANNOT_DIR" "$GTF")"; then
      echo "[ERROR] Cannot uniquely determine ORFquant *_Rannot in: $ORFQ_ANNOT_DIR" >&2
      find "$ORFQ_ANNOT_DIR" -maxdepth 1 -type f -name "*_Rannot" -print >&2 || true
      exit 2
    fi
  fi
fi

ORFQUANT_ARGS=()
if [[ -n "$ORFQUANT_PKG" ]]; then
  ORFQUANT_ARGS+=(--orfquant-pkg "$ORFQUANT_PKG")
fi
if [[ -n "$ORFQUANT_SIF" ]]; then
  ORFQUANT_ARGS+=(--container "$ORFQUANT_SIF")
fi

if [[ "$RESUME" == "TRUE" ]] && find "$OUT_ORFQUANT" -maxdepth 1 -name "${SAMPLE}_final_ORFquant_results*" -print -quit | grep -q .; then
  echo "[INFO] Step 4/4: ORFquant skipped (resume): ${SAMPLE}_final_ORFquant_results*"
else
  bash "$SCRIPT_DIR/04_orfquant_run.sh" \
    --sample "$SAMPLE" \
    --for-orfquant "$FOR_ORFQUANT" \
    --annotation "$ORFQUANT_RANNOT" \
    --fasta "$FASTA" \
    --cpus "$CPUS" \
    --outdir "$OUT_ORFQUANT" \
    "${ORFQUANT_ARGS[@]}"
fi

echo "[OK] Done. Outputs:"
echo "  - Filtered BAM:         $BAM_FOR_DOWNSTREAM"
echo "  - RiboseQC annotation: $OUT_ANNOT"
echo "  - RiboseQC analysis:   $OUT_ANALYSIS"
echo "  - ORFquant:            $OUT_ORFQUANT"
