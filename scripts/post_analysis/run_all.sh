#!/bin/bash
# ─── ORF Post-Analysis Pipeline ──────────────────────────────────────
# Usage: bash run_all.sh /path/to/project_config.yaml
#   Step 1: Compute P-site purity (Python)
#   Step 2: Filtering analysis + ORF selection (R/Quarto)
#   Step 3: Generate ggRibo coverage plots (R/Quarto)
#
# Requirements:
#   - Python 3 with numpy, pyarrow or pandas
#   - R with tidyverse, ggRibo, yaml, R6, furrr
#   - Quarto CLI
#
# Output:
#   {output_dir}/psite_purity.tsv
#   {output_dir}/01_filtering_analysis.html
#   {output_dir}/final_orfs_for_ggribo.tsv
#   {output_dir}/ggribo_plots/index.html

set -euo pipefail

CONFIG="${1:-}"
if [ -z "$CONFIG" ]; then
  echo "Usage: bash run_all.sh /path/to/project_config.yaml"
  exit 1
fi

if [ ! -f "$CONFIG" ]; then
  echo "ERROR: Config file not found: $CONFIG"
  exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

# Read output_dir from YAML (simple grep, works for unquoted paths)
OUTPUT_DIR=$(grep -E '^\s+output_dir:' "$CONFIG" | head -1 | sed 's/.*output_dir:\s*"\(.*\)"/\1/' | sed 's/.*output_dir:\s*//' | xargs)
RESULT_DIR=$(grep -E '^\s+result_dir:' "$CONFIG" | head -1 | sed 's/.*result_dir:\s*"\(.*\)"/\1/' | sed 's/.*result_dir:\s*//' | xargs)
# Expand {result_dir}
OUTPUT_DIR="${OUTPUT_DIR//\{result_dir\}/$RESULT_DIR}"
mkdir -p "$OUTPUT_DIR"
mkdir -p "$OUTPUT_DIR/logs"

echo "============================================"
echo "ORF Post-Analysis Pipeline"
echo "============================================"
echo "Config:      $CONFIG"
echo "Output dir:  $OUTPUT_DIR"
echo ""

# ─── Step 1: P-site Purity ────────────────────────────────────────────
echo "=== Step 1: P-site Purity ==="
PURITY_OUT="$OUTPUT_DIR/psite_purity.tsv"
BED=$(grep -E '^\s+unified_bed:' "$CONFIG" | head -1 | sed 's/.*unified_bed:\s*"\(.*\)"/\1/' | sed 's/.*unified_bed:\s*//' | xargs)
BED="${BED//\{result_dir\}/$RESULT_DIR}"
RIBOSEQC_DIR=$(grep -E '^\s+riboseqc_dir:' "$CONFIG" | head -1 | sed 's/.*riboseqc_dir:\s*"\(.*\)"/\1/' | sed 's/.*riboseqc_dir:\s*//' | xargs)
RIBOSEQC_DIR="${RIBOSEQC_DIR//\{result_dir\}/$RESULT_DIR}"

if [ -f "$PURITY_OUT" ]; then
  echo "  psite_purity.tsv already exists — skipping"
else
  echo "  Computing P-site purity for: $BED"
  echo "  RiboseQC directory: $RIBOSEQC_DIR"
  python3 "$SCRIPT_DIR/compute_psite_purity.py" \
    --bed "$BED" \
    --riboseqc-dir "$RIBOSEQC_DIR" \
    --output "$PURITY_OUT" \
    --workers 8 \
    2>&1 | tee "$OUTPUT_DIR/logs/step1_purity.log"
  echo "  → $PURITY_OUT"
fi

# ─── Step 2: Filtering Analysis ───────────────────────────────────────
echo ""
echo "=== Step 2: ORF Filtering Analysis ==="
quarto render "$SCRIPT_DIR/01_filtering_analysis.qmd" \
  -P config:"$CONFIG" \
  --execute \
  --output-dir "$OUTPUT_DIR" \
  2>&1 | tee "$OUTPUT_DIR/logs/step2_filtering.log"

echo "  → $OUTPUT_DIR/01_filtering_analysis.html"

# ─── Step 3: ggRibo Plots ─────────────────────────────────────────────
ORF_LIST="$OUTPUT_DIR/final_orfs_for_ggribo.tsv"
if [ -f "$ORF_LIST" ]; then
  N_ORFS=$(tail -n +2 "$ORF_LIST" | wc -l)
  echo ""
  echo "=== Step 3: ggRibo Coverage Plots ($N_ORFS ORFs) ==="
  quarto render "$SCRIPT_DIR/02_generate_ggribo.qmd" \
    -P config:"$CONFIG" \
    --execute \
    --output-dir "$OUTPUT_DIR" \
    2>&1 | tee "$OUTPUT_DIR/logs/step3_ggribo.log"
  echo "  → $OUTPUT_DIR/ggribo_plots/"
else
  echo ""
  echo "=== Step 3: SKIPPED ($ORF_LIST not found) ==="
fi

echo ""
echo "============================================"
echo "All done. Output in: $OUTPUT_DIR"
echo "============================================"
