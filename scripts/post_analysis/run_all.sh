#!/bin/bash
# ─── ORF Post-Analysis Pipeline (Two-Stage Filtering) ─────────────────
# Usage: bash run_all.sh /path/to/project_config.yaml
#
#   Step 1: Preliminary analysis → prelim_orfs_for_psite.bed (R/Quarto)
#   Step 2: P-site purity (Python) — runs on PRELIMINARY ORFs only
#   Step 3: P-site filtering + final ORF selection (R/Quarto)
#   Step 4: Generate ggRibo coverage plots (R/Quarto)
#
# Requirements:
#   - Python 3 with numpy
#   - R with tidyverse, ggRibo, yaml, R6, furrr, txdbmaker, GenomicFeatures
#   - Quarto CLI
#
# Output:
#   {output_dir}/01_prelim_analysis.html
#   {output_dir}/prelim_orfs_for_psite.bed
#   {output_dir}/psite_purity.tsv
#   {output_dir}/02_psite_filtering.html
#   {output_dir}/final_orfs_for_ggribo.tsv
#   {output_dir}/ggribo_plots/

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

# Read output_dir and result_dir from YAML
OUTPUT_DIR=$(grep -E '^\s+output_dir:' "$CONFIG" | head -1 | sed 's/.*output_dir:\s*"\(.*\)"/\1/' | sed 's/.*output_dir:\s*//' | xargs)
RESULT_DIR=$(grep -E '^\s+result_dir:' "$CONFIG" | head -1 | sed 's/.*result_dir:\s*"\(.*\)"/\1/' | sed 's/.*result_dir:\s*//' | xargs)
# Expand {result_dir}
OUTPUT_DIR="${OUTPUT_DIR//\{result_dir\}/$RESULT_DIR}"
mkdir -p "$OUTPUT_DIR"
mkdir -p "$OUTPUT_DIR/logs"

# Read riboseqc dir for P-site computation
RIBOSEQC_DIR=$(grep -E '^\s+riboseqc_dir:' "$CONFIG" | head -1 | sed 's/.*riboseqc_dir:\s*"\(.*\)"/\1/' | sed 's/.*riboseqc_dir:\s*//' | xargs)
RIBOSEQC_DIR="${RIBOSEQC_DIR//\{result_dir\}/$RESULT_DIR}"

echo "============================================"
echo "ORF Post-Analysis Pipeline (Two-Stage)"
echo "============================================"
echo "Config:      $CONFIG"
echo "Output dir:  $OUTPUT_DIR"
echo ""

# ─── Step 1: Preliminary Analysis ──────────────────────────────────────
echo "=== Step 1: Preliminary Expression-Based Analysis ==="
quarto render "$SCRIPT_DIR/01_prelim_analysis.qmd" \
  -P config:"$CONFIG" \
  --execute \
  --output-dir "$OUTPUT_DIR" \
  2>&1 | tee "$OUTPUT_DIR/logs/step1_prelim.log"
echo "  → $OUTPUT_DIR/01_prelim_analysis.html"
echo ""

# ─── Step 2: P-site Purity (on PRELIMINARY ORFs) ──────────────────────
echo "=== Step 2: Real P-site Purity (bedgraph backtracking) ==="
PURITY_OUT="$OUTPUT_DIR/psite_purity.tsv"
PRELIM_BED="$OUTPUT_DIR/prelim_orfs_for_psite.bed"

if [ ! -f "$PRELIM_BED" ]; then
  echo "ERROR: prelim_orfs_for_psite.bed not found — did Step 1 succeed?"
  exit 1
fi

N_PRELIM=$(wc -l < "$PRELIM_BED")
echo "  ORFs to scan: $N_PRELIM (preliminary filtered)"
echo "  RiboseQC dir: $RIBOSEQC_DIR"

if [ -f "$PURITY_OUT" ]; then
  echo "  psite_purity.tsv already exists — skipping (delete to recompute)"
else
  python3 "$SCRIPT_DIR/compute_psite_purity.py" \
    --bed "$PRELIM_BED" \
    --riboseqc-dir "$RIBOSEQC_DIR" \
    --output "$PURITY_OUT" \
    --workers 8 \
    2>&1 | tee "$OUTPUT_DIR/logs/step2_purity.log"
  echo "  → $PURITY_OUT"
fi
echo ""

# ─── Step 3: P-site Filtering + Final ORF Selection ────────────────────
echo "=== Step 3: Real P-site Filtering ==="
quarto render "$SCRIPT_DIR/02_psite_filtering.qmd" \
  -P config:"$CONFIG" \
  --execute \
  --output-dir "$OUTPUT_DIR" \
  2>&1 | tee "$OUTPUT_DIR/logs/step3_psite_filter.log"
echo "  → $OUTPUT_DIR/02_psite_filtering.html"
echo ""

# ─── Step 4: ggRibo Coverage Plots ─────────────────────────────────────
ORF_LIST="$OUTPUT_DIR/final_orfs_for_ggribo.tsv"
if [ -f "$ORF_LIST" ]; then
  N_ORFS=$(tail -n +2 "$ORF_LIST" | wc -l)
  echo "=== Step 4: ggRibo Coverage Plots ($N_ORFS ORFs) ==="
  quarto render "$SCRIPT_DIR/03_generate_ggribo.qmd" \
    -P config:"$CONFIG" \
    --execute \
    --output-dir "$OUTPUT_DIR" \
    2>&1 | tee "$OUTPUT_DIR/logs/step4_ggribo.log"
  echo "  → $OUTPUT_DIR/ggribo_plots/"
else
  echo ""
  echo "=== Step 4: SKIPPED ($ORF_LIST not found) ==="
fi

echo ""
echo "============================================"
echo "All done. Output in: $OUTPUT_DIR"
echo "  Step 1: 01_prelim_analysis.html"
echo "  Step 2: psite_purity.tsv"
echo "  Step 3: 02_psite_filtering.html"
echo "  Step 4: ggribo_plots/"
echo "============================================"
