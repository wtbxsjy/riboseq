#!/bin/bash
# ─── Pipeline 各阶段输出体检 ─────────────────────────────────────────────
# Usage:
#   check_pipeline_outputs.sh <result_dir> [--samples-file samplesheet.csv] [--expected-orfs N]
#
# 检查 result/ 下各阶段关键文件的存在性、行数与一致性，输出 PASS/WARN/FAIL。
# 退出码：0 全部通过；1 有警告；2 有失败。
#
# 基准（rice 实测）：23 样本 → 369,939 unified ORFs；
#   summary 行数 == metadata 行数；gencode biotype 分布 lncRNA 不应 >90%。
#
# Examples:
#   check_pipeline_outputs.sh run/rice/result
#   check_pipeline_outputs.sh run/rice/result --samples-file run/rice/scripts/samplesheet.csv --expected-orfs 370000

set -euo pipefail

[ $# -ge 1 ] || { echo "Usage: check_pipeline_outputs.sh <result_dir> [--samples-file F] [--expected-orfs N]"; exit 2; }

RESULT_DIR="$1"; shift
SAMPLES_FILE=""
EXPECTED_ORFS=""

while [ $# -gt 0 ]; do
  case "$1" in
    --samples-file)  SAMPLES_FILE="$2"; shift 2 ;;
    --expected-orfs) EXPECTED_ORFS="$2"; shift 2 ;;
    *) echo "ERROR: unknown option: $1"; exit 2 ;;
  esac
done

[ -d "$RESULT_DIR" ] || { echo "FAIL: result dir not found: $RESULT_DIR"; exit 2; }

warns=0; fails=0
ok()   { echo "PASS: $1"; }
warn() { echo "WARN: $1"; warns=$((warns+1)); }
fail() { echo "FAIL: $1"; fails=$((fails+1)); }

echo "=== checking $RESULT_DIR ==="

# ─── 样本数 ──────────────────────────────────────────────────────────────
N_SAMPLES=""
if [ -n "$SAMPLES_FILE" ] && [ -f "$SAMPLES_FILE" ]; then
  N_SAMPLES=$(tail -n +2 "$SAMPLES_FILE" | grep -c . || true)
fi

# ─── 1. 阶段目录 ─────────────────────────────────────────────────────────
for d in alignment genome riboseqc orf_unification orf_classification orf expression multiqc; do
  [ -d "$RESULT_DIR/$d" ] || fail "missing dir: $d"
done
[ -d "$RESULT_DIR/orf_predictions" ] || warn "missing dir: orf_predictions (no ORF prediction?)"

# ─── 2. riboseqc bedgraphs ───────────────────────────────────────────────
if [ -d "$RESULT_DIR/riboseqc" ]; then
  N_PS=$(find "$RESULT_DIR/riboseqc" -maxdepth 1 -name "*_P_sites_plus.bedgraph" | wc -l)
  [ "$N_PS" -gt 0 ] || fail "riboseqc: no _P_sites_plus.bedgraph found"
  if [ -n "$N_SAMPLES" ]; then
    [ "$N_PS" -eq "$N_SAMPLES" ] || warn "riboseqc: $N_PS P_sites bedgraphs vs $N_SAMPLES samples"
  fi
  ok "riboseqc: $N_PS samples with P_sites bedgraph"
fi

# ─── 3. 各 ORF 预测工具 per-sample 输出 ──────────────────────────────────
for tool in ribotish ribotricer orfquant price ribocode; do
  dir="$RESULT_DIR/orf_predictions/$tool"
  [ -d "$dir" ] || { echo "SKIP: orf_predictions/$tool (not run)"; continue; }
  n_files=$(find "$dir" -type f ! -name "*.log" | wc -l)
  n_empty=$(find "$dir" -type f ! -name "*.log" -size 0 | wc -l)
  if [ "$n_files" -eq 0 ]; then
    fail "$tool: dir exists but no output files"
  elif [ -n "$N_SAMPLES" ] && [ "$n_files" -lt "$N_SAMPLES" ]; then
    warn "$tool: $n_files outputs vs $N_SAMPLES samples"
  else
    ok "$tool: $n_files outputs"
  fi
  [ "$n_empty" -eq 0 ] || warn "$tool: $n_empty empty output files"
done

# ─── 4. unification ──────────────────────────────────────────────────────
META="$RESULT_DIR/orf_unification/unified_orfs.metadata.tsv"
if [ -f "$META" ]; then
  N_ORFS=$(tail -n +2 "$META" | wc -l)
  ok "unified_orfs.metadata.tsv: $N_ORFS ORFs"
  if [ -n "$EXPECTED_ORFS" ]; then
    LO=$((EXPECTED_ORFS * 50 / 100)); HI=$((EXPECTED_ORFS * 150 / 100))
    if [ "$N_ORFS" -lt "$LO" ] || [ "$N_ORFS" -gt "$HI" ]; then
      warn "ORF count $N_ORFS outside ±50% of expected $EXPECTED_ORFS"
    fi
  fi
else
  fail "unified_orfs.metadata.tsv missing"
fi

if [ -f "$RESULT_DIR/orf_unification/unified_orfs.stats.txt" ]; then
  ok "stats: $(tr '\n' ' ' < "$RESULT_DIR/orf_unification/unified_orfs.stats.txt" | cut -c1-160)"
else
  warn "unified_orfs.stats.txt missing"
fi

# ─── 5. expression 一致性 ────────────────────────────────────────────────
SUMMARY=$(ls "$RESULT_DIR/orf_unification/unified_orfs_expression_summary.tsv" \
              "$RESULT_DIR/expression/expression_quant_expression_summary.tsv" 2>/dev/null | head -1)
if [ -n "${SUMMARY:-}" ] && [ -f "${SUMMARY}" ]; then
  N_SUM=$(tail -n +2 "$SUMMARY" | wc -l)
  if [ -n "${N_ORFS:-}" ] && [ "$N_SUM" -eq "$N_ORFS" ]; then
    ok "expression summary rows == metadata rows ($N_SUM)"
  else
    warn "expression summary rows ($N_SUM) != metadata rows (${N_ORFS:-?})"
  fi
else
  warn "expression summary not found (expression_quant not run?)"
fi

# ─── 6. ORF_QC confidence ────────────────────────────────────────────────
CONF="$RESULT_DIR/orf/unified_orfs_orf_confidence.tsv"
if [ -f "$CONF" ]; then
  N_CONF=$(tail -n +2 "$CONF" | wc -l)
  if [ -n "${N_ORFS:-}" ] && [ "$N_CONF" -ne "$N_ORFS" ]; then
    warn "orf confidence rows ($N_CONF) != metadata rows ($N_ORFS)"
  else
    ok "orf confidence rows: $N_CONF"
  fi
else
  warn "unified_orfs_orf_confidence.tsv missing (ORF_QC skipped or errorStrategy ignored)"
fi

# ─── 7. GENCODE 分类 ─────────────────────────────────────────────────────
GENOUT="$RESULT_DIR/orf_classification/gencode/gencode_results.orfs.out"
CAT_GEN="cat"
if [ ! -f "$GENOUT" ] && [ -f "$GENOUT.gz" ]; then
  GENOUT="$GENOUT.gz"
  CAT_GEN="zcat"
fi
if [ -f "$GENOUT" ]; then
  N_GEN=$($CAT_GEN "$GENOUT" | tail -n +2 | wc -l)
  LNCRNA=$($CAT_GEN "$GENOUT" | cut -f7 | grep -c "^lncRNA$" || true)
  ok "gencode_results.orfs.out(.gz): $N_GEN classified"
  if [ "$N_GEN" -gt 0 ]; then
    PCT=$((LNCRNA * 100 / N_GEN))
    [ "$PCT" -le 90 ] || warn "lncRNA fraction ${PCT}% > 90% — possible protein header mismatch (transcript_id vs protein_id)"
  fi
else
  warn "gencode_results.orfs.out(.gz) missing (GENCODE classification skipped?)"
fi

# ─── 8. ORF-type 分类 ────────────────────────────────────────────────────
ORFTYPE="$RESULT_DIR/orf_classification/orf_type/orftype_classification.tsv"
[ -f "$ORFTYPE" ] && ok "orftype_classification.tsv: $(tail -n +2 "$ORFTYPE" | wc -l) rows" \
                  || warn "orftype_classification.tsv missing (orf_classify_mode != orf_type?)"

# ─── 汇总 ────────────────────────────────────────────────────────────────
echo ""
echo "=== summary: $fails fail, $warns warn ==="
[ "$fails" -eq 0 ] || exit 2
[ "$warns" -eq 0 ] || exit 1
exit 0
