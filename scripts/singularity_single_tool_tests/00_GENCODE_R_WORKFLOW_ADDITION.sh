#!/usr/bin/env bash
#
# 补充内容：GENCODE 注释和 R 分析工作流
#
# 此文件包含应该添加到 00_run_order_example.sh 的新步骤（步骤 10-12）
#
# 使用方法：
# 1. 手动将下面的内容添加到 00_run_order_example.sh 的第 10 步之前
# 2. 或者在完成步骤 1-9 后，单独运行此脚本
#
###############################################################################

###############################################################################
# 10) Step I: GENCODE 格式转换（多工具统一注释前的准备）
###############################################################################

# 10.1 准备 Ensembl 注释目录（一次性，适用于所有样本）
# 注意：此步骤需要下载大量数据，建议提前执行

OUT_13="${BASE_OUT}/out_13_ensembl_m35"
ENSEMBL_DIR="${ENSEMBL_DIR:-${OUT_13}}"

echo ""
echo "[INFO] 准备 GENCODE M35 (小鼠) Ensembl 注释目录..."
echo "- 如果已准备好，请设置环境变量：ENSEMBL_DIR=/path/to/Ens_GENCODE_M35"
echo "- 否则将下载到：${OUT_13}"

if [[ ! -d "${ENSEMBL_DIR}/PROTEOME_FASTA" ]]; then
  run "./prepare_gencode_m35_ensembl.sh '${OUT_13}'"
else
  echo "[INFO] Ensembl 目录已存在，跳过下载"
fi

# 10.2 转换 Ribo-TISH 预测结果为 GENCODE 格式

OUT_14="${BASE_OUT}/out_14_ribotish_gencode"

RIBOTISH_PRED="${RIBOTISH_PRED:-${OUT_07}/${SAMPLE}_pred.txt}"

run "./14_ribotish_to_gencode.sh \
  --sample '${SAMPLE}' \
  --predict '${RIBOTISH_PRED}' \
  --fasta '${FASTA}' \
  --min-length 16 \
  --outdir '${OUT_14}'"

# 10.3 转换 Ribotricer 预测结果为 GENCODE 格式

OUT_15="${BASE_OUT}/out_15_ribotricer_gencode"

RIBOTRICER_TSV="${RIBOTRICER_TSV:-${OUT_09}/${SAMPLE}_translating_ORFs.tsv}"

run "./15_ribotricer_to_gencode.sh \
  --sample '${SAMPLE}' \
  --tsv '${RIBOTRICER_TSV}' \
  --fasta '${FASTA}' \
  --min-length 16 \
  --min-phase-score 0.5 \
  --outdir '${OUT_15}'"

# 10.4 ORFquant 格式转换（如需要）
# 注意：ORFquant 的输出格式取决于具体版本，可能需要自定义转换脚本
# 此处作为占位符，请根据实际输出格式调整

OUT_16="${BASE_OUT}/out_16_orfquant_gencode"

echo ""
echo "[INFO] ORFquant 格式转换..."
echo "- ORFquant 输出格式因版本而异，请检查 ${OUT_05} 中的文件"
echo "- 如果已是 GENCODE 兼容格式，可跳过此步骤"
echo "- 否则需要自定义转换脚本（参考 14/15 脚本）"

mkdir -p "${OUT_16}"

###############################################################################
# 11) Step J: GENCODE ORF Mapper（统一注释 + 去冗余）
###############################################################################

# 注意：此步骤是多工具比较的关键！
# - 合并所有工具的 ORF 预测
# - 统一映射到 GENCODE/Ensembl 注释
# - 去除冗余（相似序列折叠）
# - 生成跨工具的检测矩阵

OUT_17="${BASE_OUT}/out_17_gencode_unified"

# 11.1 合并所有工具的 GENCODE 格式文件

MERGED_FA="${OUT_17}/merged_all_tools.fa"
MERGED_BED="${OUT_17}/merged_all_tools.bed"

mkdir -p "${OUT_17}"

echo ""
echo "[INFO] 合并所有工具的 ORF 预测..."

# 合并 FASTA（如果文件存在）
run "cat ${OUT_14}/*.gencode.fa \
     ${OUT_15}/*.gencode.fa \
     ${OUT_16}/*.gencode.fa \
     2>/dev/null > '${MERGED_FA}' || touch '${MERGED_FA}'"

# 合并 BED（如果文件存在）
run "cat ${OUT_14}/*.gencode.bed \
     ${OUT_15}/*.gencode.bed \
     ${OUT_16}/*.gencode.bed \
     2>/dev/null > '${MERGED_BED}' || touch '${MERGED_BED}'"

# 11.2 运行 GENCODE ORF mapper

PROJECT_ID="${PROJECT_ID:-Mouse_AllTools_Comparison}"

# 重要参数说明：
# --collapse-threshold 0.78: 宽松阈值，适合多工具比较（保留工具特异性）
# --collapse-method longest_string: 基于序列相似度去冗余
# --min-length 16: 最小 ORF 长度（氨基酸）

run "./16_gencode_orf_mapper.sh \
  --project '${PROJECT_ID}' \
  --fasta '${MERGED_FA}' \
  --bed '${MERGED_BED}' \
  --ensembl-dir '${ENSEMBL_DIR}' \
  --min-length 16 \
  --collapse-threshold 0.78 \
  --collapse-method longest_string \
  --outdir '${OUT_17}'"

# 关键输出文件：
UNIFIED_GTF="${OUT_17}/${PROJECT_ID}.orfs.gtf"         # 用于定量 ⭐
UNIFIED_OUT="${OUT_17}/${PROJECT_ID}.orfs.out"         # 用于工具比较 ⭐

###############################################################################
# 12) Step K: R 语言分析（定量 + 工具比较）
###############################################################################

# 12.1 P-site 定量（使用 RiboseQC 的 bedgraph 文件）

OUT_18="${BASE_OUT}/out_18_quantification"

# RiboseQC 输出的 P-site bedgraph 文件位置
# 格式：*_P_sites_unique_plus.bedgraph / *_P_sites_unique_minus.bedgraph
BEDGRAPH_DIR="${BEDGRAPH_DIR:-${OUT_03}}"  # 使用 postfilter QC 的结果

echo ""
echo "[INFO] 运行 R 脚本：P-site 定量..."
echo "- 输入 GTF: ${UNIFIED_GTF}"
echo "- Bedgraph 目录: ${BEDGRAPH_DIR}"
echo "- 注意：需要安装 R 包（见 scripts/R/README_R_WORKFLOW.md）"

run "Rscript ../R/quantify_orfs_from_psites.R \
  --gtf '${UNIFIED_GTF}' \
  --bedgraph-dir '${BEDGRAPH_DIR}' \
  --sample-pattern '(.+)_P_sites_unique_(plus|minus)\\.bedgraph$' \
  --outdir '${OUT_18}' \
  --min-count 5 \
  --threads '${CPUS}'"

# 关键输出：
# - orf_counts_matrix.csv: 用于 DESeq2 差异分析 ⭐
# - orf_counts_tpm.tsv: 用于可视化
# - orf_summary_stats.tsv: ORF 质量过滤参考

# 12.2 工具比较分析

OUT_19="${BASE_OUT}/out_19_tool_comparison"

# 工具列表（逗号分隔，需匹配 .orfs.out 文件中的列名）
TOOLS="${TOOLS:-ORFquant,RiboTISH,Ribotricer}"

echo ""
echo "[INFO] 运行 R 脚本：工具比较分析..."
echo "- 输入文件: ${UNIFIED_OUT}"
echo "- 工具列表: ${TOOLS}"

run "Rscript ../R/analyze_tool_comparison.R \
  --input '${UNIFIED_OUT}' \
  --tools '${TOOLS}' \
  --outdir '${OUT_19}' \
  --min-samples 10 \
  --min-tools 2"

# 关键输出：
# - high_confidence_orfs.tsv: 高置信度 ORF（≥2工具+≥10样本）⭐
# - upset_plot.pdf: 工具交集可视化
# - *_specific_orfs.tsv: 各工具特异性 ORF

###############################################################################
# 13) 完成提示（更新版）
###############################################################################

echo ""
echo "=========================================="
echo "