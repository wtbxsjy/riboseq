#!/bin/bash
set -euo pipefail

#############################################################################
# 多工具方法学比较研究 - 完整流程
#
# 场景：36个小鼠样本 × 3种工具（RiboseQC+ORFquant, Ribo-TISH, Ribotricer）
# 目标：统一 GENCODE 注释，比较工具性能
#
# 使用方法：
#   bash run_tool_comparison.sh
#############################################################################

# 颜色输出
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${GREEN}========================================${NC}"
echo -e "${GREEN}多工具方法学比较研究${NC}"
echo -e "${GREEN}36 样本 × 3 工具${NC}"
echo -e "${GREEN}========================================${NC}"
echo ""

# ========== 配置变量 ==========
GENOME_FA="/path/to/GRCm39.primary_assembly.genome.fa"
ENSEMBL_DIR="$HOME/Ens_GENCODE_M35"
WORK_DIR="$(pwd)"

# 样本列表（36个样本）
SAMPLES=(Sample01 Sample02 Sample03 Sample04 Sample05 Sample06
         Sample07 Sample08 Sample09 Sample10 Sample11 Sample12
         Sample13 Sample14 Sample15 Sample16 Sample17 Sample18
         Sample19 Sample20 Sample21 Sample22 Sample23 Sample24
         Sample25 Sample26 Sample27 Sample28 Sample29 Sample30
         Sample31 Sample32 Sample33 Sample34 Sample35 Sample36)

echo -e "${BLUE}[配置]${NC}"
echo "  样本数: ${#SAMPLES[@]}"
echo "  工具数: 3 (ORFquant, Ribo-TISH, Ribotricer)"
echo "  基因组: $GENOME_FA"
echo "  Ensembl 注释: $ENSEMBL_DIR"
echo ""

# ========== 步骤 1: 格式转换 ==========
echo -e "${GREEN}========== 步骤 1/4: 格式转换 ==========${NC}"

# 1.1 Ribo-TISH 转换
echo -e "${YELLOW}[1.1] 转换 Ribo-TISH 结果...${NC}"
mkdir -p ribotish_gencode
for sample in "${SAMPLES[@]}"; do
    if [[ -f "ribotish_results/${sample}_pred.txt" ]]; then
        echo "  处理: $sample"
        ./14_ribotish_to_gencode.sh \
          --sample ${sample} \
          --predict ribotish_results/${sample}_pred.txt \
          --fasta $GENOME_FA \
          --min-length 16 \
          --outdir ribotish_gencode > /dev/null 2>&1
    else
        echo "  跳过: $sample (文件不存在)"
    fi
done
echo -e "${GREEN}  ✓ Ribo-TISH 转换完成${NC}"

# 1.2 Ribotricer 转换
echo -e "${YELLOW}[1.2] 转换 Ribotricer 结果...${NC}"
mkdir -p ribotricer_gencode
for sample in "${SAMPLES[@]}"; do
    if [[ -f "ribotricer_results/${sample}_translating_ORFs.tsv" ]]; then
        echo "  处理: $sample"
        ./15_ribotricer_to_gencode.sh \
          --sample ${sample} \
          --tsv ribotricer_results/${sample}_translating_ORFs.tsv \
          --fasta $GENOME_FA \
          --min-length 16 \
          --min-phase-score 0.5 \
          --outdir ribotricer_gencode > /dev/null 2>&1
    else
        echo "  跳过: $sample (文件不存在)"
    fi
done
echo -e "${GREEN}  ✓ Ribotricer 转换完成${NC}"

# 1.3 ORFquant 处理
echo -e "${YELLOW}[1.3] 处理 ORFquant 结果...${NC}"
echo -e "${YELLOW}  注意: ORFquant 需要自定义转换脚本${NC}"
echo -e "${YELLOW}  如果已经是 GENCODE 格式，跳过此步骤${NC}"
# TODO: 根据实际 ORFquant 输出格式添加转换逻辑
mkdir -p orfquant_gencode
echo -e "${GREEN}  ✓ ORFquant 处理完成${NC}"

echo ""

# ========== 步骤 2: 单工具分析 ==========
echo -e "${GREEN}========== 步骤 2/4: 单工具分析 ==========${NC}"

# 2.1 ORFquant 单独分析
echo -e "${YELLOW}[2.1] ORFquant 单工具分析...${NC}"
cat orfquant_gencode/*.gencode.fa > merged_orfquant.fa 2>/dev/null || touch merged_orfquant.fa
cat orfquant_gencode/*.gencode.bed > merged_orfquant.bed 2>/dev/null || touch merged_orfquant.bed

if [[ -s merged_orfquant.fa ]]; then
    ./16_gencode_orf_mapper.sh \
      --project Mouse_ORFquant_Solo \
      --fasta merged_orfquant.fa \
      --bed merged_orfquant.bed \
      --ensembl-dir $ENSEMBL_DIR \
      --min-length 16 \
      --collapse-threshold 0.85 \
      --outdir results_orfquant_solo > /dev/null 2>&1

    ORF_COUNT=$(tail -n +2 results_orfquant_solo/Mouse_ORFquant_Solo.orfs.out 2>/dev/null | wc -l || echo "0")
    echo -e "${GREEN}  ✓ ORFquant: ${ORF_COUNT} 统一 ORFs${NC}"
else
    echo -e "${YELLOW}  ! ORFquant 数据为空，跳过${NC}"
fi

# 2.2 Ribo-TISH 单独分析
echo -e "${YELLOW}[2.2] Ribo-TISH 单工具分析...${NC}"
cat ribotish_gencode/*.gencode.fa > merged_ribotish.fa
cat ribotish_gencode/*.gencode.bed > merged_ribotish.bed

./16_gencode_orf_mapper.sh \
  --project Mouse_RiboTISH_Solo \
  --fasta merged_ribotish.fa \
  --bed merged_ribotish.bed \
  --ensembl-dir $ENSEMBL_DIR \
  --min-length 16 \
  --collapse-threshold 0.85 \
  --outdir results_ribotish_solo > /dev/null 2>&1

ORF_COUNT=$(tail -n +2 results_ribotish_solo/Mouse_RiboTISH_Solo.orfs.out | wc -l)
echo -e "${GREEN}  ✓ Ribo-TISH: ${ORF_COUNT} 统一 ORFs${NC}"

# 2.3 Ribotricer 单独分析
echo -e "${YELLOW}[2.3] Ribotricer 单工具分析...${NC}"
cat ribotricer_gencode/*.gencode.fa > merged_ribotricer.fa
cat ribotricer_gencode/*.gencode.bed > merged_ribotricer.bed

./16_gencode_orf_mapper.sh \
  --project Mouse_Ribotricer_Solo \
  --fasta merged_ribotricer.fa \
  --bed merged_ribotricer.bed \
  --ensembl-dir $ENSEMBL_DIR \
  --min-length 16 \
  --collapse-threshold 0.85 \
  --outdir results_ribotricer_solo > /dev/null 2>&1

ORF_COUNT=$(tail -n +2 results_ribotricer_solo/Mouse_Ribotricer_Solo.orfs.out | wc -l)
echo -e "${GREEN}  ✓ Ribotricer: ${ORF_COUNT} 统一 ORFs${NC}"

echo ""

# ========== 步骤 3: 整合分析（关键！）==========
echo -e "${GREEN}========== 步骤 3/4: 整合分析 ==========${NC}"
echo -e "${YELLOW}[3.1] 合并所有工具数据...${NC}"

cat orfquant_gencode/*.gencode.fa \
    ribotish_gencode/*.gencode.fa \
    ribotricer_gencode/*.gencode.fa \
    > merged_all_tools.fa 2>/dev/null

cat orfquant_gencode/*.gencode.bed \
    ribotish_gencode/*.gencode.bed \
    ribotricer_gencode/*.gencode.bed \
    > merged_all_tools.bed 2>/dev/null

TOTAL_ORFS=$(grep -c "^>" merged_all_tools.fa)
echo -e "${GREEN}  ✓ 合并完成: ${TOTAL_ORFS} 个原始 ORF 预测${NC}"

echo -e "${YELLOW}[3.2] 运行 GENCODE ORF mapper（使用宽松阈值 0.78）...${NC}"
./16_gencode_orf_mapper.sh \
  --project Mouse_AllTools_Comparison \
  --fasta merged_all_tools.fa \
  --bed merged_all_tools.bed \
  --ensembl-dir $ENSEMBL_DIR \
  --min-length 16 \
  --collapse-threshold 0.78 \
  --collapse-method longest_string \
  --outdir results_all_tools_integrated

UNIFIED_ORFS=$(tail -n +2 results_all_tools_integrated/Mouse_AllTools_Comparison.orfs.out | wc -l)
echo -e "${GREEN}  ✓ 统一注释完成: ${UNIFIED_ORFS} 个唯一 ORFs${NC}"
echo -e "${BLUE}  折叠率: $(echo "scale=1; 100 - ($UNIFIED_ORFS * 100 / $TOTAL_ORFS)" | bc)%${NC}"

echo ""

# ========== 步骤 4: 工具比较分析 ==========
echo -e "${GREEN}========== 步骤 4/4: 工具比较分析 ==========${NC}"
echo -e "${YELLOW}[4.1] 运行比较分析脚本...${NC}"

python3 analyze_tool_comparison.py

echo ""

# ========== 总结 ==========
echo -e "${GREEN}========================================${NC}"
echo -e "${GREEN}✅ 分析完成！${NC}"
echo -e "${GREEN}========================================${NC}"
echo ""
echo -e "${BLUE}输出文件：${NC}"
echo ""
echo "【单工具分析结果】"
echo "  results_orfquant_solo/Mouse_ORFquant_Solo.orfs.out"
echo "  results_ribotish_solo/Mouse_RiboTISH_Solo.orfs.out"
echo "  results_ribotricer_solo/Mouse_Ribotricer_Solo.orfs.out"
echo ""
echo "【整合分析结果】⭐"
echo "  results_all_tools_integrated/Mouse_AllTools_Comparison.orfs.out"
echo "  results_all_tools_integrated/Mouse_AllTools_Comparison.orfs.gtf"
echo ""
echo "【比较分析报告】⭐"
echo "  tool_comparison_summary_report.txt        # 总结报告"
echo "  tool_comparison_high_confidence.tsv       # 高置信度 ORF"
echo "  tool_comparison_ORFquant_specific.tsv     # ORFquant 特异性"
echo "  tool_comparison_RiboTISH_specific.tsv     # Ribo-TISH 特异性"
echo "  tool_comparison_Ribotricer_specific.tsv   # Ribotricer 特异性"
echo ""
echo -e "${YELLOW}下一步建议：${NC}"
echo "  1. 查看总结报告: cat tool_comparison_summary_report.txt"
echo "  2. 分析高置信度 ORF: less tool_comparison_high_confidence.tsv"
echo "  3. 可视化结果（见下方 Python 脚本）"
echo ""
