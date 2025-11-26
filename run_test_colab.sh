#!/bin/bash
#===============================================================================
# Google Colab Pro 测试脚本
#===============================================================================
# 运行环境: Google Colab Pro
# 硬件配置: 8 CPU, 50GB RAM
# 
# 特点:
#   - 使用 STAR 比对器 (Ribo-seq 黄金标准)
#   - RPBP 双链 MCMC 并行运行，更快更稳健
#   - 完整 ORF 预测流程 (Ribotish + RPBP + RiboCode)
#
# 注意: 测试数据较小，RiboCode 可能因周期性检测失败而报错，这是数据质量问题
#===============================================================================

set -euo pipefail

echo "=============================================="
echo "  Google Colab Pro Test Run"
echo "  Hardware: 8 CPU, 50GB RAM"
echo "  Aligner: STAR (gold standard)"
echo "=============================================="

# 清理旧的工作目录（可选，取消注释以启用）
# rm -rf work/ results_colab/

nextflow run . \
    -profile test_colab,singularity \
    --outdir results_colab \
    -resume

echo ""
echo "=============================================="
echo "  Run completed!"
echo "  Results: results_colab/"
echo "=============================================="
