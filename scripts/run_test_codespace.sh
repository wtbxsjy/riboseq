#!/bin/bash
#===============================================================================
# GitHub Codespace 测试脚本
#===============================================================================
# 运行环境: GitHub Codespaces
# 硬件配置: 2 CPU, 8GB RAM
# 
# 特点:
#   - 使用 HISAT2 比对器 (内存需求低，~2GB vs STAR 的 30GB+)
#   - RPBP 单链 MCMC 运行，降低峰值内存
#   - 包含 RiboCode ORF 预测
#
# 注意: 测试数据较小，RiboCode 可能因周期性检测失败而报错，这是数据质量问题
#===============================================================================

set -euo pipefail

echo "=============================================="
echo "  GitHub Codespace Test Run"
echo "  Hardware: 2 CPU, 8GB RAM"
echo "  Aligner: HISAT2 (low memory)"
echo "=============================================="

# 清理旧的工作目录（可选，取消注释以启用）
# rm -rf work/ results_codespace/

nextflow run . \
    -profile test_codespace,singularity \
    --outdir results_codespace \
    -resume

echo ""
echo "=============================================="
echo "  Run completed!"
echo "  Results: results_codespace/"
echo "=============================================="
