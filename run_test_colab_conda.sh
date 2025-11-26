#!/bin/bash
#===============================================================================
# Google Colab Pro 测试脚本 (Conda 版本)
#===============================================================================
# 运行环境: Google Colab Pro (不支持 Docker/Singularity)
# 硬件配置: 8 CPU, 50GB RAM
# 依赖管理: Conda/Mamba
#
# 特点:
#   - 使用 Conda 替代容器化方案 (Colab 不支持 Docker/Singularity)
#   - 使用 STAR 比对器 (Ribo-seq 黄金标准)
#   - RPBP 双链 MCMC 并行运行
#   - 完整 ORF 预测流程 (Ribotish + Ribotricer + RPBP + RiboCode)
#   - RiboseQC 质控分析
#
# 使用方法 (在 Colab notebook 中):
#   !git clone https://github.com/your-repo/riboseq.git
#   %cd riboseq
#   !bash run_test_colab_conda.sh
#
# 注意:
#   - 首次运行需要较长时间下载和创建 Conda 环境
#   - 测试数据较小，RiboCode 可能因周期性检测失败而报错
#===============================================================================

set -euo pipefail

echo "=============================================="
echo "  Google Colab Pro Test Run (Conda)"
echo "  Hardware: 8 CPU, 50GB RAM"
echo "  Aligner: STAR (gold standard)"
echo "  Package Manager: Conda/Mamba"
echo "=============================================="

# 检查 Nextflow 是否已安装
if ! command -v nextflow &> /dev/null; then
    echo "Installing Nextflow..."
    curl -fsSL get.nextflow.io | bash
    chmod +x nextflow
    export PATH=$PATH:$(pwd)
fi

# 检查 Mamba/Conda 是否已安装
if ! command -v mamba &> /dev/null && ! command -v conda &> /dev/null; then
    echo "Installing Miniforge (Mamba)..."
    wget -q "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh" -O miniforge.sh
    bash miniforge.sh -b -p $HOME/miniforge3
    rm miniforge.sh
    export PATH=$HOME/miniforge3/bin:$PATH
    mamba init bash
    source ~/.bashrc
fi

# 显示环境信息
echo ""
echo "Environment Info:"
echo "  Nextflow: $(nextflow -version 2>&1 | head -1)"
if command -v mamba &> /dev/null; then
    echo "  Mamba: $(mamba --version | head -1)"
elif command -v conda &> /dev/null; then
    echo "  Conda: $(conda --version)"
fi
echo ""

# 清理旧的工作目录（可选，取消注释以启用）
# rm -rf work/ results_colab_conda/

# 运行 pipeline
nextflow run . \
    -profile test_colab_conda \
    --outdir results_colab_conda \
    -resume

echo ""
echo "=============================================="
echo "  Run completed!"
echo "  Results: results_colab_conda/"
echo "=============================================="
