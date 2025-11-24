#!/bin/bash

# 运行环境: Google Colab Pro / 高配机器
# 硬件配置: ~8 CPU, ~50GB RAM
# 运行策略:
#   1. 使用 STAR 比对器 (Riboseq 黄金标准)，利用大内存优势
#   2. RPBP 设置 --chains 2 (默认) 并行运行，提高稳健性和速度
#   3. 充分利用 8 核 CPU 和 40GB 内存

echo "Starting test run for Colab Pro environment (High Memory)..."

nextflow run . \
    -profile test,singularity \
    --aligner star \
    --extra_rpbp_args "--chains 2" \
    --max_memory '40.GB' \
    --max_cpus 8 \
    --outdir results_test_colab \
    -resume

echo "Run completed. Check results in results_test_colab/"
