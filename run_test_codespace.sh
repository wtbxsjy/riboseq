#!/bin/bash

# 运行环境: GitHub Codespaces / 低配机器
# 硬件配置: ~2 CPU, ~8GB RAM
# 运行策略:
#   1. 使用 HISAT2 比对器以大幅降低内存消耗 (STAR 需要 >30GB 内存建立索引)
#   2. RPBP 设置 --chains 1 串行运行 MCMC 链，降低峰值内存
#   3. 限制最大内存为 6GB，预留系统内存

echo "Starting test run for Codespace environment (Low Memory)..."

nextflow run . \
    -profile test,singularity \
    --aligner hisat2 \
    --extra_rpbp_args "--chains 1" \
    --max_memory '6.GB' \
    --max_cpus 2 \
    --outdir results_test_codespace \
    -resume

echo "Run completed. Check results in results_test_codespace/"
