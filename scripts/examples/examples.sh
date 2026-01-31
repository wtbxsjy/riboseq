#!/bin/bash
# ============================================================
# Ribo-seq Workflow Preparation Examples
# ============================================================
# 支持的数据格式: FASTQ (.fastq.gz/.fq.gz) 或 SRA (.sra)
# SRA 文件会自动转换为 FASTQ.gz
# 不指定 -r 参考目录时会自动下载并解压参考序列
# ============================================================

# ============================================================
# Human (人类) - FASTQ 输入
# ============================================================
python3 scripts/prepare_workflow.py \
    -w /path/to/human_workdir \
    -d /path/to/fastq_files \
    --species human \
    --genome GRCh38 \
    --orfquant-container /path/to/orfquant_patched.sif \
    --rpbp-container /path/to/rpbp.sif \
    --unify-orf-container /path/to/unify_orf.sif

# ============================================================
# Human (人类) - SRA 输入（自动转换）
# ============================================================
python3 scripts/prepare_workflow.py \
    -w /path/to/human_workdir \
    -d /path/to/sra_files \
    --species human \
    --genome GRCh38 \
    --sra-threads 16 \
    --pigz-threads 8 \
    --orfquant-container /path/to/orfquant_patched.sif \
    --rpbp-container /path/to/rpbp.sif \
    --unify-orf-container /path/to/unify_orf.sif

# ============================================================
# Mouse (小鼠)
# ============================================================
python3 scripts/prepare_workflow.py \
    -w /path/to/mouse_workdir \
    -d /path/to/fastq_files \
    --species mouse \
    --genome GRCm39 \
    --orfquant-container /path/to/orfquant_patched.sif \
    --rpbp-container /path/to/rpbp.sif \
    --unify-orf-container /path/to/unify_orf.sif

# ============================================================
# Rice (水稻)
# ============================================================
python3 scripts/prepare_workflow.py \
    -w /path/to/rice_workdir \
    -d /path/to/fastq_files \
    --species rice \
    --genome IRGSP-1.0 \
    --orfquant-container /path/to/orfquant_patched.sif \
    --rpbp-container /path/to/rpbp.sif \
    --unify-orf-container /path/to/unify_orf.sif

# ============================================================
# Maize (玉米)
# ============================================================
python3 scripts/prepare_workflow.py \
    -w /path/to/maize_workdir \
    -d /path/to/fastq_files \
    --species maize \
    --genome Zm-B73-REFERENCE-NAM-5.0 \
    --orfquant-container /path/to/orfquant_patched.sif \
    --rpbp-container /path/to/rpbp.sif \
    --unify-orf-container /path/to/unify_orf.sif

# ============================================================
# Wheat (小麦)
# ============================================================
python3 scripts/prepare_workflow.py \
    -w /path/to/wheat_workdir \
    -d /path/to/fastq_files \
    --species wheat \
    --genome IWGSC \
    --orfquant-container /path/to/orfquant_patched.sif \
    --rpbp-container /path/to/rpbp.sif \
    --unify-orf-container /path/to/unify_orf.sif

# ============================================================
# Soybean (大豆)
# ============================================================
python3 scripts/prepare_workflow.py \
    -w /path/to/soybean_workdir \
    -d /path/to/fastq_files \
    --species soybean \
    --genome Glycine_max_v2.1 \
    --orfquant-container /path/to/orfquant_patched.sif \
    --rpbp-container /path/to/rpbp.sif \
    --unify-orf-container /path/to/unify_orf.sif

# ============================================================
# 使用 SRA 目录作为输入（通用模板）
# ============================================================
# python3 scripts/prepare_workflow.py \
#     -w /path/to/workdir \
#     -d /path/to/sra_directory \
#     --species <species> \
#     --genome <genome_name> \
#     --sra-threads 16 \
#     --pigz-threads 8 \
#     --orfquant-container /path/to/orfquant_patched.sif \
#     --rpbp-container /path/to/rpbp.sif

# ============================================================
# 最简配置（自动下载参考库）
# ============================================================
# python3 scripts/prepare_workflow.py \
#     -w /path/to/workdir \
#     -d /path/to/data \
#     --species human \
#     --genome GRCh38