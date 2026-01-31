# ============================================================
# Human (人类)
# ============================================================
python3 scripts/prepare_workflow.py \
    -w /path/to/human_workdir \
    -d /path/to/fastq_files \
    --species human \
    --genome GRCh38 \
    --orfquant-container /path/to/orfquant_patched.sif \
    --rpbp-container /path/to/rpbp.sif

# ============================================================
# Mouse (小鼠)
# ============================================================
python3 scripts/prepare_workflow.py \
    -w /path/to/mouse_workdir \
    -d /path/to/fastq_files \
    --species mouse \
    --genome GRCm39 \
    --orfquant-container /path/to/orfquant_patched.sif \
    --rpbp-container /path/to/rpbp.sif

# ============================================================
# Rice (水稻)
# ============================================================
python3 scripts/prepare_workflow.py \
    -w /path/to/rice_workdir \
    -d /path/to/fastq_files \
    --species rice \
    --genome IRGSP-1.0 \
    --orfquant-container /path/to/orfquant_patched.sif \
    --rpbp-container /path/to/rpbp.sif

# ============================================================
# Maize (玉米)
# ============================================================
python3 scripts/prepare_workflow.py \
    -w /path/to/maize_workdir \
    -d /path/to/fastq_files \
    --species maize \
    --genome Zm-B73-REFERENCE-NAM-5.0 \
    --orfquant-container /path/to/orfquant_patched.sif \
    --rpbp-container /path/to/rpbp.sif

# ============================================================
# Wheat (小麦)
# ============================================================
python3 scripts/prepare_workflow.py \
    -w /path/to/wheat_workdir \
    -d /path/to/fastq_files \
    --species wheat \
    --genome IWGSC \
    --orfquant-container /path/to/orfquant_patched.sif \
    --rpbp-container /path/to/rpbp.sif

# ============================================================
# Soybean (大豆)
# ============================================================
python3 scripts/prepare_workflow.py \
    -w /path/to/soybean_workdir \
    -d /path/to/fastq_files \
    --species soybean \
    --genome Glycine_max_v2.1 \
    --orfquant-container /path/to/orfquant_patched.sif \
    --rpbp-container /path/to/rpbp.sif