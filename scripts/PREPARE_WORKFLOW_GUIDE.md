# Prepare Workflow Script - 使用指南

## 概述

`prepare_workflow.py` 是一个自动化工作流准备脚本，用于快速建立 Ribo-seq 分析的标准化工作环境。

## 功能特性

### 1. 目录结构创建
自动创建标准化的项目目录结构：
```
workdir/
├── data/          # 原始测序数据 (FASTQ 文件)
├── reference/     # 参考基因组、GTF、转录本文件
├── containers/    # Singularity 容器镜像
├── process/       # Nextflow 工作目录和日志
├── result/        # 流程输出结果
└── scripts/       # 分析脚本和配置文件
    ├── samplesheet.csv          # 样本表
    ├── workflow_config.json     # 配置摘要
    └── run_pipeline.sh          # 执行脚本
```

### 2. 符号链接管理
- 自动为数据文件创建软链接（避免复制大文件）
- 支持参考文件和容器镜像的符号链接
- 智能文件类型识别

### 3. 样本表生成
- 调用 `get_sample_sheet.py` 自动生成样本表
- 自动识别单端/双端数据
- 支持自定义 strandedness 和样本类型

### 4. Nextflow 脚本生成
- 生成可执行的 Nextflow 运行脚本
- 包含推荐的参数配置
- 自动配置容器路径和参考文件路径

## 使用方法

### 基础用法

```bash
# 最简单的用法（仅指定数据目录）
python3 scripts/prepare_workflow.py \
    -w /path/to/workdir \
    -d /path/to/fastq_files
```

### 完整配置

```bash
python3 scripts/prepare_workflow.py \
    -w /home/user/riboseq_project \
    -d /data/fastq \
    -r /data/reference \
    -c /data/containers \
    --genome GRCh38 \
    --species human \
    --orfquant-container /path/to/orfquant_patched.sif \
    --rpbp-container /path/to/rpbp.sif \
    --run-prefilter-qc \
    --unify-orf-min-len 6
```

### 干运行模式（查看将要执行的操作）

```bash
python3 scripts/prepare_workflow.py \
    -w /path/to/workdir \
    -d /path/to/fastq_files \
    --dry-run
```

## 参数说明

### 必需参数

| 参数 | 说明 |
|------|------|
| `-w, --workdir` | 工作目录（将自动创建） |
| `-d, --data-dir` | 包含 FASTQ 文件的目录 |

### 可选参数 - 数据源

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `-r, --reference-dir` | 参考文件目录 | None |
| `-c, --container-dir` | 容器镜像目录 | None |

### 可选参数 - 基因组设置

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `--genome` | 基因组名称 | GRCh38 |
| `--species` | 物种 (human/mouse/rice/maize/wheat) | human |

### 可选参数 - 样本表

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `--strandedness` | 链特异性 | auto |
| `--sample-type` | 样本类型 (riboseq/tiseq) | riboseq |

### 可选参数 - 容器

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `--orfquant-container` | ORFquant 容器路径 | None |
| `--rpbp-container` | RPBP 容器路径 | None |

### 可选参数 - 流程配置

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `--run-prefilter-qc` | 启用 prefilter QC 分析 | False |
| `--unify-orf-min-len` | 统一 ORF 预测的最小长度 | 6 |
| `--profile` | Nextflow profile | singularity |

### 可选参数 - 其他

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `--script-name` | 生成的执行脚本名称 | run_pipeline.sh |
| `--dry-run` | 干运行模式 | False |

## 使用示例

### 示例 1: 本地数据快速准备

```bash
# 场景：本地有 FASTQ 文件，需要快速建立分析环境
python3 scripts/prepare_workflow.py \
    -w ~/riboseq_analysis \
    -d /data/raw_fastq \
    --genome GRCh38 \
    --species human
```

执行后：
```
~/riboseq_analysis/
├── data/ -> /data/raw_fastq/*.fastq.gz (软链接)
├── scripts/
│   ├── samplesheet.csv (自动生成)
│   └── run_pipeline.sh (可执行)
└── [其他标准目录]
```

运行流程：
```bash
cd ~/riboseq_analysis
bash run_pipeline.sh
```

### 示例 2: 使用自定义容器和参考基因组

```bash
python3 scripts/prepare_workflow.py \
    -w /work/riboseq_project \
    -d /data/fastq \
    -r /data/gencode_v49 \
    --orfquant-container /containers/orfquant_patched.sif \
    --rpbp-container /containers/rpbp.sif \
    --genome GRCh38 \
    --run-prefilter-qc \
    --unify-orf-min-len 6
```

### 示例 3: 小鼠数据分析

```bash
python3 scripts/prepare_workflow.py \
    -w /work/mouse_riboseq \
    -d /data/mouse_fastq \
    -r /data/gencode_mouse_M38 \
    --genome GRCm39 \
    --species mouse \
    --orfquant-container /containers/orfquant_patched.sif
```

### 示例 4: 植物（水稻）数据分析

```bash
python3 scripts/prepare_workflow.py \
    -w /work/rice_riboseq \
    -d /data/rice_fastq \
    -r /data/ensembl_rice \
    --genome IRGSP-1.0 \
    --species rice \
    --strandedness reverse
```

### 示例 5: 干运行模式（预览操作）

```bash
# 先预览将要执行的操作
python3 scripts/prepare_workflow.py \
    -w /work/test_project \
    -d /data/test_fastq \
    --dry-run

# 确认无误后，去掉 --dry-run 正式执行
python3 scripts/prepare_workflow.py \
    -w /work/test_project \
    -d /data/test_fastq
```

## 输出文件说明

### 1. 目录结构

所有标准目录会自动创建，带有描述性日志说明每个目录的用途。

### 2. samplesheet.csv

自动生成的样本表，格式如下：
```csv
sample,fastq_1,fastq_2,strandedness,type
sample1,/path/to/sample1_R1.fastq.gz,/path/to/sample1_R2.fastq.gz,auto,riboseq
sample2,/path/to/sample2_R1.fastq.gz,/path/to/sample2_R2.fastq.gz,auto,riboseq
```

### 3. run_pipeline.sh

可执行的 Nextflow 脚本，包含：
- 环境设置
- 完整的 Nextflow 命令（带所有参数）
- 报告生成配置
- 结果路径说明

示例内容：
```bash
#!/bin/bash
nextflow run main.nf \
    -profile singularity \
    -w /work/riboseq/process/work \
    --input /work/riboseq/scripts/samplesheet.csv \
    --outdir /work/riboseq/result \
    --genome GRCh38 \
    --orfquant_container /work/riboseq/containers/orfquant_patched.sif \
    --skip_unify_orf_predictions false \
    --skip_orf_classification false \
    --unify_orf_min_len 6 \
    -resume \
    -with-report result/pipeline_report.html
```

### 4. workflow_config.json

JSON 格式的配置摘要，记录所有设置：
```json
{
  "workflow_setup": {
    "created_at": "2026-01-31T10:30:00",
    "working_directory": "/work/riboseq"
  },
  "data": {
    "source_directory": "/data/fastq",
    "sample_sheet": "/work/riboseq/scripts/samplesheet.csv"
  },
  "containers": {
    "orfquant": "/work/riboseq/containers/orfquant_patched.sif"
  },
  "pipeline_options": {
    "run_prefilter_qc": false,
    "unify_orf_min_len": 6
  }
}
```

## 工作流程

1. **准备阶段**（使用本脚本）
   ```bash
   python3 scripts/prepare_workflow.py -w workdir -d data_dir
   ```

2. **执行阶段**
   ```bash
   cd workdir
   bash run_pipeline.sh
   ```

3. **结果查看**
   ```
   workdir/result/
   ├── orf_predictions/
   │   ├── ribotish/postfilter/
   │   ├── ribotricer/postfilter/
   │   └── orfquant/postfilter/
   ├── orf_unification/
   └── orf_classification/
   ```

## 常见问题

### Q1: 数据文件很大，会被复制吗？
**A**: 不会。脚本使用符号链接（symlink），不会复制数据文件，仅创建指向原文件的链接。

### Q2: 如何添加更多样本？
**A**: 
1. 将新的 FASTQ 文件放入原数据目录
2. 重新运行 `prepare_workflow.py`（会跳过已存在的链接）
3. 或手动编辑 `samplesheet.csv`

### Q3: 可以在不同服务器上使用吗？
**A**: 可以，但需要确保：
- 符号链接指向的源文件可访问
- 容器镜像路径正确
- Nextflow 和依赖工具已安装

### Q4: 如何更改流程参数？
**A**: 两种方式：
1. 重新运行 `prepare_workflow.py` 并指定新参数
2. 直接编辑生成的 `run_pipeline.sh` 脚本

### Q5: 干运行模式有什么用？
**A**: 
- 预览将要创建的目录和链接
- 检查参数是否正确
- 验证源文件是否存在
- 查看将要生成的脚本内容

## 高级用法

### 自定义执行脚本

如果需要更多控制，可以：

1. 使用 `--dry-run` 查看生成的脚本内容
2. 正式运行生成脚本
3. 手动编辑 `run_pipeline.sh` 添加自定义参数

### 批量处理多个项目

创建批处理脚本：
```bash
#!/bin/bash
# batch_prepare.sh

projects=(
    "project1:/data/project1_fastq"
    "project2:/data/project2_fastq"
    "project3:/data/project3_fastq"
)

for project in "${projects[@]}"; do
    name="${project%%:*}"
    data="${project##*:}"
    
    python3 scripts/prepare_workflow.py \
        -w "/work/${name}" \
        -d "${data}" \
        --genome GRCh38 \
        --species human
done
```

### 与版本控制集成

```bash
# 在工作目录初始化 git（可选）
cd workdir
git init
git add scripts/samplesheet.csv scripts/workflow_config.json
git commit -m "Initial workflow setup"
```

## 依赖要求

- Python 3.6+
- `scripts/get_sample_sheet.py` (自动调用)
- Nextflow (运行生成的脚本时需要)
- Singularity/Docker (运行流程时需要)

## 许可证

与主流程相同的许可证。

## 更新日志

### v1.0 (2026-01-31)
- 初始版本
- 支持标准目录结构创建
- 自动样本表生成
- Nextflow 脚本生成
- 配置摘要输出
