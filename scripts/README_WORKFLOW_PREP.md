# Workflow Preparation Scripts

## 概览

这个目录包含用于准备 Ribo-seq 分析工作流的脚本集合。

## 脚本清单

### 1. prepare_workflow.py ⭐ (主要脚本)

**功能**: 自动化工作流准备，一键建立完整分析环境

**核心特性**:
- ✅ 创建标准化目录结构 (data, reference, containers, process, result, scripts)
- ✅ 智能符号链接管理（避免数据复制）
- ✅ 自动生成样本表 (samplesheet.csv)
- ✅ 生成可执行的 Nextflow 脚本
- ✅ 保存配置摘要 (JSON)
- ✅ 支持干运行模式预览

**快速开始**:
```bash
# 基础用法
python3 scripts/prepare_workflow.py \
    -w /path/to/workdir \
    -d /path/to/fastq_files

# 完整配置
python3 scripts/prepare_workflow.py \
    -w ~/riboseq_project \
    -d /data/fastq \
    -r /data/reference \
    -c /data/containers \
    --genome GRCh38 \
    --species human \
    --orfquant-container /path/to/orfquant_patched.sif \
    --run-prefilter-qc
```

**详细文档**: 查看 [PREPARE_WORKFLOW_GUIDE.md](PREPARE_WORKFLOW_GUIDE.md)

---

### 2. get_sample_sheet.py

**功能**: 从 FASTQ 文件目录生成样本表

**用法**:
```bash
python3 scripts/get_sample_sheet.py \
    -i /path/to/fastq_dir \
    -o samplesheet.csv \
    --strandedness auto \
    --type riboseq
```

**特性**:
- 自动识别单端/双端数据
- 智能提取样本名称
- 支持多种文件命名格式 (_R1/_R2, _1/_2)

**输出示例**:
```csv
sample,fastq_1,fastq_2,strandedness,type
sample1,/path/sample1_R1.fastq.gz,/path/sample1_R2.fastq.gz,auto,riboseq
sample2,/path/sample2_R1.fastq.gz,/path/sample2_R2.fastq.gz,auto,riboseq
```

---

### 3. prepare_reference_db_v2.2.py

**功能**: 准备参考数据库（rRNA, tRNA, 基因组等）

**用法**:
```bash
python3 scripts/prepare_reference_db_v2.2.py \
    --species human \
    --output-dir /path/to/reference \
    --download-genome \
    --download-rrna
```

**支持物种**:
- 动物: human, mouse
- 植物: rice, maize, wheat

**输出**:
- 基因组 FASTA
- GTF 注释文件
- 转录本序列
- rRNA/tRNA 序列
- 污染物序列

---

### 4. sra2fq.sh

**功能**: 批量转换 SRA 文件为 FASTQ.gz

**用法**:
```bash
# 单个文件
bash scripts/sra2fq.sh -t 16 -p 8 -o /output/dir SRR1234567.sra

# 多个文件
bash scripts/sra2fq.sh -t 16 -p 8 -o /output/dir *.sra
```

**参数**:
- `-t`: fasterq-dump 线程数
- `-p`: pigz 压缩线程数
- `-o`: 输出目录

---

### 5. quick_start_example.sh

**功能**: 快速开始示例脚本

**用法**:
1. 编辑脚本设置数据目录:
   ```bash
   nano scripts/quick_start_example.sh
   # 修改 DATA_DIR="/path/to/your/fastq_files"
   ```

2. 运行:
   ```bash
   bash scripts/quick_start_example.sh
   ```

---

## 完整工作流示例

### 场景 1: 从原始数据开始

```bash
# Step 1: 准备参考数据库
python3 scripts/prepare_reference_db_v2.2.py \
    --species human \
    --output-dir ~/references/human_gencode_v49

# Step 2: （可选）从 SRA 下载数据
bash scripts/sra2fq.sh -o ~/data/fastq SRR*.sra

# Step 3: 准备工作流
python3 scripts/prepare_workflow.py \
    -w ~/riboseq_analysis \
    -d ~/data/fastq \
    -r ~/references/human_gencode_v49 \
    --genome GRCh38 \
    --species human

# Step 4: 运行分析
cd ~/riboseq_analysis
bash run_pipeline.sh
```

### 场景 2: 使用已有数据和参考

```bash
# 直接准备工作流
python3 scripts/prepare_workflow.py \
    -w /work/my_project \
    -d /data/existing_fastq \
    -r /data/existing_reference \
    -c /data/containers \
    --genome GRCh38 \
    --orfquant-container /data/containers/orfquant_patched.sif \
    --run-prefilter-qc \
    --unify-orf-min-len 24

# 运行
cd /work/my_project
bash run_pipeline.sh
```

### 场景 3: 测试运行（干运行）

```bash
# 先预览
python3 scripts/prepare_workflow.py \
    -w /work/test \
    -d /data/test_fastq \
    --dry-run

# 确认无误后执行
python3 scripts/prepare_workflow.py \
    -w /work/test \
    -d /data/test_fastq
```

---

## 目录结构说明

使用 `prepare_workflow.py` 后会创建以下结构：

```
workdir/
├── data/                    # FASTQ 数据（符号链接）
│   ├── sample1_R1.fastq.gz -> /original/path/sample1_R1.fastq.gz
│   └── sample1_R2.fastq.gz -> /original/path/sample1_R2.fastq.gz
│
├── reference/               # 参考文件（符号链接或实际文件）
│   ├── genome.fa.gz
│   ├── annotation.gtf.gz
│   └── contaminant.fa
│
├── containers/              # 容器镜像（符号链接）
│   ├── orfquant_patched.sif -> /containers/orfquant_patched.sif
│   └── rpbp.sif -> /containers/rpbp.sif
│
├── process/                 # Nextflow 工作目录
│   └── work/               # (运行时创建)
│
├── result/                  # 流程输出
│   ├── orf_predictions/    # (运行后生成)
│   ├── orf_unification/
│   └── pipeline_report.html
│
└── scripts/                 # 脚本和配置
    ├── samplesheet.csv          # 样本表
    ├── workflow_config.json     # 配置摘要
    └── run_pipeline.sh          # 执行脚本 ⭐
```

---

## 快速参考

### 常用命令组合

```bash
# 1. 准备工作流（最小配置）
python3 scripts/prepare_workflow.py -w workdir -d fastq_dir

# 2. 准备工作流（完整配置）
python3 scripts/prepare_workflow.py \
    -w workdir -d fastq_dir -r ref_dir -c container_dir \
    --genome GRCh38 --run-prefilter-qc

# 3. 生成样本表
python3 scripts/get_sample_sheet.py -i fastq_dir -o samplesheet.csv

# 4. 准备参考数据库
python3 scripts/prepare_reference_db_v2.2.py --species human --output-dir ref_dir

# 5. SRA 转 FASTQ
bash scripts/sra2fq.sh -t 16 -o output_dir *.sra
```

### 参数速查

| 脚本 | 关键参数 | 说明 |
|------|---------|------|
| prepare_workflow.py | `-w, -d` | 工作目录, 数据目录 (必需) |
| | `-r, -c` | 参考目录, 容器目录 (可选) |
| | `--genome` | 基因组名称 |
| | `--run-prefilter-qc` | 启用 prefilter QC |
| get_sample_sheet.py | `-i, -o` | 输入目录, 输出文件 |
| prepare_reference_db_v2.2.py | `--species` | 物种名称 |
| sra2fq.sh | `-t, -p, -o` | 线程, 压缩线程, 输出 |

---

## 故障排除

### 问题 1: 找不到 FASTQ 文件
```bash
# 检查数据目录
ls -lh /path/to/fastq_dir/*.fastq.gz

# 确保文件扩展名正确 (.fastq.gz 或 .fq.gz)
```

### 问题 2: 符号链接失败
```bash
# 检查源文件权限
ls -l /source/file

# 使用绝对路径
python3 scripts/prepare_workflow.py -w $(pwd)/workdir -d /full/path/to/data
```

### 问题 3: 样本表为空
```bash
# 手动测试 get_sample_sheet.py
python3 scripts/get_sample_sheet.py -i fastq_dir -o test.csv

# 检查输出
cat test.csv
```

### 问题 4: 容器路径错误
```bash
# 验证容器文件存在
ls -lh /path/to/orfquant_patched.sif

# 使用绝对路径
--orfquant-container $(realpath /path/to/orfquant_patched.sif)
```

---

## 更多信息

- **详细指南**: [PREPARE_WORKFLOW_GUIDE.md](PREPARE_WORKFLOW_GUIDE.md)
- **Pipeline 文档**: [../docs/usage.md](../docs/usage.md)
- **输出说明**: [../docs/output.md](../docs/output.md)

---

## 更新日志

### 2026-01-31
- ✨ 新增 `prepare_workflow.py` - 主要工作流准备脚本
- 📝 新增详细使用指南
- 🚀 新增快速开始示例
- 🔧 整合现有脚本功能

---

## 贡献

如有问题或建议，请提交 Issue 或 Pull Request。

## 许可证

与主项目相同。
