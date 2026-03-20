# Lab Tutorial / 实验室教学版 Vignette

This tutorial is written for lab members who may have used Jupyter before but are new to bioinformatics pipelines, Nextflow, or Ribo-seq analysis.

这份教学版文档面向组内会一点 Jupyter、但不熟悉生信流程、Nextflow 或 Ribo-seq 分析的同学。目标不是讲清所有参数，而是帮助你快速回答下面几个问题：

1. 这个流程在做什么？
2. 我要准备什么输入？
3. 公共数据库数据怎么找、怎么下？
4. 流程怎么跑、怎么跟踪？
5. 跑完之后先看哪些结果？

## Start Here / 从哪里开始

建议按这个顺序使用：

1. 先打开交互式 notebook：
   [docs/notebooks/riboseq_guided_tutorial.ipynb](notebooks/riboseq_guided_tutorial.ipynb)
2. 再回到这份文档看背景说明和结果解释。
3. 真正运行前，再查正式参考文档：
   [usage.md](usage.md)
   [output.md](output.md)

## What The Pipeline Does / 这个流程做什么

`nf-core/riboseq` 是一个 Nextflow DSL2 流程，用来分析核糖体测序数据（Ribo-seq），也兼容部分 RNA-seq 配套输入。它把“原始 reads 或 BAM”一路处理成：

- 预处理后的 clean reads / filtered BAM
- 比对结果
- Ribo-seq 质量控制结果
- 多个 ORF calling 工具的预测结果
- 统一 ORF 集合（unified ORFs）
- ORF 分类结果
- MultiQC 总览报告

对初学者来说，可以把它理解成 5 个层次：

1. `Data preparation`
   - 识别样本、准备 samplesheet、确认 `type` 和 `group`
2. `Preprocessing and alignment`
   - 做 read QC、去接头、去污染、比对、排序
3. `Ribo-seq specific QC`
   - 看长度分布、P-site、周期性、frame bias
4. `ORF prediction and post-processing`
   - 不同工具预测 ORF，然后做统一和分类
5. `Reporting`
   - MultiQC、结果目录、统计摘要

## Input Modes / 输入模式

这个流程常见有两种入口：

### 1. FASTQ mode

适合你只有原始测序数据时使用。流程会从 raw FASTQ 开始做预处理、比对和后续分析。

你需要准备：

- `samplesheet.csv`
- 参考基因组和注释文件
- 一个可用运行环境（Docker / Singularity / Conda）

### 2. BAM mode

适合你已经有对参考基因组比对好的 BAM 文件时使用。

你需要准备：

- BAM 和 BAI
- `samplesheet_bam.csv`
- 明确指定 strandedness

## Public Data Retrieval / 从 ENA 和 NCBI 获取公开数据

这次教学材料新增了一层“公开数据库数据发现与 metadata 整理”。核心思路是：

- 优先查 ENA，因为它常常直接给出 `fastq_ftp`
- ENA 不完整时，用 NCBI SRA 补字段或作为下载兜底
- `type` 和 `group` 只做保守推断，必须人工确认

### 常见 accession 类型

- Run: `SRR`, `ERR`, `DRR`
- Experiment: `SRX`, `ERX`, `DRX`
- Study / Project: `SRP`, `ERP`, `DRP`, `PRJNA`, `PRJEB`
- Sample / BioSample: `SRS`, `ERS`, `DRS`, `SAMN`, `SAMEA`

### 新增脚本

仓库新增了：

- [scripts/fetch_public_metadata.py](../scripts/fetch_public_metadata.py)

它负责：

- accession 展开到 run-level metadata
- 汇总 ENA / NCBI 字段
- 推断 `inferred_type`
- 给出 `suggested_group`
- 输出下载清单和候选 samplesheet

### 推荐命令

```bash
python scripts/fetch_public_metadata.py \
  --accession SRR15480782 \
  --accession SRR15480788 \
  --accession SRR11192680 \
  --source-strategy ena-first \
  --output-prefix tutorial_outputs/public_metadata \
  --emit-download-manifest \
  --emit-samplesheet-template
```

### 关键输出文件

- `*.metadata_raw.tsv`
  - 原始来源字段，适合排错和追溯
- `*.metadata_raw.json`
  - 原始结构化记录
- `*.metadata_curated.tsv`
  - notebook 和后续人工整理的主表
- `*.downloads.tsv`
  - 下载清单
- `*.samplesheet.csv`
  - 候选 samplesheet
- `*.warnings.tsv`
  - 哪些 accession 缺字段、失败或需要注意

### 怎样理解 metadata 关键列

- `library_strategy`
  - 最重要的 assay 提示之一。`RNA-Seq` 通常直接支持判为 `rnaseq`
- `library_layout`
  - `SINGLE` 或 `PAIRED`
- `sample_title`, `experiment_title`, `study_title`
  - 用来辅助判断 Ribo-seq / RNA-seq，以及识别 control、treated、replicate 等线索
- `inferred_type`
  - 脚本的保守推断结果：`riboseq`, `rnaseq`, `tiseq`, `unknown`
- `suggested_group`
  - 脚本根据标题和样本名给出的候选分组
- `needs_manual_review`
  - 为 `true` 时，说明你不应该直接把它静默带入正式分析

### 为什么不能完全自动化

公共数据库的命名很不统一，所以 `type` 和 `group` 不能黑盒自动化。常见问题包括：

- Ribo-seq 的 `library_strategy` 经常写成 `OTHER`
- study title 里可能既有 RNA-seq 也有 Ribo-seq 信息
- 分组可能藏在自由文本里
- 数字后缀可能代表 replicate，也可能代表剂量或时间点

### 人工确认 checklist

正式跑流程前，请确认：

1. 这个样本真的是 `riboseq` 还是 `rnaseq`？
2. 是否混入了别的 assay（例如 amplicon、ChIP、其他转录组数据）？
3. `suggested_group` 是否真的对应你的 biological condition？
4. 数字后缀到底是 replicate、dose、timepoint，还是纯编号？

## Build A Samplesheet / 生成 samplesheet

你可以用两种方式得到 samplesheet：

### 1. 从公开 metadata 生成

最适合刚从 ENA / NCBI 挑完数据的时候。先用 notebook 或 `fetch_public_metadata.py` 生成候选表，再手工确认。

### 2. 从本地 FASTQ 目录扫描

仓库已有：

- [scripts/get_sample_sheet.py](../scripts/get_sample_sheet.py)

适合你已经下载好 FASTQ，想从目录自动拼一个样本表：

```bash
python scripts/get_sample_sheet.py \
  -i /path/to/fastq_dir \
  -o samplesheet.csv \
  --strandedness auto \
  --type riboseq
```

## Recommended Onboarding Workflow / 推荐的新手上手顺序

### Step 1. 用 notebook 查公开数据

打开：

- [docs/notebooks/riboseq_guided_tutorial.ipynb](notebooks/riboseq_guided_tutorial.ipynb)

在 notebook 里：

1. 输入 accession
2. 选择 `ena-first` 或别的源策略
3. 查看 metadata 表
4. 手工确认 `type` 和 `group`

### Step 2. 生成下载命令

notebook 会给出几种常见下载方式：

- `wget`
- `curl`
- `ascp`
- `prefetch + fasterq-dump`

默认建议先生成命令，检查无误后再执行。

### Step 3. 生成 samplesheet

当 `unknown` 行都修正后，再生成 `samplesheet.csv`。

### Step 4. 生成并运行 Nextflow 命令

最小示例：

```bash
nextflow run . \
  -profile docker \
  --input samplesheet.csv \
  --outdir results
```

常见变体：

```bash
nextflow run . -profile singularity --input samplesheet.csv --outdir results
nextflow run . -profile docker --input samplesheet.csv --aligner hisat2 --outdir results
nextflow run . -profile docker --input samplesheet.csv --outdir results -resume
```

### Step 5. 跟踪运行状态

新手最值得看的地方：

- 终端里的 Nextflow 进度
- `work/`
- `results/`
- `pipeline_info/`
- MultiQC

如果失败：

1. 看报错的 process 名称
2. 修参数或环境问题
3. 用 `-resume` 重跑

## Key Pipeline Concepts / 流程里最容易混淆的概念

### 1. `type`

`type` 是 samplesheet 里最重要的字段之一。

- `riboseq`
  - 会走 Ribo-seq 专属 QC 和 ORF 预测路径
- `rnaseq`
  - 主要作为配套或未来 TE 分析输入
- `tiseq`
  - Translation initiation related 数据

### 2. `group`

`group` 用来表示 biological/technical grouping。它尤其影响：

- replicate merge
- 结果解释
- 后续条件比较

### 3. `merge_replicates`

如果多个样本共享同一个 `group`，可通过：

```bash
--merge_replicates
```

在 ORF calling 前合并过滤后的 BAM，提高信号强度。

## Result Overview / 结果怎么看

对初学者来说，最有价值的结果通常分 4 类。

### 1. MultiQC

这是全流程总览入口，适合先快速判断：

- reads 质量是否正常
- 比对率是否离谱
- 是否有明显污染或文库问题

### 2. Ribo-seq QC

重点关注：

- read length distribution
- P-site offset
- periodicity
- frame bias

这些结果能帮助你判断样本是否像“真正的 Ribo-seq”。

### 3. ORF prediction outputs

不同工具会给出不同视角的 ORF 候选。最重要的是理解：

- 哪些 ORF 被多个工具共同支持
- 哪些 ORF 只有单一工具支持

### 4. Unified ORF results

这是教学材料里最推荐初学者阅读的一层，因为它已经把多个工具结果做了统一整理。

本仓库 demo 数据位于：

- [test_data/tutorial_demo_public_data/unified_orfs_demo.tsv](../test_data/tutorial_demo_public_data/unified_orfs_demo.tsv)
- [test_data/tutorial_demo_public_data/unified_orfs.stats.txt](../test_data/tutorial_demo_public_data/unified_orfs.stats.txt)

`unified_orfs_demo.tsv` 里最适合先看的列：

- `orf_id`
- `tools`
- `samples`
- `length_aa`
- `start_codon`
- `unique_psites`
- `pN`
- `overlapping_genes`

你可以把它们粗略理解为：

- `tools`
  - 哪些工具支持这个 ORF
- `samples`
  - 哪些样本检测到它
- `unique_psites`
  - 支持信号强不强
- `pN`
  - 一个统计支持指标

`unified_orfs.stats.txt` 则更像摘要：

- 原始 ORF 数量有多少
- 每个工具各自贡献多少
- 合并后还剩多少

## Demo Assets / 这次新增的教学示例

新增的轻量教学数据位于：

- [test_data/tutorial_demo_public_data/accessions.txt](../test_data/tutorial_demo_public_data/accessions.txt)
- [test_data/tutorial_demo_public_data/metadata_raw.tsv](../test_data/tutorial_demo_public_data/metadata_raw.tsv)
- [test_data/tutorial_demo_public_data/metadata_raw.json](../test_data/tutorial_demo_public_data/metadata_raw.json)
- [test_data/tutorial_demo_public_data/metadata_curated.tsv](../test_data/tutorial_demo_public_data/metadata_curated.tsv)
- [test_data/tutorial_demo_public_data/downloads.tsv](../test_data/tutorial_demo_public_data/downloads.tsv)
- [test_data/tutorial_demo_public_data/samplesheet_from_metadata.csv](../test_data/tutorial_demo_public_data/samplesheet_from_metadata.csv)
- [test_data/tutorial_demo_public_data/unified_orfs_demo.tsv](../test_data/tutorial_demo_public_data/unified_orfs_demo.tsv)

它们覆盖了：

- `rnaseq`
- `riboseq`
- 一个会落到 `unknown` 的样本
- 带有 group 提示的样本名

## Quick Checklist / 快速清单

第一次上手时，请按这个清单走：

1. 打开 notebook，输入 accession
2. 检查 `inferred_type`
3. 检查 `suggested_group`
4. 生成下载命令
5. 生成 samplesheet
6. 生成 `nextflow run` 命令
7. 运行流程
8. 先看 MultiQC，再看 unified ORF demo 对应的真实结果

## Related References / 相关正式文档

- [usage.md](usage.md)
- [output.md](output.md)
- [README.md](../README.md)
