---
name: riboseq-data-prep
description: >
  为 nf-core/riboseq（本仓库 fork）准备分析项目与输入数据：用 scripts/prepare_workflow.py
  一键创建项目目录/样本表/run 脚本，准备或校验参考基因组与 Ensembl/GENCODE 分类目录
  （5 个标准 symlink），按 PRJ/GEO accession 下载公共数据并交互式标注样本
  （fetch_public_metadata.py + guided tutorial notebook），SRA 转 FASTQ，HQ 样本过滤。
  当用户提到准备样本表(samplesheet)、参考基因组/GTF/FASTA、Ensembl 目录、新建
  riboseq/rice/maize/wheat/soybean/mouse 等项目、PRJEB/PRJNA/GSE/SRR 等公共数据
  下载、SRA 转换、测序数据整理等数据准备任务时使用本 skill——即使没有明确说"数据准备"。
---

# riboseq 数据准备

本仓库的数据准备核心是 `scripts/prepare_workflow.py`（一键准备）+ 若干配套脚本。
所有脚本通过仓库内相对路径引用，不复制进本 skill。

## 决策树

- **全新项目（多样本）** → 用 `prepare_workflow.py` 一键准备（第 1 节）
- **已有项目补组件** → 跳到对应小节：样本表（第 2 节）、参考文件（第 3 节）、Ensembl 目录（第 4 节）、SRA 转换（第 5 节）
- **要快速了解各脚本** → 读 `scripts/README_WORKFLOW_PREP.md` 和 `scripts/PREPARE_WORKFLOW_GUIDE.md`（仓库自带完整文档）

## 1. 一键准备：prepare_workflow.py

**位置**：`scripts/prepare_workflow.py`

**功能**：建 6 目录 → data/ 内建 FASTQ symlink（源不动，避免复制大文件）→ 从 `-r` 复制/链接参考文件（自动解压 .gz）→ 自动准备 Ensembl 目录 → 复制容器 sif → 调 `get_sample_sheet.py` 生成 samplesheet.csv → 生成 `run_pipeline.sh` + `workflow_config.json`。

**标准项目布局**（`run/{project}/`）：

```
run/{project}/
├── data/          # FASTQ symlink → 真实数据目录
├── reference/     # 基因组/GTF/转录本/污染序列 + ensembl/
├── containers/    # *.sif 容器镜像（每项目一份）
├── process/       # Nextflow work 目录 + run_pipeline.sh + 日志
├── result/        # pipeline 输出
└── scripts/       # samplesheet.csv + workflow_config.json + 自定义 config
```

**实际项目完整示例**（rice，源自真实运行）：

```bash
python3 scripts/prepare_workflow.py \
  -w ~/riboseq/run/rice \
  -d ~/riboseq/data/rice \
  -r ~/riboseq/reference/rice \
  --contaminant-dir ~/riboseq/reference/rice \
  --species rice \
  --orfquant-container /path/to/containers/orfquant_patched.sif \
  --unify-orf-container /path/to/containers/unify_orf.sif \
  --gencode-orf-mapper-container /path/to/containers/gencode_orf_mapper.sif \
  --ribowaltz-container /path/to/containers/ribowaltz.sif \
  --price-container /path/to/containers/gedi.sif \
  --skip-collect-qc-stats
```

**最简用法**（不确定时先干跑预览）：

```bash
python3 scripts/prepare_workflow.py -w /path/to/workdir -d /path/to/fastq --dry-run
```

**关键参数速查**（完整列表见 `scripts/PREPARE_WORKFLOW_GUIDE.md`）：

| 参数 | 说明 | 默认 |
|---|---|---|
| `-w/--workdir` | 项目目录（自动创建） | 必需 |
| `-d/--data-dir` | FASTQ 目录（支持 .sra，自动转 FASTQ） | 必需 |
| `-r/--reference-dir` | 参考目录；未指定时自动下载 | None |
| `--contaminant-dir` | 污染索引目录（maize 用 bowtie `contamination_indices`） | None |
| `--species` | human/mouse/rice/maize/wheat/soybean/arabidopsis，驱动 Ensembl 自动准备 | human |
| `--ensembl-release/--ensembl-assembly` | 覆盖默认（plants release 58、ensembl 110） | None |
| `--fasta/--gtf/--transcript-fasta` | 显式指定参考文件（跳过自动匹配） | None |
| `--orf-classify-ensembl-dir` | 已有 Ensembl 目录直接指定；`--skip-orf-classify-ensembl` 关闭自动准备 | None |
| `--strandedness/--sample-type` | 样本表列默认值 | auto/riboseq |
| `--sorf-read-len-min/max` | sORF 读长过滤（用 0/0 关闭长度过滤） | None |
| `--ribotricer-phase-score-cutoff` | 植物低深度建议 0.1 | None |
| `--pathogen-fasta/--pathogen-gtf` | 双基因组（见 CLAUDE.md gotcha 16：建议预先拼接后走 --fasta/--gtf） | None |
| `--*-container` 系列 | 6 个自定义容器路径（orfquant/unify/gencode-mapper/ribowaltz/price/rpbp） | None |
| `--aligner/--max-memory/--max-cpus/--max-time` | 写入 run_pipeline.sh 的资源参数 | star/150.GB/42/48.h |
| `--unify-orf-min-len` | 统一 ORF 最小长度（真实项目用 6 或 24） | 6 |
| `--unify-orf-frame-merge-min-overlap` | frame 合并最小重叠（真实项目用 0.9） | 0.9 |
| `--profile` | nextflow profile | singularity |
| `--sra-threads/--pigz-threads` | SRA 转换线程 | 8/8 |

**两个必知坑**：

1. **生成的 `run_pipeline.sh` 不接受额外参数**——`bash run_pipeline.sh -resume` 的 `-resume` 会被丢弃、从头重跑。需要 resume 时手动复制脚本里的 nextflow 命令加上 `-resume`。
2. **脚本内多行反斜杠续行禁止用 `#` 注释**——注释会吞掉后续所有参数。要禁用某个参数就整行删除。

## 2. 样本表（samplesheet.csv）

**标准 5 列格式**（单端 `fastq_2` 留空）：

```csv
sample,fastq_1,fastq_2,strandedness,type
SRR14340942_1,/home/25119231r/riboseq/run/rice/data/SRR14340942_1.fastq.gz,,auto,riboseq
```

**自动生成**：

```bash
python3 scripts/get_sample_sheet.py -i /path/to/fastq_dir -o samplesheet.csv \
  --strandedness auto --type riboseq
```

**type 列语义**：`riboseq`（完整 RPF 流程）/ `tiseq` / `rnaseq` / `lncrna`。
⚠️ `tiseq` 在 main.nf 中是空壳分支——BAM 被丢弃、无下游消费者。TI-seq 数据要标 `riboseq` 才会走完整流程（GSE157490 的 QTI 样本即如此处理）。

**6 列变体**（含分组信息，用于 --group-map / TE 分析）：`sample,fastq_1,fastq_2,strandedness,type,group`。

**HQ 样本过滤**（低质量样本先剔除再跑 pipeline，否则 RIBOSEQC 会失败）：salmon 唯一比对率 < 20% 的样本应移除。用本 skill 自带的 `scripts/filter_samplesheet_by_rate.py`：

```bash
python3 skills/riboseq-data-prep/scripts/filter_samplesheet_by_rate.py \
  --rates /tmp/alignment_rates.txt \
  --samplesheet samplesheet.csv \
  --output samplesheet_hq.csv \
  --min-rate 20
```

（maize 实际案例：100 → 79 样本；率文件格式为 `rate ... sample` 三列，默认取第 1 列 rate、第 3 列 sample，可用 `--rate-col/--sample-col` 调整。）

⚠️ 样本表变更后 resume 前，必须先清理被移除样本的残留 work dir（详见 riboseq-pipeline-run skill）。

## 3. 参考文件约定

`run/{project}/reference/` 下四个文件（prepare_reference_db_v2.2.py 的输出命名）：

```
{species}.genome.fa               # 基因组
{species}.gtf                     # 注释
{species}.transcripts.fa          # 转录本
{species}_final_contamination.fasta  # 污染物（rRNA/tRNA 等）
```

- **索引不用手工建**：STAR index 由 pipeline 运行时自动构建，`--save_reference true` 会保存到 `result/genome/`（含 index/、.fai、.sizes、filtered.gtf）。
- 污染过滤：`--contaminant_fasta` 传序列文件；或用 bowtie 预建索引目录走 `--contaminant-dir`。
- 参考文件命名不规范时用 prepare_workflow.py 的 `--fasta/--gtf/--transcript-fasta` 显式指定。

## 4. Ensembl/GENCODE 分类目录

GENCODE 分类器（CLASSIFY_ORFS_GENCODE）要求 `--orf_classify_ensembl_dir` 指向一个含 **5 个标准 symlink** 的目录：

```
TRANSCRIPTOME_FASTA       → *.transcriptome.fa / cdna
SORTED_TRANSCRIPTOME_GTF  → 按 chrom + start 排序的 GTF
PROTEOME_FASTA            → *.pep.all.fa
TRANSCRIPT_SUPPORT        → transcript_support_level.txt (TSL/APPRIS)
PSITES_BED                → 每条转录本 1 个 P-site 位置（start_codon 派生）
```

**路线 A（自动，首选）**：prepare_workflow.py 按 `--species` 自动映射 Ensembl 名/division/release/assembly（plants release-58、ensembl 110），下载到 `reference/ensembl/Ens58_{species}/` 并建 symlink。BioMart 的 TSL/APPRIS 下载失败时会自动降级为从 GTF 生成（正常现象）。也可单独用 `scripts/retrieve_ensembl_data.sh`。

**路线 B（本地 GTF 手工构建）**：无 Ensembl 物种或需离线时用 gffread 从本地 GTF+基因组构建 → 完整步骤与坑见 `references/ensembl_dir_manual.md`。

⚠️ 配套容器 `gencode_orf_mapper.sif` 必须自带 **bedtools + BioPython**（unify_orf 容器没有 bedtools，不能混用）。

## 5. 公共数据下载（PRJ ID / GEO accession）

按 PRJ ID（PRJEB/PRJNA/SRP/ERP）、GEO（GSE/GDS）或 run ID（SRR/ERR/DRR）从
公共数据库下载数据并构建样本表。两条路径：

**路径 A：CLI（Claude 直接执行，首选）**——核心引擎是 `scripts/fetch_public_metadata.py`：

```bash
python3 scripts/fetch_public_metadata.py \
  --accession PRJEB26593 --accession GSE149973 \
  --output-prefix run/{project}/scripts/public_metadata \
  --source-strategy ena-first \
  --emit-download-manifest --emit-samplesheet-template
```

产出 `*.metadata_curated.tsv`（28 列，含 inferred_type/suggested_group/
fastq_ftp_1/2/needs_manual_review）→ 人工审阅 type/group 标注 → 按
`fastq_ftp_*` 生成 wget/curl/ascp 下载命令（无 FTP 时 fallback 到
prefetch + fasterq-dump）→ 生成 8 列样本表。

**路径 B：交互式 notebook（用户手动）**——`docs/notebooks/riboseq_guided_tutorial.ipynb`
是同一引擎的 ipywidgets 包装（Step 0-4：accession 输入 → 元数据审阅 →
每 run 的 type 下拉框/group 文本框标注 → 下载脚本 → 样本表 + nextflow 命令组装）。
存在 `inferred_type == "unknown"` 时会阻断样本表生成，保证先审阅再导出。

⚠️ **needs_manual_review=true 的 run 必须人工确认**，unknown 类型清零后才能跑 pipeline。

**本地 .sra 文件转换**：

```bash
bash scripts/sra2fq.sh -t 16 -p 8 -o /output/dir SRR1234567.sra
# -t: fasterq-dump 线程, -p: pigz 压缩线程, 支持多个 .sra 文件
```

下载后校验 `fastq_md5_1/2`。完整细节（元数据列说明、下载三选一命令、
样本表生成规则、notebook 依赖/内核、4 个坑）→ `references/public_data_download.md`。

## 6. 关键坑速查（数据准备范围）

1. **gffread 蛋白 FASTA header 不匹配**（最惨案例）：gffread `-y` 输出 header 是 `transcript_id`（`>Zm..._T001`），GTF 里却是 `protein_id`（`Zm..._P001`）→ protein_seq_map 为空 → 660K ORF 只分类出 2190 个且几乎全是 lncRNA。**必须把蛋白 FASTA header 重命名为 GTF 的 protein_id**。
2. **BED6 未排序** → GENCODE mapper 分类错误，sort 后重跑。
3. **prepare 生成脚本不接受参数 / 行内 `#` 注释吞参数**（见第 1 节）。
4. **reference 目录 symlink 污染**：多物种共目录时残留其他物种的 symlink，逐个 `ls -la` 检查删除。
5. **data/ 自指坏链**：清理时 `find data/ -type l -exec rm {} \;`（symlink 指向自身）。
6. **低质量样本**：salmon 唯一比对率 2.98%/5.80% 级别的样本会拖垮 RIBOSEQC，先 HQ 过滤（>20%）。
7. **workflow_config.json 的 genome 字段可能是默认残留值 "GRCh38"**，不反映真实基因组，别信它。

## 完成数据准备后

下一步交给 **riboseq-pipeline-run** skill（运行/resume/排错）。若用户要跳过运行直接做 ORF 分析，见 **riboseq-orf-analysis** skill 的手动 bypass 路径。
