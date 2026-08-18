# 公共数据下载（PRJ/GEO accession）与交互式样本标注

## 组件关系

```
docs/notebooks/riboseq_guided_tutorial.ipynb   ← 交互式 UI（ipywidgets，Step 0-4）
        └─ riboseq_guided_tutorial_helpers.py ← widget 逻辑
              └─ scripts/fetch_public_metadata.py ← 核心引擎（CLI，可脱离 notebook 使用）
```

- **核心引擎**：`scripts/fetch_public_metadata.py` 查询 ENA/NCBI 元数据并生成所有辅助表。
- **notebook** 是给用户在 Jupyter/VS Code 里交互式审阅标注用的包装（每 run 一个 type 下拉框 + group 文本框）。
- Claude Code 场景下**直接走 CLI 路径**即可，notebook 留给用户手动交互时使用。

## 支持的 accession

- **Project ID**：PRJEB/PRJNA/SRP/ERP（按 project 枚举全部 run）
- **GEO**：GSE/GDS（`resolve_geo_to_sra` 先解析到 SRA study，再经 ENA 枚举 run）
- **Run**：SRR/ERR/DRR 单个或多个
- 批量：`--accession-file` 每行一个（空行和 # 注释忽略）

## 元数据查询（CLI）

```bash
python3 scripts/fetch_public_metadata.py \
  --accession PRJEB26593 --accession GSE149973 \
  --output-prefix run/human_PRJEB26593/scripts/public_metadata \
  --source-strategy ena-first \
  --emit-download-manifest \
  --emit-samplesheet-template
```

参数：

| 参数 | 说明 | 默认 |
|---|---|---|
| `--accession` | 可多次传入，或 `--accession-file` 批量 | — |
| `--source-strategy` | ena-first / ncbi-first / ena-only / ncbi-only | ena-first |
| `--output-prefix` | 输出文件前缀 | public_metadata |
| `--emit-download-manifest` | 生成下载清单 TSV | off |
| `--emit-samplesheet-template` | 生成候选样本表 CSV | off |
| `--strandedness-default` | 候选样本表的 strandedness 默认值 | auto |
| `--email/--api-key` | NCBI 请求可选凭据（避免限流） | None |

**输出**：`{prefix}.metadata_curated.tsv`（主表，28 列）、`{prefix}.metadata_raw.{tsv,json}`、`{prefix}.downloads.tsv`、`{prefix}.samplesheet_from_metadata.csv`、`{prefix}.warnings.tsv`。

**metadata_curated.tsv 关键列**（审阅/标注就看这几列）：

```
input_accession, run_accession, sample_title, library_strategy,
fastq_ftp_1/2, fastq_md5_1/2, sra_download_accession,
inferred_type, type_confidence, type_evidence,
suggested_group, group_evidence, replicate_hint, needs_manual_review
```

`inferred_type`（riboseq/rnaseq/tiseq/unknown）与 `suggested_group` 是启发式推断，
**needs_manual_review=true 的 run 必须人工确认**——这一步不能跳过。

## 样本标注（两种方式）

**方式 A：notebook 交互式**（用户手动，推荐给新用户）：

```bash
cd docs/notebooks && jupyter lab riboseq_guided_tutorial.ipynb
# 或 VS Code 打开，内核选 Python (riboseq-notebook)，环境见 docs/notebooks/README.md
```

Step 0 输入 accession（文本框支持 GSE 解析）→ 元数据表格 + 每 run 的
type Dropdown / group Text → "Apply review edits" → "Export review CSV"。
**设计约束**：存在 `inferred_type == "unknown"` 时 `build_samplesheet()` 会
抛错阻断，保证先审阅后生成。

**方式 B：CLI + 人工审阅**（Claude 代跑）：跑完 fetch 后展示
`metadata_curated.tsv` 的审阅列，由用户（或按用户指令）修正
inferred_type/suggested_group，再进入下载和样本表步骤。

## 下载命令生成

metadata 里 `fastq_ftp_1/2` 存在 → **ena-ftp**（推荐，直接下 FASTQ）；
缺失 → **ncbi-sra**（prefetch + fasterq-dump 兜底，需 SRA Toolkit）。

```bash
# ena-ftp（三选一；ascp 最快但需 Aspera 客户端）
wget -c ftp://{fastq_ftp_1} -O data/{文件名}
curl -L ftp://{fastq_ftp_1} -o data/{文件名}
ascp -QT -l 300m -P33001 era-fasp@fasp.sra.ebi.ac.uk:{fastq_ftp_1} data/

# ncbi-sra fallback
prefetch {acc} && fasterq-dump --split-files --threads 8 {acc} && pigz -p 8 {acc}*.fastq
```

下载后验证 `fastq_md5_1/2`（md5sum 比对）。已下载为 .sra 的本地文件用
`scripts/sra2fq.sh` 转换（见 SKILL.md §5）。

## 样本表生成规则（与 notebook 一致）

- `sample` 名取自 `sample_title`（空则用 run_accession），非字母数字替换为 `_`，重名加 `_2/_3` 后缀
- 列：`sample, fastq_1, fastq_2, strandedness, type, group, run_accession, input_accession`
- `fastq_1` 填 `fastq_ftp_1` 的文件名（下载后的本地名）或 SRA accession
- type 用审阅后的 inferred_type；strandedness 默认 auto（BAM 输入模式不能 auto，见 riboseq-pipeline-run）

## notebook 依赖与运行环境

- 依赖：pandas、ipywidgets、matplotlib（缺一个 `launch_tutorial()` 会报缺依赖清单）
- 需在仓库根内打开（helpers 靠 repo root 定位 `fetch_public_metadata.py`）
- 教程输出默认写到 `tutorial_outputs/`；演示数据在 `test_data/tutorial_demo_public_data/`
  （含 accessions.txt、metadata_curated.tsv、samplesheet_from_metadata.csv 等成品示例）

## 已知坑

1. **unknown type 阻断样本表**：CLI 生成的候选样本表不会阻断（模板含 unknown），
   但 notebook 的 build_samplesheet 会抛错。两条路都要保证 unknown 清零后再跑 pipeline。
2. **ENA FTP 偶发不可达**：个别 run 的 FTP URL 失效时用 ncbi-sra fallback
   （manifest 里每个文件都带 `fallback_sra_accession`）。
3. **NCBI 限流**：大量查询时传 `--email/--api-key`。
4. **type=tiseq 是空壳分支**（见 SKILL.md §2）：TI-seq 数据标 `riboseq` 才会走完整流程。
