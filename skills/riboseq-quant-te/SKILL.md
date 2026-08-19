---
name: riboseq-quant-te
description: >
  ORF 定量（RPKM/TPM、P-site 计数）、post_analysis 两阶段过滤与 ggRibo 绘图、
  以及 TE/ΔTE（翻译效率）分析流程。
  涵盖 pipeline 产出的表达量文件（unified_orfs_expression_summary/rpkm_tpm.tsv 列与公式）、
  P-site 数据语义（_P_sites_ 整数 bedgraph vs _coverage_ RPM 浮点、pN 一词三义）、
  scripts/post_analysis/run_all.sh 四步流程（初筛→P-site 纯度→过滤→ggRibo 批量图）、
  单基因 ggRibo 绘图（bin/plot_orf_ggribo.R + 最小 GTF + 隐式全局变量坑）、
  ΔTE 模块（DESeq2 交互模型 condition+type+condition:type、--contrasts、预过滤）。
  当用户提到 ORF 定量/表达量/RPKM/TPM、pN、P-site 纯度/过滤、ggRibo/三框图/翻译图、
  TE/ΔTE/翻译效率/DESeq2/差异翻译、post_analysis、final_orfs_for_ggribo 时使用本 skill。
---

# riboseq 定量与 TE 分析（quantification → post_analysis → ggRibo → ΔTE）

## 决策树

- **要理解表达量文件（列含义/公式）** → 第 1-2 节 → `references/quantification.md`
- **要跑两阶段过滤 + 批量 ggRibo 图（post_analysis 流程）** → 第 3 节 → `references/post_analysis_workflow.md`
- **要手动定量（pipeline 之外）** → 第 4 节 → `references/quantification.md`
- **要做 TE/ΔTE（差异翻译）** → 第 5 节 → `references/te_analysis.md`
- **要画单个 ORF 的 ggRibo 图** → 第 6 节 → `references/ggribo_plotting.md`
- **pN / p_site_pct / reads 语义不清** → 第 7 节（数据语义速查）

## 1. pipeline 产出的表达量文件（UNIFY 发布）

位置 `{result}/orf_unification/`（输入是 RiboseQC 的 coverage/P-sites bedgraph）：

| 文件 | 列 |
|---|---|
| `unified_orfs_expression_summary.tsv` | `orf_id, chrom, start, end, strand` + 每样本 `{s}_reads`、`{s}_pN` + `total_reads, n_expressed_samples` |
| `unified_orfs_expression_rpkm_tpm.tsv` | `orf_id, ..., orf_length, orf_length_kb` + 每样本 `{s}_coverage_rpm, {s}_rpkm, {s}_tpm` |

**公式**（coverage bedgraph 值已是 RPM，无需 library size）：
- `rpkm = coverage_rpm / orf_length_kb`
- `tpm = rpkm / Σ(该样本所有 ORF 的 rpkm) × 1e6`
- 每样本 `pN = max_psite × count_at_max / Σ(psite)`（P-site 峰/均值比，4 位小数）

⚠️ **`{s}_reads` 与真实 P-site 计数不是一回事**：`reads` 来自 coverage bedgraph
（RPKM 换算回整数，与 `p_site_GSE` 相关性 r=0.99），`pN` 是周期性分数（与
`p_site_pct` 相关性仅 r=-0.20）。两阶段过滤 Stage 1 用 reads+pN 初筛，Stage 2 才用
真实 P-site。

## 2. P-site 数据语义（RiboseQC 产物）

`{result}/riboseqc/` 下每样本 8 类文件（plus/minus × P_sites/coverage × 普通/uniq）：

- `{s}_P_sites_{plus,minus}.bedgraph`：value = **整数 P-site 计数**
- `{s}_coverage_{plus,minus}.bedgraph`：value = **RPM 浮点数**（不能当整数读！）
- `{s}_P_sites_calcs`：RiboseQC 自己的 offset/周期性表（frame_preference 等）
- `{s}_ggribo.tsv`：pipeline 自动生成的 ggRibo 4 列表
  （`count \t chrom \t position(1-based) \t strand`，由 P_sites bedgraph awk 转换）

**`p_site_pct` 口径**：`= p_site_GSE / expression_summary 的 {s}_reads`（两个整数同单位，
在 `_functions.R` 里重算）——**不要用** compute_psite_*.py 输出的 `reads_GSE/pct` 列
（那是 coverage bedgraph 求和，RPM 单位，做分母无意义，已弃用）。

## 3. post_analysis 两阶段流程（批量 ggRibo 的标准入口）

```bash
bash scripts/post_analysis/run_all.sh <project_config.yaml>
# 模板：scripts/post_analysis/config_template.yaml（复制后填 {result_dir}）
# 真实例子：run/rice/post_analysis/project_config_with_orfquant.yaml
```

四步：Step1 初筛（quarto render 01_prelim_analysis.qmd → `prelim_orfs_for_psite.bed`）
→ Step2 P-site 纯度（`compute_psite_purity.py`，纯 Python 精确回溯，**无需 bedtools**，
2026-08-19 起为 run_all.sh 默认）→ Step3 P-site 过滤（02_psite_filtering.qmd →
`final_orfs_for_ggribo.tsv`）→ Step4 ggRibo 批量图（03_generate_ggribo.qmd →
`ggribo_plots/{biotype}/{orf_id}.png`）。

**两级阈值**（config `filtering:` 段，均 per-sample 且需 ≥ N 样本同时满足）：
Stage1 `prelim_reads_min: 9`、`prelim_pN_min: 0.5`、`prelim_cross_sample: 1`；
Stage2 `p_site_gse_min: 9`、`p_site_pct_min: 0.5`、`p_site_pos_min: 2`、
`final_cross_sample: 2`；最后再叠加 `length_aa ≤ 50` 与 biotype 白名单。

**rice 实测**（output_with_orfquant_strandfix，final_cross_sample=1）：
394,053 unified → 207,413 prelim → **10,255 final**（ggribo_plots/ 10,255 张 PNG）。
⚠️ 目录内 README 记录的旧数字（86,375 → 35,073 → 5,843 → 1,590）已过时。

注意：`compute_psite_fast.py` 的 `p_site_pos` 恒为 0.00（只适合快速预检）；需要真实
位置信息用 `compute_psite_purity.py`。2026-08-19 起 purity.py 的 overlap 查找已 numpy
向量化（~20-50x，strand 过滤保留、语义与旧循环完全一致）——旧基准 rice 3.1h（394K
ORF×23，向量化之前）已大幅过时；本机 rice（394K ORF）与 maize（660K ORF×97 样本）
现有 psite_purity.tsv 均已是 purity.py 产出，run_all.sh 的 fast→purity 切换对它们
重跑无行为变化。

## 4. 手动定量（pipeline 之外）

- **`scripts/R/quantify_orfs_from_psites.R`**（optparse CLI）：`--gtf` + `--bedgraph-dir`
  + `--sample-pattern`（默认 `(.+)_P_sites_unique_(plus|minus)\.bedgraph$`）+
  `--min-count 5` → 输出 `orf_counts_raw.tsv / orf_counts_tpm.tsv /
  orf_counts_matrix.csv（for DESeq2）/ sample_summary_stats.tsv / orf_summary_stats.tsv`
- **`scripts/quant_analysis.R`**：CHX 小鼠单对比模板（DESeq2 `~ condition`，非 CLI，
  样本注释硬编码在脚本顶部）——不是通用工具，抄改时改 L26-38 与 bedgraph 命名
- **`bin/calc_orf_rpkm_tpm.py`**：`--expression expression_summary.tsv --coverage-dir
  riboseqc/ --sample-pattern "*_coverage_plus.bedgraph" --output ...`（独立重算 RPKM/TPM）

细节 → `references/quantification.md`。

## 5. TE/ΔTE（pipeline 内，⚠️ 从未在真实项目启用过）

开关 `--skip_te_analysis`（默认 false）；子流程 `subworkflows/local/te_analysis.nf`：
`QUANTIFY_ORFS`（featureCounts，BED12→SAF）→ `MERGE_COUNTS`（merged_counts.tsv +
sample_sheet.csv）→ `DESEQ2_DELTATE`。

**启用前提**：samplesheet 必须有 `type` 之外的 treatment/group 列（如 `condition`），
并传 `--contrasts contrasts.csv`（列：`id, variable, reference, target` [+ `batch`]；
`variable`/`batch` 必须是样本表列名）。

**核心模型**：DESeq2 交互设计 `~ condition + type + condition:type`
（Chothani 2019 deltaTE）；`type` 列自动识别 ribo/rna
（grep `ribo|rp|fp` / `rna|mrna|total|lncrna`）；`sfType="poscounts"`；
交互项系数即 TE。

**关键参数**（nextflow.config L196-215）：
`te_lfc_threshold 0.2630344`、`te_prefilter_min_nonzero 2`、`te_prefilter_min_frac 0.2`
（预过滤防 `estimateSizeFactors: every gene contains at least one zero` 崩溃）、
`extra_deltate_args null`（`"--shrinkage_type normal --shrink_lfc true"` 等）、
`deseq2_container`（默认 wave 拉取，可 `apptainer build deseq2_deltate.sif
containers/Singularity.r_te_analysis.def`）。

输出 `{result}/translational_efficiency/deltate/`：`*.translation.deltate.results.tsv`
（TE）、`*.translated_mRNA...`（Ribo）、`*.total_mRNA...`（RNA）、anota2seq 风格分类
基因表（intensified/buffering/translation/mRNA_abundance）、PCA/heatmap/volcano 图。

细节与坑 → `references/te_analysis.md`。

## 6. 单个 ORF 的 ggRibo 图

**数据准备**（二选一）：
1. pipeline 产出的 `{sample}_ggribo.tsv`（推荐，RiboseQC 模块自动生成）
2. 直接给 bedgraph：`list(plus="..._P_sites_plus.bedgraph", minus="..._P_sites_minus.bedgraph")`

**批量**：跑第 3 节的 run_all.sh（自动），或 `run/maize/scripts/run_ggribo_batch.sh`
（每批 50 ORF、top 10 样本；maize 实测 13,959 张 PNG）、
`run/rice/scripts/run_ggribo_parallel.sh`（GNU parallel 12 workers）。

**单基因**：`bin/plot_orf_ggribo.R`（optparse CLI，与 maize/rice 版同源）：
`--orf-meta unified_orfs.metadata.tsv --expression expression_summary.tsv
--psites-dir riboseqc/ --gtf unified_orfs.gtf --orf-ids ORF_1_gene1
--output-dir plots --extend 200 --n-samples-per-orf 3`。

**三个必踩的坑**（详见 `references/ggribo_plotting.md`）：
1. ggRibo 函数隐式依赖全局变量 → 先 `assign("Txome_Range", ..., .GlobalEnv)` 和
   `assign("inputs_full", ..., .GlobalEnv)`
2. unified_orfs.gtf 缺 gene/transcript 行，TxDb 解析失败 → 每 ORF 生成最小 GTF
   （1 gene + 1 transcript + CDS/exon，多外显子用 exon_blocks 展开）
3. BiocGenerics 命名空间冲突（与 ORFquant 同源）→ 不要调 ggRibo 的 `gtf_import`，
   用自定义 `gtf_import_custom`（txdbmaker::makeTxDbFromGFF + 最小 Range_info 类）

## 7. 数据语义速查（pN 一词三义）

| 出现位置 | 含义 | 公式 |
|---|---|---|
| metadata.tsv 聚合 `pN` | P-site 密度/nt | `total_psites / length_nt`（`unique_pN = unique_psites / length_nt`） |
| expression_summary 每样本 `{s}_pN` | 周期性分数（峰/均值比） | `max × count_at_max / Σ(psite)` |
| quant_analysis.R 每样本 `pN` | 样本级 P-site 密度/nt | `unique_psites / length_nt`（RiboseQC unique bedgraph 统计） |

两阶段过滤 Stage1 的 `prelim_pN_min: 0.5` 用的是第 2 义；ggRibo 标题里的 pN 也是第 2 义。
旧 strict sORF 流程的 `max_pn ≥ 5` 是同一列但阈值量级完全不同（旧体系）。

## 完成后

- 上游输入（unified/分类/riboseqc 目录）→ riboseq-orf-analysis skill
- 四级评分框架（S_purity/S_expression/S_confidence/S_biology 加权）→
  `~/riboseq/post_analysis/design_integrated_scoring.md` + `integrate_orf_scores.R`
