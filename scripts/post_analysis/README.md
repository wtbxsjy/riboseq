# ORF Post-Analysis Module

Reusable analysis pipeline for nf-core/riboseq ORF outputs.

## Quick Start

```bash
# 1. Copy config template to your project
cp scripts/post_analysis/config_template.yaml run/YOUR_PROJECT/post_analysis/project_config.yaml

# 2. Edit the config with your paths

# 3. Run
bash scripts/post_analysis/run_all.sh run/YOUR_PROJECT/post_analysis/project_config.yaml
```

## Input Requirements

All files are standard pipeline outputs:

| File | Source Module | Description |
|------|--------------|-------------|
| `unified_orfs.metadata.tsv` | UNIFY_ORF_PREDICTIONS | ORF coordinates + metadata |
| `unified_orfs_orf_confidence.tsv` | ORF_QC | OCS scores + tier |
| `unified_orfs_expression_summary.tsv` | EXPRESSION_QUANT | Per-sample reads + pN |
| `gencode_results.orfs.out.gz` | CLASSIFY_ORFS_GENCODE | ORF biotype classification |
| `unified_orfs.bed.gz` | UNIFY_ORF_PREDICTIONS | BED12 ORF coordinates |
| `riboseqc/*_P_sites_{plus,minus}.bedgraph` | RiboseQC | P-site density |

## Two-Stage Filtering Pipeline

The pipeline uses a **two-stage** strategy to balance accuracy and performance:

### Stage 1: Preliminary (expression-based, fast)
`01_prelim_analysis.qmd` — Uses `{sample}_reads` and `{sample}_pN` from the expression
summary to eliminate ORFs with insufficient signal. Reads correlate strongly with
real P-site counts (r=0.99); pN measures frame periodicity.

**Output:** `prelim_orfs_for_psite.bed` — reduced ORF set for P-site computation.

### Stage 2: Real P-site (bedgraph-based, exact)
`compute_psite_purity.py` — Backtracks through actual RiboseQC bedgraph files
(`_P_sites_{plus,minus}.bedgraph` and `_coverage_{plus,minus}.bedgraph`),
computing exact `p_site_GSE`, `p_site_pct`, and `p_site_pos` per ORF per sample.
Runs on the Stage 1 reduced set, dramatically cutting runtime.

`02_psite_filtering.qmd` — Applies real P-site thresholds and selects final candidates.

### Stage 3: ggRibo Visualization
`03_generate_ggribo.qmd` — Generates per-ORF ribosome footprint coverage plots
with reading-frame shading. Top N samples sorted by a configurable metric
(default: `p_site_GSE`).

## Output

```
{output_dir}/
├── 01_prelim_analysis.html       # Stage 1: preliminary analysis
├── prelim_orfs_for_psite.bed      # ORFs for P-site computation
├── psite_purity.tsv               # Real P-site purity data
├── 02_psite_filtering.html        # Stage 2: P-site filtering report
├── final_orfs_for_ggribo.tsv      # Final selected ORFs
├── 03_generate_ggribo.html        # ggRibo plots + index
├── ggribo_plots/                  # Per-biotype ggRibo PNGs
│   ├── uORF/
│   ├── dORF/
│   ├── intergenic/
│   └── ...
├── tmp_gtf/                       # Temp GTF files (can delete)
└── logs/
```

## Intergenic ORFs

ORFs that don't overlap any annotated transcript in the GENCODE reference are
labelled **"intergenic"**. They are inherently non-CDS and are included in
downstream analysis alongside classified non-CDS ORFs (uORF, dORF, etc.).
This recovers the ~53% of unified ORFs that would otherwise be silently dropped.

## Adding a New Project

```bash
mkdir -p run/NEW_PROJECT/post_analysis/output/{logs,ggribo_plots}
cp config_template.yaml run/NEW_PROJECT/post_analysis/project_config.yaml
# Edit paths in project_config.yaml
bash run_all.sh run/NEW_PROJECT/post_analysis/project_config.yaml
```

---

## Manual two-stage chain（PRJEB26593 / GSE120762 / rice / maize 实际使用的路径）

`run_all.sh` 是 Quarto 打包版；PRJEB26593 与 mouse_GSE120762 实际走的是下面这条**脚本链**，
判据相同但完全可控、可复现、可中途检查。新项目建议按这个顺序做：

```bash
P=<project>                       # e.g. mouse_GSE120762
R=run/$P/result
O=post_analysis/$P                # 输出目录

# 1. CDS 排除 + Stage-1（reads > 9，≥1 样本）
python3 scripts/post_analysis/stage1_expression_filter.py \
  --orftype  $R/orf_classification/orf_type/orftype_classification.tsv \
  --gencode  $R/orf_classification/gencode/gencode_results.orfs.out.gz \
  --expression $R/orf_unification/unified_orfs_expression_summary.tsv \
  --out-dir $O --prefix $P

# 2. 从 unified BED 抽出 stage-1 子集（**加速关键**：purity 只扫这一步的 ORF）
awk -F'\t' 'NR>1{print $1}' $O/${P}_stage1_passed.tsv | sort > $O/stage1_ids.txt
zcat $R/orf_unification/unified_orfs.bed.gz \
  | awk -F'\t' 'NR==FNR{k[$1];next} ($4 in k)' $O/stage1_ids.txt - > $O/stage1_orfs.bed

# 3. P-site 纯度（务必用仓库版本，见下方坑 1）
python3 scripts/post_analysis/compute_psite_purity.py \
  --bed $O/stage1_orfs.bed --riboseqc-dir $R/riboseqc \
  --output $O/${P}_psite_purity.tsv --workers 12

# 4. Stage-2 过滤 + 长度 + 每 biotype top 10%
python3 scripts/post_analysis/stage2_psite_filter.py \
  --purity $O/${P}_psite_purity.tsv --stage1 $O/${P}_stage1_passed.tsv \
  --noncds $O/step1_nonCDS.tsv.gz --expression $R/orf_unification/unified_orfs_expression_summary.tsv \
  --out-dir $O --prefix $P

# 5. FASTA 两套（每个集合各一份原始序列 + 严格 CDS 的 NA/AA）
python3 scripts/post_analysis/extract_orfs_fasta.py \
  --metadata $R/orf_unification/unified_orfs.metadata.tsv \
  --ids $O/${P}_stage2_passed.tsv --out $O/${P}_stage2_orfs.fa
#     ... 同样对 ${P}_filtered_len50aa.tsv 跑一次
#     CDS 版：extract_cds_fasta.py（PRJEB26593 版按项目拷一份改路径即可）
```

**加速**：第 3 步的 purity 是唯一的重计算（PRJEB26593 全量 120 万 ORF 要 9.5 h）。只跑
stage-1 子集（约 5–15 万 ORF）可降到 15–55 分钟，**结果完全等价** ——
因为 stage2 只消费 stage1 的子集，子集外的纯度值无人读取。

### 已知坑（都会静默给出错误数字，务必对照）

1. **`compute_psite_purity.py` 必须用仓库版本。** 仓库版含
   `find_overlapping_orfs(..., strand=)`：bedgraph 的 P-site 只匹配**同链** ORF。
   曾经的工作副本缺这个参数（见 commit `cd60f79`），会把正链 bedgraph 的 P-site
   同时记到反链 ORF 上，在重叠密集区显著高估计数 —— 而且**不会报错**。
   拿不准时：`grep -n "strand=" compute_psite_purity.py` 应能看到调用点带 `strand_char`。
2. **`extract_cds_fasta.py` 的终止密码子剥离要求长度已按密码子对齐**
   （CLAUDE.md gotcha 31）。原始 `sequence` 带 1–2 nt 边界噪声，`len(seq) % 3 != 0`
   是常态；对一条 37 nt / 12 aa 的 ORF，末 3 字符 `TAA` **不是密码子**，误删会让
   NA 短 3 nt。判据：`len(NA) == 3*len(AA)` 必须 100% 成立，任何 mismatch 都要查。
   该 bug 会被 PRICE 来源的 ORF 集中触发（修复前 PRICE 的 `sequence` 全是 `N`，
   所以旧结果 mismatch 恒为 0 —— 那是假象，不是更干净）。
3. **不要用 stage2 的 `total_psites` 与表达量做比值**：purity 的 `_p_site_pct` 分母是
   RPM 归一的 coverage bedgraph，量纲不匹配；正确比值是
   `{sample}_p_site_GSE / {sample}_reads`（见 CLAUDE.md gotcha 28）。
4. **pN 不能当阈值**：结构性 ≥ 1，`pN > 0.5` 是空判据（gotcha 28）。
   Stage-1 的有效门就是 `reads > 9`。
