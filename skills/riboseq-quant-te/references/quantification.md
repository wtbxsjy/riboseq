# 定量：表达量文件、公式与手动脚本

## pipeline 产出（UNIFY 发布，位置 `{result}/orf_unification/`）

**`unified_orfs_expression_summary.tsv`**：
`orf_id, chrom, start, end, strand` + 每样本 `{s}_reads`、`{s}_pN` +
`total_reads, n_expressed_samples`

- `{s}_reads`：coverage bedgraph（RPM 浮点）换算的整数 reads（≈ RPKM×bp），
  与真实 P-site 计数 `p_site_GSE` 相关 r=0.99 —— Stage1 初筛用它
- `{s}_pN`：P-site 周期性分数（峰/均值比），不是 P-site 比例的代理
  （与 p_site_pct 相关 r=-0.20）——Stage1 用 `pN > 0.5` 过滤

**`unified_orfs_expression_rpkm_tpm.tsv`**：
`orf_id, chrom, start, end, strand, orf_length, orf_length_kb` +
每样本 `{s}_coverage_rpm, {s}_rpkm, {s}_tpm`

## 公式（两处实现一致，`scripts/unify_orf_predictions.py` L1993-2013 与 `bin/calc_orf_rpkm_tpm.py`）

```python
length_kb = length_nt / 1000.0
rpkm = coverage_rpm / length_kb          # coverage bedgraph 已是 RPM，无需 library size
tpm  = rpkm / Σ(rpkm per sample) * 1e6
pN   = max_psite * count_at_max / Σ(psite)   # 每样本周期性分数
```

## P-site bedgraph 语义（riboseqc 目录）

- `{s}_P_sites_{plus,minus}.bedgraph`：value = **整数** P-site 计数
- `{s}_coverage_{plus,minus}.bedgraph`：value = **RPM 浮点**（0.00502 这种）
- `{s}_P_sites_uniq_*` / `{s}_coverage_uniq_*`：unique（去重）版本
- `{s}_P_sites_calcs`：RiboseQC offset/周期性表
- `{s}_ggribo.tsv`：4 列无 header（count/chrom/position(1-based)/strand），
  由 P_sites bedgraph awk 转换（bedgraph 0-based start → ggRibo 要 1-based）

## pN 一词三义（解析不同文件时先确认是哪个）

| # | 位置 | 含义 | 公式 |
|---|---|---|---|
| 1 | `unified_orfs.metadata.tsv` 的 `pN`/`unique_pN` | P-site 密度/nt | `total_psites / length_nt` |
| 2 | `expression_summary.tsv` 的 `{s}_pN` | 周期性分数 | `max_psite × count_at_max / Σ(psite)` |
| 3 | `quant_analysis.R` 样本级 pN | unique P-site 密度/nt | `unique_psites / length_nt` |

RiboseQC 自身输出**没有**叫 pN 的字段（周期性在 `_P_sites_calcs` 的
frame_preference/gain_codons）。

## 手动定量脚本

### scripts/R/quantify_orfs_from_psites.R（optparse CLI，最通用的独立工具）

```bash
Rscript scripts/R/quantify_orfs_from_psites.R \
  --gtf results/Mouse_Unified.orfs.gtf \
  --bedgraph-dir riboseqc_results/ \
  --sample-pattern "(.+)_P_sites_unique_(plus|minus)\\.bedgraph$" \
  --min-count 5 \
  --outdir quantification_results
```

- 输入：GENCODE ORF GTF（gencode-riboseqORFs 产出）+ RiboseQC unique P-site bedgraph
- 内部：regex 配对 plus/minus 样本 → strand-specific `findOverlaps(type="within")`
  计数（bedgraph 0-based→1-based `start+1`）→ `total_counts >= --min-count` 过滤 →
  TPM 归一
- 输出：`orf_counts_raw.tsv`（计数+注释列）、`orf_counts_tpm.tsv`、
  `orf_counts_matrix.csv`（纯计数矩阵，注释说 for DESeq2）、
  `sample_summary_stats.tsv`、`orf_summary_stats.tsv`、`session_info.txt`
- 注意：`--threads` 参数存在但循环实际串行

### scripts/quant_analysis.R（CHX 小鼠一次性模板，非 CLI）

- 顶部硬编码（L26-38）：metadata 路径、riboseqc 目录、outdir、阈值
  （padj<0.05、|lfc|≥1）、4 样本 2 条件的 `sample_annotations` data.table
- 计数：`data.table::foverlaps` 把 `{s}_P_sites_uniq_{plus|minus}.bedgraph`
  按 strand 重叠到 ORF exon blocks，`value × overlap_width` 加权求和
- DESeq2：`design = ~ condition`，`results(dds, contrast=c("condition","CHX_LPS","CHX_NT"))`
  —— 本脚本**无 lfcShrink**
- 预过滤：`length_aa >= 16` 且聚合 `pN >= 1`，重算 per-sample pN 后
  `max_sample_pN >= 1`
- 输出：long/wide 矩阵、`deseq2_*.orf_results.tsv`（含 significant/direction 列）、
  volcano.png、ComplexHeatmap 热图
- 抄改要点：改顶部路径与样本注释，确认 bedgraph 命名模式

### bin/calc_orf_rpkm_tpm.py（独立重算 RPKM/TPM）

```bash
python3 bin/calc_orf_rpkm_tpm.py \
  --expression expression_summary.tsv \
  --coverage-dir riboseqc/ \
  --sample-pattern "*_coverage_plus.bedgraph" \
  --output expression_rpkm_tpm.tsv --workers 4
```

### 历史产物（~/riboseq/post_analysis/ 旧体系，谨慎混用）

`rice_expression_summary.tsv`（16,668 ORF）、`rice_rpkm_tpm.tsv`、
`rice_orf_bedgraph_all.tsv`（473MB）等是旧版脚本（extract_amp_expression*.py、
calc_orf_rpkm*.py）产物，列定义与新 pipeline 不完全一致，不要和新 unified 文件交叉使用。
