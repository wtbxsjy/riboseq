# 各阶段结果检查清单

用自带脚本一键体检：`scripts/check_pipeline_outputs.sh <result_dir>`。
以下按 result/ 目录结构逐阶段说明"正常输出长什么样"和"异常意味着什么"。
基准数字来自 rice（23 样本，369,939 unified ORFs）和 maize（97 样本，~660K）实测。

## 逐阶段检查表

| 阶段 | 关键文件 | 正常形态 | 异常信号与解读 |
|---|---|---|---|
| alignment | `alignment/star/*.bam(.bai)` | 每样本 1 对 BAM/BAI | 缺样本 → 比对失败或样本被过滤；查 .nextflow.log |
| genome | `genome/index/`、`*.fa.fai`、`*.fa.sizes`、`*.filtered.gtf` | `--save_reference true` 时存在 | 缺 filtered.gtf → 下游 unify/classify 用的是原始 GTF |
| riboseqc | `{sample}_P_sites_{plus,minus}.bedgraph`、`_coverage_*`、`_for_ORFquant`、`_P_sites_calcs` | 每样本 8-9 个文件 | 某样本全缺 → RIBOSEQC 失败（常见原因：低质量样本，salmon 比对率 <20%）；`_for_ORFquant` 缺 → ORFquant 自动跳过 |
| orf_predictions/ribotish | `{sample}_pred.txt` | rice 每样本 5.7万–12.4万 行 | 行数骤减 → 深度不足或过滤参数变了 |
| orf_predictions/ribotricer | `{sample}_translating_ORFs.tsv` | 行数显著少于 ribotish | 全空 → phase_score_cutoff 过严（植物用 0.1） |
| orf_predictions/orfquant | `{sample}_Detected_ORFs.gtf.gz` | 每样本 1 个；5-11 分钟/样本 | 全缺 → 跳过或 BSgenome 问题（非模式生物 saga，见 riboseq-orf-analysis） |
| orf_predictions/price | `{sample}.orfs.tsv` | 修复后 ~18,091 ORFs/样本 | 行数异常 → GEDI 参数问题（-genomic 需绝对路径） |
| orf_predictions/ribocode | `{sample}_collapsed.gtf.gz` 等 7 类 | 仅 `--run_ribocode` 时存在 | 缺 → 未启用或转录组 BAM 缺失 |
| orf_unification | `unified_orfs.bed.gz/.gtf.gz/metadata.tsv/stats.txt` + `*_expression_summary.tsv` + `*_expression_rpkm_tpm.tsv` | metadata 行数：rice 369,939；summary 行数 == metadata 行数 | summary 行数 ≠ metadata 行数 → EXPRESSION_QUANT 未完成或部分失败 |
| orf_unification | `unified_orfs.stats.txt` | 各工具贡献计数（如 Ribo-TISH 767,610 / ORFquant 181,630 / Ribotricer 37,561 为 human 62 样本） | 某工具计数为 0 → 该工具结果没进 unify（检查输入路径/后缀解析） |
| orf | `unified_orfs_orf_confidence.tsv`、`_psite_harmonized.tsv`、`_tool_agreement.tsv`、`_qc_report.html` | 行数 ≈ metadata 行数 | 缺 → ORF_QC 被 skip 或 errorStrategy='ignore' 静默失败 |
| orf_classification/gencode | `gencode_results.orfs.out.gz`（实际为 gz 压缩，human 会话里有未压缩 .out） | biotype 分布正常（见下） | `lncRNA > 90%` → gffread 蛋白 header 不匹配坑（660K→2190 案例）；BED 未排序 → 分类错乱 |
| orf_classification/orf_type | `orftype_classification.tsv` | 行数 ≥ unified ORF 数（含 frame 变体） | 缺 → orf_classify_mode 被设为非 orf_type |
| expression | `expression_quant_expression_{summary,rpkm_tpm}.tsv` | 与 orf_unification 下同名文件内容一致 | 文件极小 → 量化脚本早期版本（awk O(n×m) 27h 版）或未完成 |
| multiqc | MultiQC 报告 | 存在 | — |

## biotype 分布检查（GENCODE 分类后必做）

```bash
# ⚠️ orf_biotype 是第 10 列（文件有 header 行，第 7 列是 trans，cut -f7 会得到转录本 ID）
# 实际输出多为 .out.gz；未压缩的 .out 把 zcat 换成 cat 即可
zcat gencode_results.orfs.out.gz | tail -n +2 | cut -f10 | sort | uniq -c
```

正常分布：CDS/dORF/uORF/doORF/uoORF/intORF/lncRNA 各类均有，非 CDS 类占大头
（unified ORFs 本身就是找非常规 ORF 的）。异常：**lncRNA 占比 >90%** =
蛋白 FASTA header 用 transcript_id 而非 protein_id（详见
riboseq-data-prep skill 的 ensembl_dir_manual.md 坑 1），或 **2 类塌缩**
（只剩 lncRNA/CDS）= FASTA header 版本后缀与 GTF protein_id 不匹配，2026-08-19
已修复（commit 7281ecd，双 key 索引），见 riboseq-orf-analysis classifiers.md。

## 快速健康数字（跨项目参考）

| 项目 | 样本数 | unified ORFs | metadata.tsv | summary.tsv | rpkm_tpm.tsv |
|---|---|---|---|---|---|
| rice | 23 | 369,939（最终） | 271MB | 47MB | 89MB |
| maize | 97 | ~660K | — | — | — |
| human_GSE157490 | 62 | 956,971 | — | — | — |

行数差几个数量级（如 rice 只有几百个 ORF）说明 unify 输入不全或
`--unify_orf_min_len` 设置错误（wheat/soybean 用 24，rice 用 6）。
