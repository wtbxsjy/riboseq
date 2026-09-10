# unify_orf_predictions 输出数据模型

脚本：`scripts/unify_orf_predictions.py`（pipeline 内 UNIFY_ORF_PREDICTIONS 模块，
容器 unify_orf.sif ~190MB，默认 fallback biopython:1.79）。

## 输出文件（publishDir: result/orf_unification/）

| 文件 | 内容 |
|---|---|
| `unified_orfs.bed.gz` | BED12 坐标 |
| `unified_orfs.gtf.gz` | GTF 注释 |
| `unified_orfs.metadata.tsv` | 主表（26 列，见下） |
| `unified_orfs.stats.txt` | 各工具输入统计（行数、过滤数） |
| `unified_orfs_expression_summary.tsv` | per-sample reads + pN（EXPRESSION_QUANT 合并进 UNIFY 后产出） |
| `unified_orfs_expression_rpkm_tpm.tsv` | per-sample RPM/RPKM/TPM |
| `versions.yml` | 版本信息 |

## metadata.tsv 列（26 列，代码原文）

```
orf_id, chrom, strand, start, end, length_aa, exon_blocks,
gene_id, transcript_id, tools, samples, tool_scores, tool_pvalues,
total_reads, unique_reads, total_psites, unique_psites, pN, unique_pN,
num_subset_orfs, subset_orfs, sequence, start_codon, aa_sequence,
is_cds_overlap, overlapping_genes
```

- `orf_id = ORF_{i}_{gid}`（全局递增序号 + 基因 id）
- `tools` 列是逗号分隔的**大写工具名**：`Ribo-TISH` / `Ribotricer` / `ORFquant` /
  `PRICE` / `RiboCode`——过滤、计数、解析 sources 集合时必须用这些精确字符串
  （不是小写）
- `samples` 列逗号分隔样本 id
- `exon_blocks` 格式：`start-end,start-end`（解析多外显子用）

## 合并规则（三步，代码证实）

1. **exact-match**：坐标完全相同 → 去重
2. **frame-aware**：同读框 ORF 合并（human 实测贡献 18.3% 缩减）
3. **overlap grouping**：重叠聚类，每组选一个代表（`subset_orfs` 记录被合并成员）

## 各工具解析行为

| 工具 | 输入文件 | 分数来源 | 过滤 |
|---|---|---|---|
| Ribo-TISH | `{sample}_pred.txt`（每样本 5.7万–12.4万 行，human 单文件可达 160 万行） | score + pvalue | TisType / start codon（unify 日志：`Filtered: N by TisType`） |
| Ribotricer | `{sample}_translating_ORFs.tsv` | score | — |
| ORFquant | `{sample}_Detected_ORFs.gtf.gz` | GTF 属性 `P_sites`（整数） | — |
| PRICE | `{sample}.orfs.tsv` | score_from_pvalue | — |
| RiboCode | `{sample}_collapsed.gtf.gz`（优先解析 collapsed，后缀检测顺序） | score_from_pvalue | — |

- 各工具按文件后缀推断样本归属（`infer_sample_id_from_prediction_path`）——
  文件名必须符合 `{sample}_xxx` 约定
- 旧坑（已修）：`.gtf.gz` 的 `endswith('.gtf')` 判断漏匹配；bash 先 gunzip 后
  Python 仍按 `.gz` 路径读

**contig 命名归一化（2026-09-10 新增）**：PRICE 用 Ensembl 风格 contig 名
（`1` … `X`,`Y`,`MT`），其余工具与参考 FASTA/BAM 用 `chr` 前缀。未归一化时 PRICE
ORF 会全 `N` 序列、零 psite/counts、且无法与其他工具 exact-match 合并（`id_key`
含 chrom）。`main()` 在所有文件解析完成后、合并/GTF 查表/取序列之前调用
`_make_chrom_normalizer(gtf_index)`（基于 `GTFIndex.chrom_names` + `_chrom_aliases()`；
幂等；GL/KI 脚手架名两侧一致故不动），并重建 `cand.id_key`（该键在 `__init__` 缓存
了 chrom）。修复后需重跑 UNIFY，**ORF ID 会整体重排**。

## stats.txt 解读

```
=== Input Statistics (raw, per tool and per sample) === By Tool:
Ribo-TISH: 767,610 / ORFquant: 181,630 / Ribotricer: 37,561  （human 62 样本示例）
```

某工具计数为 0 → 该工具结果没进 unify（检查输入路径/后缀解析）。
`stats_mode=auto`：≤512MB 且 ≤20000 ORF 用 preload，否则 stream
（rice 752MiB / 171,549 ORFs → stream 模式、4 workers、184 个 bedgraph 文件）。

## 实测资源

- rice 23 样本：35m26s、65.7GB 内存、metadata.tsv 271MB、369,938 ORFs
- human 62 样本手动 unify：8307s（~2.3h）、956,971 ORFs
- `--duckdb-db`/`--duckdb-memory-limit`（orfont 路径）：持久化 DB 与内存上限
