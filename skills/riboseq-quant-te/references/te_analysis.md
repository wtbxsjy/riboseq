# TE/ΔTE 分析（deseq2_deltate）

⚠️ **现状**：代码完整但**从未在任何真实项目启用过**——所有 run 目录的
samplesheet 只有 5 列（无 treatment/group），run_pipeline.sh 均无 `--contrasts`，
`.nextflow.log` 无 DESEQ2_DELTATE 进程记录。首次启用时按本文件配置，跑通后回填实测数字。

## 启用条件

1. samplesheet 增加分组列（如 `condition`：treated/control）
2. 制备 contrasts CSV 并传 `--contrasts`：

```csv
id,variable,reference,target
lps_vs_nt,condition,CHX_NT,CHX_LPS
```

- 列要求：`id`（代码必需）、`variable`（= 样本表列名）、`reference`、`target`
  （variable 列的两个取值），可选 `batch`（也必须是样本表列名）
- workflow 解析：`.splitCsv(header:true).map{ [meta, row.variable, row.reference, row.target] }`

## 流程（subworkflows/local/te_analysis.nf）

```
QUANTIFY_ORFS（featureCounts，subread 2.1.1）→ MERGE_COUNTS → DESEQ2_DELTATE
```

- **QUANTIFY_ORFS**（modules/local/quantify_orfs/）：BED12→SAF（awk：GeneID=$4、
  Start=$2+1、End=$3、Strand=$6），
  `featureCounts -a annotation.saf -F SAF -t exon -g GeneID -s 0 -T {cpus} --minOverlap 1`
  → `{outdir}/translational_efficiency/counts/{sample}_counts.tsv`
- **MERGE_COUNTS**：内嵌 R 把 `*_counts.tsv`（Geneid + 末列 count）按 Geneid 合并，
  NA→0 → `merged_counts.tsv`；samplesheet 复制为 `sample_sheet.csv`
- **DESEQ2_DELTATE**：R 模板 `templates/deseq2_deltate.R`（`template` 渲染，
  `task.ext.args`/`task.ext.prefix` 注入）

## 核心模型

- 列名约定：`sample_id_col="sample"`、`seq_type_col="type"`、`gene_id_col="gene_id"`
- seq_type 自动识别：ribo = grep `ribo|rp|fp`；rna = grep `rna|mrna|total|lncrna`
  （`type` 列的值决定哪些样本算 RPF、哪些算 mRNA）
- 设计公式：`~ [batch +] condition + type + condition:type`
  （Chothani et al. 2019 deltaTE 方法）；**交互项系数 = TE**
- DESeq 调用：`DESeq(dds_combined, fitType="parametric", sfType="poscounts")`
- 交互项系数：`grep(paste0(contrast_var, "_", target), resultsNames(dds))` 取第一个匹配
- lfcShrink：默认 `type="apeglm"`，正确写法是**同时传 coef 与 res**：
  `lfcShrink(dds, coef=coef, res=res, type=shrink_type)`（规避 apeglm
  "coef does not uniquely specify" 报错）
- 分类（anota2seq 风格，lfc_threshold_te/rna/ribo 默认全 0）：
  `intensified`（TE+ribo+RNA 全显著同向）/ `buffering`（异向）/
  `translation`（TE+ribo 显著、RNA 不显著）/ `mRNA_abundance`（仅 ribo+RNA）/
  `dteg_other` / `other`

## 参数（nextflow.config L196-215）

| 参数 | 默认 | 说明 |
|---|---|---|
| `skip_te_analysis` | false | 关闭整个 TE 子流程 |
| `te_method` | deltate | deltate \| anota2seq |
| `te_lfc_threshold` | 0.2630344 | anota2seq 风格效应量阈值 |
| `rna_lfc_threshold` / `ribo_lfc_threshold` | 0 | 分类用 |
| `te_prefilter_min_nonzero` | 2 | 每 seq_type 最小非零样本数（绝对） |
| `te_prefilter_min_frac` | 0.2 | 每 seq_type 最小非零样本比例 |
| `extra_deltate_args` | null | `"--shrinkage_type normal --shrink_lfc true"` 等，进 R 模板 parse_args |
| `deseq2_container` | null | 自定义容器（默认 wave 拉 bioconductor-deseq2 镜像） |

**预过滤**（防崩溃，commit d9bb5d3）：`min_ribo_nonzero = max(prefilter_min_nonzero,
ceiling(n_ribo × prefilter_min_frac))`，ribo 与 rna 两侧都达标才保留；
`if (sum(keep) < 10) stop("Too few genes after filtering")`。这是
`estimateSizeFactors: every gene contains at least one zero` 崩溃的修复
（Lishuqi 数据实测 3,122/491,581 基因过预过滤）。

## 输出（`{result}/translational_efficiency/deltate/`，图在 `plots/`）

- `{p}.translation.deltate.results.tsv`（TE）、`{p}.translated_mRNA.deltate.results.tsv`
  （Ribo）、`{p}.total_mRNA.deltate.results.tsv`（RNA）
  列：gene_id, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj
- 基因列表：`{p}.{dtegs|mRNA_abundance|translation|intensified|buffering}.deltate.genes.tsv`
- 图：fold_change.png、interaction_p_distribution.png、pca_ribo/rna.{png,tsv}、
  heatmap.{png,tsv}；另 DESeqDataSet.rds、R_sessionInfo.log、versions.yml

## 坑

1. **apeglm→normal**（会话经验，仓库无记录）：模板默认 apeglm 至今未改；若 apeglm
   报错/不收敛，用 `--extra_deltate_args "shrinkage_type=normal --shrink_lfc=true"`
   覆盖（仓库模板里没有从 apeglm 切 normal 的历史注释，此条来自对话经验）
2. **estimateSizeFactors 崩溃**：见上，靠预过滤 + `sfType="poscounts"` 解决
3. **safe_vst**：样本数少时 vst() 抛 "less than 'nsub'" → 模板自动回退
   `varianceStabilizingTransformation(fitType="mean")`
4. **merge_counts fallback**：`.collect()` 的 List<Path> 在 R 端 Sys.glob 失败时
   用 `list.files(pattern="_counts\\.tsv$", recursive=TRUE)` 兜底
5. **容器**：`containers/Singularity.r_te_analysis.def`（R 4.5.3 + DESeq2, apeglm,
   ComplexHeatmap）；构建 `apptainer build deseq2_deltate.sif
   containers/Singularity.r_te_analysis.def`

## 手工单对比模板

`scripts/quant_analysis.R`：CHX_LPS vs CHX_NT 的 `~ condition` 单因子 DESeq2
（无交互项、无 lfcShrink），从 RiboseQC unique P-site bedgraph 用 foverlaps 计数。
适合快速试点（改 L26-38 路径与样本注释），正式分析仍走 ΔTE 模块。
