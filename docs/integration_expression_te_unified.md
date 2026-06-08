# 表达量定量 × TE 分析 × ggRibo: 统一整合评估

## 1. 现有 TE 模块分析

### 1.1 架构

```
TE_ANALYSIS subworkflow:
  ├── QUANTIFY_ORFS      → featureCounts 对每个样本 BAM × unified ORF BED 做计数
  ├── MERGE_COUNTS        → 合并 per-sample counts 为 gene×sample 矩阵
  └── DESEQ2_DELTATE      → DESeq2 interaction model: ~ contrast + seq_type + contrast:seq_type
```

### 1.2 关键设计特征

| 方面 | 实现 |
|------|------|
| 定量方法 | **featureCounts** (read-level count, 从 BAM 统计) |
| 输入 | RNA-seq + Ribo-seq BAM 全部样本 |
| 注释 | unified ORF BED12 (来自 UNIFY_ORF_PREDICTIONS) |
| 分析方法 | DESeq2 interaction model (deltaTE) |
| 输出 | DTEGs 分类 (translation / mRNA_abundance / intensified / buffering) |
| 依赖 | subread 容器 (QUANTIFY_ORFS), DESeq2 容器 (DESEQ2_DELTATE) |

### 1.3 样本选择逻辑

```
ch_te_bams = RNA-seq BAMs (unfiltered) + Ribo-seq BAMs (post-sORF-filtered)
              ↓
         meta.sample_type = "riboseq" | "rnaseq"  (用于 DESeq2 中的 interaction term)
              ↓
         meta.group = treatment/control             (来自 samplesheet)
```

---

## 2. 两个定量体系的对比

### 2.1 核心区别

| 维度 | TE 分析 (featureCounts) | post_analysis (bedgraph 查询) |
|------|------------------------|-------------------------------|
| **定量对象** | Total reads 覆盖 ORF 区域 (raw count) | P-site reads 精确定位于 P-site 位置 |
| **数据源** | BAM 文件 | RiboseQC 预处理的 bedgraph (已计算 P-site 偏移) |
| **归一化** | DESeq2 内部 (median-of-ratios) | RiboseQC 已做 RPM; 再算 RPKM/TPM |
| **信号类型** | 全 read 覆盖 (含 rRNA/tRNA 污染可能) | 仅 P-site 信号 (更精确的翻译信号) |
| **用途** | 统计推断: 差异翻译 (deltaTE) | 表达量描述: 绝对值、pN、排序 |
| **输出行列** | gene × sample (raw count matrix) | ORF × sample (reads/pN/RPKM/TPM matrix) |
| **统计框架** | DESeq2 GLM + interaction term | 无; 纯定量 |

### 2.2 互补性分析

```
                        unified ORF set
                       /                \
                      /                  \
         featureCounts (BAM)      RiboseQC bedgraph
              │                         │
              ▼                         ▼
       raw count matrix           P-site reads + pN
              │                         │
              ▼                         ▼
       DESeq2 deltaTE              RPKM / TPM / pN
       (差异翻译推断)              (绝对表达量描述)
              │                         │
              └─────────┬───────────────┘
                        ▼
              合并注释表: OCS + tier + 分类 + counts + RPKM + pN
```

**核心互补性**:
- TE 告诉你"哪些 ORF 在条件间翻译效率有差异"
- Expression quantification 告诉你"每个 ORF 在每个样本中的绝对表达水平"
- 两者从不同角度回答不同问题，数据完全互补

---

## 3. 整合方案: 统一 ORF 表达量模块

### 3.1 方案: 扩展现有 TE 流程，新增 P-site 定量

保持现有 TE 分析不变，新增一个平行的 **QUANTIFY_ORFS_PSITE** 模块:

```
QUANTIFY_ORFS (featureCounts on BAM)     QUANTIFY_ORFS_PSITE (bedgraph query)
       │                                           │
       ▼                                           ▼
  raw count matrix                          P-site reads + pN matrix
       │                                           │
       ▼                                           ▼
  MERGE_COUNTS                              MERGE_PSITE_COUNTS
       │                                           │
       ▼                                           ▼
  DESeq2 deltaTE                             RPKM/TPM calculation
       │                                           │
       ▼                                           ▼
  DTEG classification                       expression_summary.tsv
                                                    │
       └─────────────────┬──────────────────────────┘
                         ▼
                COMBINE_ORF_ANNOTATIONS
                         │
                         ▼
         orf_expression_combined.tsv
         (OCS + tier + 分类 + counts + pN + RPKM + TPM + DTEG status)
```

### 3.2 关键数据流

**共享的输入**:
1. Unified ORF BED — 两个模块都需要
2. 样本列表 — 两个模块共享 (来自 samplesheet)
3. RiboseQC bedgraph — 仅新的 P-site 定量需要 (已在 pipeline 中)

**TE 已经产出的**:
- `merged_counts.tsv`: gene × sample raw count matrix
- 这个可以直接作为输入，避免重复 featureCounts 调用

**新的 P-site 定量产出**:
- `expression_summary.tsv`: ORF × sample (reads, mean_cov, max_cov, pN)
- `expression_rpkm_tpm.tsv`: ORF × sample (RPKM, TPM)

---

## 4. 具体实现路径

### 4.1 新增模块

```
modules/local/quantify_orf_psites/
├── main.nf                    # 查询 RiboseQC bedgraph, 输出 per-ORF per-sample 表达量
├── meta.yml
└── environment.yml            # Python + pandas

modules/local/combine_orf_annotations/
├── main.nf                    # 合并 counts + expression + OCS + classification → 完整表
└── meta.yml
```

### 4.2 改造后的数据流 (workflow 内)

```nextflow
// Step 1: featureCounts — 已有
QUANTIFY_ORFS(ch_te_bams, ch_unified_bed)

// Step 2: P-site quantification — 新增
QUANTIFY_ORF_PSITES(
    ch_unified_meta,           // unified ORF metadata
    ch_psites_bedgraph,        // 所有样本 P-site bedgraph (来自 RiboseQC)
    ch_coverage_bedgraph       // 所有样本 coverage bedgraph
)

// Step 3: 合并两套定量结果
COMBINE_ORF_ANNOTATIONS(
    MERGE_COUNTS.out.counts,          // raw count matrix
    QUANTIFY_ORF_PSITES.out.psites,   // P-site reads/pN matrix
    QUANTIFY_ORF_PSITES.out.rpkm_tpm, // RPKM/TPM matrix
    ORF_QC.out.confidence,            // OCS + tier
    CLASSIFY_ORFS.out.results          // classification results
)

// Step 4: TE 分析保持不变
DESEQ2_DELTATE(...)
```

### 4.3 QUANTIFY_ORF_PSITES 设计

```nextflow
process QUANTIFY_ORF_PSITES {
    tag "${meta.id}"  // or "quantify_psites"
    label 'process_medium'

    input:
    path unified_meta              // ORF 元数据
    path psites_bedgraph_plus      // 所有样本的 P_sites_plus.bedgraph (收集)
    path psites_bedgraph_minus     // 所有样本的 P_sites_minus.bedgraph
    path coverage_bedgraph_plus    // 所有样本的 coverage_plus.bedgraph
    path coverage_bedgraph_minus   // 所有样本的 coverage_minus.bedgraph
    val sample_list                // 样本名称列表 (来自 meta)

    output:
    path "*_expression_summary.tsv"  , emit: psites
    path "*_expression_rpkm_tpm.tsv" , emit: rpkm_tpm

    script:
    """
    # Phase 1: Extract per-ORF P-site reads + pN from bedgraph
    quantify_orf_expression.py \
        --orf-meta ${unified_meta} \
        --psites-plus ${psites_bedgraph_plus} \
        --psites-minus ${psites_bedgraph_minus} \
        --samples ${sample_list.join(',')} \
        --output ${prefix}_expression_summary.tsv

    # Phase 2: Calculate RPKM/TPM from coverage bedgraph
    calc_orf_rpkm_tpm.py \
        --expression ${prefix}_expression_summary.tsv \
        --coverage-plus ${coverage_bedgraph_plus} \
        --coverage-minus ${coverage_bedgraph_minus} \
        --samples ${sample_list.join(',')} \
        --output ${prefix}_expression_rpkm_tpm.tsv
    """
}
```

### 4.4 与 TE 分析的重用点

| 组件 | TE 使用 | 新的 P-site 定量 | 是否可共享 |
|------|---------|-----------------|-----------|
| Unified ORF BED | featureCounts 注释 | ORF 坐标查询 | ✅ 同一文件 |
| 样本列表 | BAM channel 中的 meta | bedgraph channel 中的 meta | ✅ 同一来源 |
| 容器 | subread | python:pandas | ❌ 不同 |
| 输出 | raw count matrix | expression matrix | ✅ 互补合并 |

---

## 5. ggRibo 可视化与 TE 的关系

ggRibo 与 TE 是完全独立的功能层:

| 层级 | 功能 | 输入 | 输出 |
|------|------|------|------|
| **定量层** | TE (featureCounts) + Expression (bedgraph) | BAM + bedgraph | 数值矩阵 |
| **统计层** | DESeq2 deltaTE | count matrix | DTEG 分类 |
| **可视化层** | ggRibo | P-site bedgraph + GTF | ORF 覆盖图 |

ggRibo 图对于验证 TE 结果特别有价值:

> 例如: TE 报告某个 ORF 在 treatment 中 translationally upregulated，ggRibo 图可以直接展示该 ORF 在 treatment vs control 样本中的 P-site 覆盖差异，提供视觉证据。

### ggRibo 集成建议

保持为独立的 post-hoc 工具 (option 2 from assessment)，但:
1. 脚本放置于 `bin/plot_orf_ggribo.R`
2. 接受 `--orf-ids` (ORF ID 列表, 可从 TE DTEgs 或 high-OCS ORFs 中选择)
3. 接受 `--samples` (指定哪些样本绘图)
4. 接受 `--group-by` (按分组上色, 如 treatment vs control)

这样用户可以在 pipeline 完成后，选取感兴趣的 DTEGs 或 high-confidence ORFs 批量出图:

```bash
# 从 TE 结果中取 top DTEGs 的 ORF ID
cut -f1 *_translation.deltate.genes.tsv | head -20 > top_dtegs.txt

# 用 ggRibo 绘图
plot_orf_ggribo.R \
    --orf-ids top_dtegs.txt \
    --psites-dir results/riboseqc/ \
    --gtf reference.gtf \
    --group-by treatment \
    --output ggribo_plots/
```

---

## 6. 实施优先级汇总

```
Phase 1 ── QUANTIFY_ORF_PSITES + expression matrix
           │  (新建 module, 2-3 天)
           │  输入: unified ORF meta + RiboseQC bedgraphs
           │  输出: expression_summary.tsv + rpkm_tpm.tsv
           ▼
Phase 2 ── COMBINE_ORF_ANNOTATIONS
           │  (新建 module, 0.5 天)
           │  输入: counts + psites + rpkm_tpm + OCS + classification
           │  输出: orf_expression_combined.tsv
           ▼
Phase 3 ── 改进 TE 的输入
           │  (轻量改动: TE 可选使用 per-ORF expression 矩阵替代重跑 featureCounts)
           │  如果 featureCounts 的 raw count 与 P-site quantitative 高度相关,
           │  TE 可直接使用 expression_summary 中的 reads 列
           ▼
Phase 4 ── ggRibo 独立脚本
              (bin/plot_orf_ggribo.R, 1 天)
              无 pipeline 集成压力, post-hoc 使用
```

---

## 7. 结论

**TE 分析和 post_analysis 表达量定量天然互补，整合价值高**:

1. **TE 分析** 回答: "条件 A vs B 中哪些 ORF 翻译效率有差异?" (统计推断)
2. **表达量定量** 回答: "每个 ORF 在每个样本中表达多少?" (绝对定量)
3. **ggRibo 可视化** 提供: "感兴趣的 ORF 在基因组上的 P-site 覆盖模式" (视觉验证)

三者在数据层面共享: unified ORF 注释 + RiboseQC bedgraph + BAM。建议以 `QUANTIFY_ORF_PSITES` 为核心新模块，与现有 `QUANTIFY_ORFS` (featureCounts) 并行运行，在 `COMBINE_ORF_ANNOTATIONS` 中合并为一套完整的 ORF 注释输出。
