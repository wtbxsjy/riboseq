# unify_orf_predictions.py 升级说明 (v2.0)

本文档说明 `unify_orf_predictions.py` v2.0 的新功能和改进。

---

## 📊 主要改进

### 1. 新增统计字段

#### metadata.tsv 新增列
| 列名 | 含义 | 数据类型 | 示例值 |
|------|------|---------|--------|
| **tool_pvalues** | 各工具的原始 p-value | string | `ribotish:1.20e-10` |
| **total_reads** | 总读段数（所有样本） | integer | `150` |
| **unique_reads** | 唯一映射读段数 | integer | `142` |
| **total_psites** | 总 P-site 数 | integer | `45` |
| **unique_psites** | 唯一 P-site 数 | integer | `42` |
| **pN** | P-sites per nucleotide | float | `0.135135` |
| **unique_pN** | Unique P-sites per nucleotide | float | `0.126126` |

### 2. 保留原始 p-value

**以前**：Ribo-TISH 的 TisPvalue 被转换为 `-log10(p)`  
**现在**：
- `tool_scores`：仍然包含转换后的 `-log10(p)`（便于排序）
- `tool_pvalues`：**新增**，保留原始 p-value（便于统计分析）

**示例**：
```tsv
tool_scores: Ribo-TISH:9.92,ORFquant:23.000,Ribotricer:0.950
tool_pvalues: Ribo-TISH:1.20e-10
```

### 3. RiboseQC Bedgraph 整合

通过 `--bedgraph-dir` 和 `--sample-list` 参数，脚本可以自动从 RiboseQC bedgraph 文件计算统一的统计量。

---

## 🚀 使用方法

### 基本用法（不含 bedgraph）

```bash
python3 unify_orf_predictions.py \
  --gtf gencode.gtf \
  --fasta genome.fa \
  --ribotish "s1_pred.txt s2_pred.txt" \
  --ribotricer "s1_orfs.tsv s2_orfs.tsv" \
  --orfquant "s1_ORFs.gtf s2_ORFs.gtf" \
  --output unified \
  --min_len 10
```

**输出**：
- `unified.metadata.tsv`：包含 `tool_scores` 和 `tool_pvalues`，但统计列为 0
- `unified.bed`
- `unified.gtf`

### 高级用法（含 RiboseQC bedgraph）⭐

```bash
python3 unify_orf_predictions.py \
  --gtf gencode.gtf \
  --fasta genome.fa \
  --ribotish "s1_pred.txt s2_pred.txt" \
  --ribotricer "s1_orfs.tsv s2_orfs.tsv" \
  --orfquant "s1_ORFs.gtf s2_ORFs.gtf" \
  --output unified \
  --min_len 10 \
  --bedgraph-dir ./riboseqc_output \
  --sample-list "sample1,sample2"
```

**输出**：
- `unified.metadata.tsv`：包含完整的统计量（total_psites, unique_psites, pN, unique_pN 等）
- `unified.bed`
- `unified.gtf`

---

## 📁 bedgraph 文件要求

### 目录结构

```
riboseqc_output/
├── sample1_P_sites_plus.bedgraph
├── sample1_P_sites_minus.bedgraph
├── sample1_P_sites_uniq_plus.bedgraph
├── sample1_P_sites_uniq_minus.bedgraph
├── sample1_coverage_plus.bedgraph       # 可选
├── sample1_coverage_minus.bedgraph       # 可选
├── sample2_P_sites_plus.bedgraph
├── sample2_P_sites_minus.bedgraph
└── ...
```

### 文件命名规范

- **P-site bedgraph**: `{sample}_P_sites_{plus|minus}.bedgraph`
- **Unique P-site**: `{sample}_P_sites_uniq_{plus|minus}.bedgraph`
- **Coverage**: `{sample}_coverage_{plus|minus}.bedgraph`（可选）
- **Unique coverage**: `{sample}_coverage_uniq_{plus|minus}.bedgraph`（可选）

---

## 📊 输出示例

### metadata.tsv (v2.0)

```tsv
orf_id	chrom	strand	start	end	length_aa	exon_blocks	gene_id	transcript_id	tools	samples	tool_scores	tool_pvalues	total_reads	unique_reads	total_psites	unique_psites	pN	unique_pN	sequence
ORF_1_ENSG00000123	chr1	+	100000	100333	111	100000-100333	ENSG00000123	ENST00000456	Ribo-TISH,ORFquant	sample1,sample2	Ribo-TISH:9.921,ORFquant:45.000	Ribo-TISH:1.20e-10	150	142	45	42	0.135135	0.126126	ATGAAGCTG...
ORF_2_ENSG00000456	chr1	-	200000	200150	50	200000-200150	ENSG00000456	ENST00000789	Ribotricer	sample1	Ribotricer:0.950	NA	80	75	38	35	0.253333	0.233333	ATGCCCTGG...
```

---

## 🎯 数据分析示例

### 1. 筛选高质量 ORF（多工具 + 高 p-value 显著性）

```bash
# 使用 R/tidyverse
library(tidyverse)

orfs <- read_tsv("unified.metadata.tsv")

# 高置信度 ORF
high_conf <- orfs %>%
  filter(
    str_count(tools, ",") >= 1 |  # 多工具检测
    str_detect(tool_pvalues, "e-0[5-9]|e-1[0-9]")  # p < 1e-5
  ) %>%
  filter(unique_psites >= 10)  # 至少 10 个 unique P-sites

# 高翻译活性 ORF
high_trans <- orfs %>%
  filter(unique_pN > 0.1)  # pN > 0.1

write_tsv(high_conf, "high_confidence_orfs.tsv")
write_tsv(high_trans, "high_translation_orfs.tsv")
```

### 2. 比较不同工具的分数分布

```python
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("unified.metadata.tsv", sep="\t")

# 解析 tool_scores
def parse_scores(row):
    scores = {}
    if row['tool_scores'] != 'NA':
        for item in row['tool_scores'].split(','):
            tool, score = item.split(':')
            scores[tool] = float(score)
    return scores

df['scores_dict'] = df.apply(parse_scores, axis=1)

# 提取各工具分数
df['ribotish_score'] = df['scores_dict'].apply(lambda x: x.get('Ribo-TISH', None))
df['orfquant_score'] = df['scores_dict'].apply(lambda x: x.get('ORFquant', None))
df['ribotricer_score'] = df['scores_dict'].apply(lambda x: x.get('Ribotricer', None))

# 绘制分布图
fig, axes = plt.subplots(1, 3, figsize=(15, 4))

df['ribotish_score'].dropna().hist(bins=50, ax=axes[0])
axes[0].set_title('Ribo-TISH: -log10(p)')
axes[0].set_xlabel('Score')

df['orfquant_score'].dropna().hist(bins=50, ax=axes[1])
axes[1].set_title('ORFquant: P_sites')
axes[1].set_xlabel('Score')

df['ribotricer_score'].dropna().hist(bins=50, ax=axes[2])
axes[2].set_title('Ribotricer: phase_score')
axes[2].set_xlabel('Score')

plt.tight_layout()
plt.savefig('tool_score_distributions.pdf')
```

### 3. pN 值与 ORF 长度的关系

```R
library(ggplot2)

orfs <- read_tsv("unified.metadata.tsv")

ggplot(orfs, aes(x = length_aa, y = unique_pN)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "loess") +
  scale_x_log10() +
  scale_y_log10() +
  labs(
    title = "ORF Translation Efficiency vs Length",
    x = "ORF Length (aa)",
    y = "Unique pN (P-sites per nt)"
  ) +
  theme_minimal()

ggsave("pN_vs_length.pdf", width = 8, height = 6)
```

### 4. 差异翻译分析（DESeq2）

```R
library(DESeq2)
library(tidyverse)

orfs <- read_tsv("unified.metadata.tsv")

# 如果有多个样本的 bedgraph，可以重新运行 unify 脚本为每个样本单独计算
# 或者使用 scripts/R/quantify_orfs_from_psites.R

# 构建计数矩阵（示例：假设有每样本的 psites）
# count_matrix <- ...

# DESeq2 分析
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = sample_metadata,
  design = ~ condition
)

dds <- DESeq(dds)
res <- results(dds)

# 筛选差异翻译 ORF
sig_orfs <- res %>%
  as.data.frame() %>%
  rownames_to_column("orf_id") %>%
  filter(padj < 0.05, abs(log2FoldChange) > 1)

write_tsv(sig_orfs, "differential_translation_orfs.tsv")
```

---

## 🔧 故障排查

### Q1: bedgraph 文件未找到
**错误**：`Warning: Error reading bedgraph ...`

**解决**：
1. 检查 `--bedgraph-dir` 路径是否正确
2. 确认文件命名符合规范
3. 检查 `--sample-list` 中的样本名是否与文件名匹配

### Q2: tool_pvalues 列为 NA
**原因**：只有 Ribo-TISH 提供 p-value，Ribotricer 和 ORFquant 无 p-value

**正常**：
```
tool_pvalues: Ribo-TISH:1.20e-10  # 仅 Ribo-TISH 有值
```

### Q3: 统计列全为 0
**原因**：未提供 `--bedgraph-dir` 参数

**解决**：添加 bedgraph 参数或接受默认值 0（表示未计算）

---

## 📚 相关文档

- [ORF_SCORING_METRICS.md](ORF_SCORING_METRICS.md)：各工具评分机制详解
- [ORF_STATISTICAL_MODELS.md](ORF_STATISTICAL_MODELS.md)：统计检验模型汇总
- [README_WORKFLOW_PREP.md](../scripts/README_WORKFLOW_PREP.md)：工作流准备脚本文档

---

## 📝 版本历史

### v2.0 (2026-02-01)
- ✅ 新增 `tool_pvalues` 列（保留原始 p-value）
- ✅ 新增 bedgraph 统计功能（total_psites, unique_psites, pN, unique_pN）
- ✅ 支持 `--bedgraph-dir` 和 `--sample-list` 参数
- ✅ 改进 ORFCandidate 类，支持统计量合并
- ✅ 更新文档（本文件 + ORF_STATISTICAL_MODELS.md）

### v1.1 (2026-01-31)
- 新增 `tool_scores` 列
- GTF 输出添加 `sources`, `samples`, `num_tools` 属性

### v1.0 (初始版本)
- 基本 ORF 合并功能
- metadata.tsv, BED, GTF 输出

---

*最后更新：2026-02-01*
