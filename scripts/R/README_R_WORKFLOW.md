# R-based GENCODE ORF Analysis Workflow

完整的R语言解决方案，用于从P-site bedgraph和GENCODE注释进行ORF定量和多工具比较分析。

## 概述

本工作流提供两个核心R脚本：

1. **`quantify_orfs_from_psites.R`** - 从RiboseQC的P-site bedgraph重新定量统一的GENCODE ORFs
2. **`analyze_tool_comparison.R`** - 多工具方法学比较分析

## 依赖包安装

```r
# CRAN packages
install.packages(c("optparse", "tidyverse", "ggplot2", "data.table",
                   "pheatmap", "RColorBrewer", "ggvenn", "UpSetR"))

# Bioconductor packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("GenomicRanges", "rtracklayer", "ComplexHeatmap"))
```

## 工作流程

```
步骤1: 运行ORF预测工具 → 步骤2: 转换为GENCODE格式 →
步骤3: GENCODE映射统一 → 步骤4: P-site定量 → 步骤5: 工具比较分析
```

### 步骤1-3: 准备统一的GENCODE ORF注释

参考 `README.md` 中的bash脚本部分，使用：
- `14_ribotish_to_gencode.sh`
- `15_ribotricer_to_gencode.sh`
- `16_gencode_orf_mapper.sh`

输出文件：
- `Mouse_AllTools.orfs.gtf` - **用于定量**
- `Mouse_AllTools.orfs.out` - **用于工具比较**

### 步骤4: P-site定量 ⭐

#### 使用方法

```bash
Rscript quantify_orfs_from_psites.R \
  --gtf results_all_tools_integrated/Mouse_AllTools_Comparison.orfs.gtf \
  --bedgraph-dir riboseqc_results/ \
  --sample-pattern "(.+)_P_sites_unique_(plus|minus)\\.bedgraph$" \
  --outdir quantification_results \
  --min-count 5 \
  --threads 4
```

#### 参数说明

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `--gtf` | GENCODE ORF GTF文件 | 必需 |
| `--bedgraph-dir` | RiboseQC bedgraph文件目录 | 必需 |
| `--sample-pattern` | 提取样本名的正则表达式 | `(.+)_P_sites_unique_(plus\|minus)\\.bedgraph$` |
| `--outdir` | 输出目录 | `quantification_results` |
| `--min-count` | 最小P-site总数（过滤低表达ORF） | 5 |
| `--threads` | 线程数（用于并行处理） | 4 |

#### 输出文件

```
quantification_results/
├── orf_counts_raw.tsv          # 原始P-site计数（含ORF注释）
├── orf_counts_tpm.tsv          # TPM归一化计数
├── orf_counts_matrix.csv       # 纯计数矩阵（用于DESeq2）⭐
├── sample_summary_stats.tsv    # 样本级统计
├── orf_summary_stats.tsv       # ORF级统计
└── session_info.txt            # R会话信息
```

#### 文件格式示例

**`orf_counts_raw.tsv`**:
```
orf_id                    seqname  start    end      strand  gene_name  gene_type  Sample01  Sample02  ...
ENSMUST00000456328_uORF1  chr1     100000   100333   +       Gapdh      uORF       123       456       ...
```

**`orf_counts_matrix.csv`** (用于DESeq2):
```
          Sample01  Sample02  Sample03  ...
orf_id_1  123       456       789       ...
orf_id_2  234       567       890       ...
```

### 步骤5: 工具比较分析 ⭐

#### 使用方法

```bash
Rscript analyze_tool_comparison.R \
  --input results_all_tools_integrated/Mouse_AllTools_Comparison.orfs.out \
  --tools "ORFquant,RiboTISH,Ribotricer" \
  --outdir tool_comparison_results \
  --min-samples 10 \
  --min-tools 2
```

#### 参数说明

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `--input` | GENCODE ORF输出文件（.orfs.out） | 必需 |
| `--tools` | 工具名称（逗号分隔） | `ORFquant,RiboTISH,Ribotricer` |
| `--outdir` | 输出目录 | `tool_comparison_results` |
| `--min-samples` | 高置信度ORF最小样本数 | 10 |
| `--min-tools` | 高置信度ORF最小工具数 | 2 |

#### 输出文件

**数据文件**:
```
tool_comparison_results/
├── high_confidence_orfs.tsv        # 高置信度ORF（≥2工具+≥10样本）⭐
├── ORFquant_specific_orfs.tsv      # ORFquant特异性ORF
├── RiboTISH_specific_orfs.tsv      # Ribo-TISH特异性ORF
├── Ribotricer_specific_orfs.tsv    # Ribotricer特异性ORF
├── orf_biotype_by_tool.tsv         # ORF类型分布表
├── tool_comparison_summary.txt     # 总结报告
└── session_info.txt
```

**可视化图表**:
```
├── venn_diagram.pdf                # Venn图（2-3工具时）
├── upset_plot.pdf                  # UpSet图（显示所有交集组合）⭐
├── biotype_heatmap.pdf             # ORF类型分布热图
└── concordance_barplot.pdf         # 工具一致性柱状图
```

## 完整流程示例

### 场景：36个小鼠样本 × 3种工具

```bash
#!/bin/bash
# 完整的GENCODE注释和R分析流程

# ========== 步骤1-3: 准备GENCODE注释（bash脚本）==========
cd scripts/singularity_single_tool_tests

# 1.1 准备Ensembl目录
bash prepare_gencode_m35_ensembl.sh ~/Ens_GENCODE_M35

# 1.2 转换工具输出为GENCODE格式
for sample in Sample{01..36}; do
    # Ribo-TISH
    ./14_ribotish_to_gencode.sh \
      --sample ${sample} \
      --predict ribotish_results/${sample}_pred.txt \
      --fasta ~/genome/GRCm39.fa \
      --outdir ribotish_gencode

    # Ribotricer
    ./15_ribotricer_to_gencode.sh \
      --sample ${sample} \
      --tsv ribotricer_results/${sample}_translating_ORFs.tsv \
      --fasta ~/genome/GRCm39.fa \
      --outdir ribotricer_gencode
done

# 1.3 合并所有工具的ORF
cat orfquant_gencode/*.gencode.fa \
    ribotish_gencode/*.gencode.fa \
    ribotricer_gencode/*.gencode.fa \
    > merged_all_tools.fa

cat orfquant_gencode/*.gencode.bed \
    ribotish_gencode/*.gencode.bed \
    ribotricer_gencode/*.gencode.bed \
    > merged_all_tools.bed

# 1.4 运行GENCODE mapper（统一注释+去冗余）
./16_gencode_orf_mapper.sh \
  --project Mouse_AllTools_Comparison \
  --fasta merged_all_tools.fa \
  --bed merged_all_tools.bed \
  --ensembl-dir ~/Ens_GENCODE_M35 \
  --min-length 16 \
  --collapse-threshold 0.78 \
  --collapse-method longest_string \
  --outdir results_all_tools_integrated

# ========== 步骤4: P-site定量（R脚本）==========
cd ../../R

Rscript quantify_orfs_from_psites.R \
  --gtf ../singularity_single_tool_tests/results_all_tools_integrated/Mouse_AllTools_Comparison.orfs.gtf \
  --bedgraph-dir ~/riboseqc_results/ \
  --sample-pattern "(.+)_P_sites_unique_(plus|minus)\\.bedgraph$" \
  --outdir quantification_results \
  --min-count 5 \
  --threads 8

# ========== 步骤5: 工具比较分析（R脚本）==========
Rscript analyze_tool_comparison.R \
  --input ../singularity_single_tool_tests/results_all_tools_integrated/Mouse_AllTools_Comparison.orfs.out \
  --tools "ORFquant,RiboTISH,Ribotricer" \
  --outdir tool_comparison_results \
  --min-samples 10 \
  --min-tools 2

echo "✅ 完整流程完成！"
echo "定量结果: quantification_results/orf_counts_matrix.csv"
echo "工具比较: tool_comparison_results/high_confidence_orfs.tsv"
```

## 下游分析建议

### 1. 差异表达分析（DESeq2）

```r
library(DESeq2)
library(tidyverse)

# 读取计数矩阵
count_matrix <- read.csv("quantification_results/orf_counts_matrix.csv",
                         row.names = 1)

# 准备样本信息（根据你的实验设计修改）
coldata <- data.frame(
  sample = colnames(count_matrix),
  condition = rep(c("Control", "Treatment"), each = 18),  # 示例
  row.names = colnames(count_matrix)
)

# 创建DESeq数据集
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = coldata,
  design = ~ condition
)

# 过滤低表达ORF
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

# 运行差异分析
dds <- DESeq(dds)
res <- results(dds)

# 提取显著差异ORF
sig_orfs <- as.data.frame(res) %>%
  rownames_to_column("orf_id") %>%
  filter(padj < 0.05, abs(log2FoldChange) > 1)

write_tsv(sig_orfs, "differential_orfs.tsv")
```

### 2. 整合工具比较和表达数据

```r
library(tidyverse)

# 读取高置信度ORF
high_conf <- read_tsv("tool_comparison_results/high_confidence_orfs.tsv")

# 读取定量数据
orf_stats <- read_tsv("quantification_results/orf_summary_stats.tsv")

# 整合
integrated <- high_conf %>%
  left_join(orf_stats, by = "orf_id") %>%
  filter(mean_tpm > 1, n_samples_detected >= 10)

# 按ORF类型统计
integrated %>%
  group_by(orf_biotype) %>%
  summarise(
    n_orfs = n(),
    mean_expression = mean(mean_tpm),
    median_samples = median(n_samples_detected)
  )
```

### 3. 可视化

```r
library(ggplot2)
library(pheatmap)

# 读取TPM数据
tpm_data <- read_tsv("quantification_results/orf_counts_tpm.tsv")

# 选择高表达ORF
top_orfs <- orf_stats %>%
  arrange(desc(mean_tpm)) %>%
  head(50) %>%
  pull(orf_id)

# 准备热图矩阵
heatmap_matrix <- tpm_data %>%
  filter(orf_id %in% top_orfs) %>%
  select(orf_id, starts_with("Sample")) %>%
  column_to_rownames("orf_id") %>%
  as.matrix()

# 对数转换（添加伪计数）
heatmap_matrix_log <- log2(heatmap_matrix + 1)

# 绘制热图
pheatmap(heatmap_matrix_log,
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         main = "Top 50 Expressed ORFs",
         filename = "top_orfs_heatmap.pdf",
         width = 12, height = 10)
```

## 故障排查

### 问题1: 找不到bedgraph文件

**错误**: `No bedgraph files found matching pattern`

**解决**:
```bash
# 检查文件名格式
ls riboseqc_results/*bedgraph | head

# 调整--sample-pattern参数
# 示例：如果文件名是 "S1_psites_plus.bedgraph"
Rscript quantify_orfs_from_psites.R \
  --sample-pattern "(.+)_psites_(plus|minus)\\.bedgraph$" \
  ...
```

### 问题2: GTF中没有transcript特征

**警告**: GTF可能只包含CDS或exon特征

**说明**: 脚本会自动处理，按transcript_id合并CDS区域

### 问题3: 内存不足

**错误**: `cannot allocate vector of size...`

**解决**:
```bash
# 增加R内存限制
R --max-mem-size=32G

# 或分批处理样本
# 修改脚本，每次处理10个样本
```

### 问题4: 工具名称不匹配

**错误**: `No columns found for tool: XXX`

**解决**:
```bash
# 检查.orfs.out文件的列名
head -1 Mouse_AllTools_Comparison.orfs.out

# 调整--tools参数以匹配实际列名
# 示例：如果列名是"RiboTish_Sample01"而不是"RiboTISH_Sample01"
Rscript analyze_tool_comparison.R \
  --tools "ORFquant,RiboTish,Ribotricer" \
  ...
```

## 性能优化

### 并行处理

```r
# quantify_orfs_from_psites.R支持多线程
--threads 8  # 根据CPU核心数调整
```

### 数据压缩

```bash
# 压缩中间文件节省空间
gzip quantification_results/orf_counts_raw.tsv
```

## 输出文件总结

| 文件 | 用途 | 格式 |
|------|------|------|
| `orf_counts_matrix.csv` | DESeq2差异分析 | 纯数值矩阵 |
| `orf_counts_tpm.tsv` | 可视化、聚类 | 含注释的TPM值 |
| `high_confidence_orfs.tsv` | 高质量ORF子集 | 完整注释表 |
| `upset_plot.pdf` | 工具交集可视化 | PDF图表 |
| `orf_summary_stats.tsv` | ORF过滤参考 | 统计表 |

## 引用

如果使用本工作流，请引用：

- gencode-riboseqORFs: [GitHub](https://github.com/jorruior/gencode-riboseqORFs)
- GenomicRanges: Lawrence et al. (2013) doi:10.1371/journal.pcbi.1003118
- DESeq2: Love et al. (2014) doi:10.1186/s13059-014-0550-8

## 联系与支持

相关文档：
- bash脚本文档: `scripts/singularity_single_tool_tests/README.md`
- 转换器文档: `bin/README_CONVERTERS.md`
- 项目总结: `GENCODE_SESSION_UPDATE_20260117.md`

问题反馈：请参考项目主README
