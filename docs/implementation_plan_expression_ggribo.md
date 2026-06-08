# Implementation Plan: Expression Quantification + ggRibo + TE Bug Fix

**Target project:** `~/projects/Lishuqi` (Lishuqi_Ribo-seq_Human+MRSA)
**Working directory for development:** `/media/home/renzhe/riboseq/riboseq`

---

## Phase 0 — 现状审计 & TE Bug 定位

### 0.1 TE 分析问题诊断

**已确认的 Bug:**
- `deltate/plots/$filename` — Nextflow `$filename` 模板变量未被正确展开，说明 `main.nf` 的 `template` 处理有问题
- 缺少核心输出: `*.translation.deltate.results.tsv`, `*.dtegs.deltate.genes.tsv` 等均未生成

**症状分析:**
- R session 成功加载所有包 (sessionInfo.log 正常)
- `heatmap_zscores.tsv` 和 `heatmap_annotations.tsv` 有内容 → 说明 DESeq2 分析跑到了 heatmap 绘制阶段
- 但 heatmap PNG 未生成 → `plot_heatmap()` 可能在 `draw()` 时崩溃
- 检查 `ComplexHeatmap` 包是否正确安装
- 检查 `$filename` 模板变量 → 在 `deseq2_deltate/main.nf` 中需要确认 `'$counts'` 和 `'$filename'` 的引用

**TE 手工测试步骤:**
1. 进入 TE 工作目录，用已有的 `merged_counts.tsv` + `sample_sheet.csv` + contrast 手动执行 `deseq2_deltate.R`
2. 追踪错误位置
3. Fix R script 模板中的问题

### 0.2 可用数据清单

| 数据类型 | 路径 | 内容 |
|---------|------|------|
| Unified ORF metadata | `result/orf_unification/unified_orfs.metadata.tsv` | 592K rows, 27 columns |
| Unified ORF BED | `result/orf_unification/unified_orfs.bed.gz` | BED12 |
| ORF confidence | `result/orf_qc/qc_orf_confidence.tsv` | 491K rows, OCS + tier |
| Merged counts (TE) | `result/translational_efficiency/counts/merged_counts.tsv` | 491K × 12 samples |
| P-site bedgraph (8 Ribo samples) | `result/riboseqc/*_P_sites_{plus,minus}.bedgraph` | per-nt P-site positions |
| Coverage bedgraph (8 Ribo samples) | `result/riboseqc/*_coverage_{plus,minus}.bedgraph` | RPM-normalized coverage |
| Samplesheet | `result/translational_efficiency/counts/sample_sheet.csv` | 12 samples |
| Contrasts | `scripts/contrasts.csv` | treatment_vs_control |

---

## Phase 1 — 表达量定量脚本: `bin/quantify_orf_expression.py`

### 1.1 功能
从 RiboseQC P-site bedgraph 为每个 unified ORF 提取 per-sample 表达量。

### 1.2 接口设计

```bash
quantify_orf_expression.py \
    --orf-meta unified_orfs.metadata.tsv \     # 输入: ORF坐标表
    --orf-confidence orf_confidence.tsv \       # 输入: OCS scores (可选)
    --psites-dir riboseqc/ \                    # 输入: P-site bedgraph 目录
    --sample-pattern "Ribo_*_P_sites_plus.bedgraph" \  # 样本匹配模式
    --output expression_summary.tsv \           # 输出: per-ORF per-sample reads+pN
    --min-ocs 0.0 \                            # 过滤: 最低 OCS
    --workers 4                                 # 并行数
```

### 1.3 改造要点 (相对于 post_analysis)

| 方面 | post_analysis | pipeline 集成版 |
|------|--------------|----------------|
| ORF 输入 | `*result.tsv` (特殊格式) | `unified_orfs.metadata.tsv` (标准格式) |
| 样本发现 | 硬编码路径扫描 | `--sample-pattern` glob 或显式列表 |
| 坐标解析 | 手动 parse `P1_11:2512743_2512850:+:1:108` | 直接读 `chrom, start, end, strand` 列 |
| 染色体名 | strip P{version}_ 前缀 | 直接使用 metadata 中的 chrom |
| 输出格式 | AMP-specific 列 | 通用格式: `orf_id | chrom | start | end | strand | sample1_reads | sample1_pN | ...` |
| 依赖 | gencode 结果做样本过滤 | 直接用 OCS + bedgraph 数据存在性 |

### 1.4 实现

在 `bin/quantify_orf_expression.py` 中:

```python
#!/usr/bin/env python3
"""
从 RiboseQC P-site bedgraph 提取每个 unified ORF 的 per-sample 表达量。
不依赖 gencode 结果或 AMP 特定格式 —— 接受标准的 unified ORF metadata 作为输入。
"""
# 核心类: ORFExpressionQuantifier
# - load_orf_metadata(): 读取 unified metadata
# - discover_samples(): 扫描 bedgraph 目录发现样本
# - query_orf_psites(orf_row, sample): awk 查询 bedgraph → reads + pN
# - run(): 遍历 ORF, 并行查询, 输出矩阵
```

---

## Phase 2 — RPKM/TPM 计算: `bin/calc_orf_rpkm_tpm.py`

### 2.1 功能
从 RiboseQC **coverage** bedgraph (已 RPM 归一化) 计算每个 ORF 的 RPKM 和 TPM。

### 2.2 接口

```bash
calc_orf_rpkm_tpm.py \
    --expression expression_summary.tsv \       # Phase 1 输出 (含 ORF 坐标)
    --coverage-dir riboseqc/ \                  # coverage bedgraph 目录
    --sample-pattern "*_coverage_plus.bedgraph" \
    --output expression_rpkm_tpm.tsv \
    --workers 4
```

### 2.3 计算逻辑 (与 post_analysis 一致)

```
RPKM = sum(coverage_RPM_over_orf) / (orf_length / 1000)
TPM  = RPKM / sum(all_RPKM) * 1e6
```

coverage bedgraph 的值已经是 RPM (reads per million mapped reads)，无需再获取 library size。

---

## Phase 3 — 合并注释: `bin/combine_orf_annotations.py`

### 3.1 功能
将 OCS、分类结果 (如有)、TE counts、表达量 合并为一张完整 ORF 注释表。

### 3.2 接口

```bash
combine_orf_annotations.py \
    --orf-meta unified_orfs.metadata.tsv \
    --orf-confidence orf_confidence.tsv \       # 可选
    --expression expression_summary.tsv \        # Phase 1
    --rpkm-tpm expression_rpkm_tpm.tsv \         # Phase 2
    --te-counts merged_counts.tsv \              # 可选: 来自 TE
    --classification-dir orf_classification/ \   # 可选: gencode/orfquant/orftype
    --output orf_expression_combined.tsv
```

输出列:
```
orf_id | chrom | start | end | strand | gene_id | orf_type |
ocs | tier | s_translation | s_agreement | s_coverage | s_periodicity | s_readlevel |
detecting_tools | n_detecting |
sample1_reads | sample1_pN | sample1_rpkm | sample1_tpm | ...
gencode_biotype | orfquant_category | orftype_class |
te_count_Ribo_11 | te_count_RNA_11 | ... (可选)
```

---

## Phase 4 — ggRibo 独立脚本: `bin/plot_orf_ggribo.R`

### 4.1 功能
接受 pipeline 输出，为指定 ORF 批量生成 ggRibo 覆盖度图。

### 4.2 接口

```bash
plot_orf_ggribo.R \
    --orf-ids ORF_4011_ENSG...,ORF_17515_ENSG... \    # ORF IDs (逗号分隔)
    --orf-ids-file top_dtegs.txt \                      # 或文件列表
    --orf-meta unified_orfs.metadata.tsv \              # ORF 坐标
    --psites-dir riboseqc/ \                            # P-site bedgraph 目录
    --gtf reference.gtf \                               # 参考 GTF (用于 exon 结构)
    --output-dir ggribo_plots/ \
    --samples Ribo_11,Ribo_12,Ribo_14 \                # 指定样本 (默认 top 3 by reads)
    --sample-group group \                              # 分组标注 (可选)
    --extend 200 \                                      # 扩展 bp
    --n-top-orfs 20 \                                   # 若无 --orf-ids, 按 OCS 取 top N
    --container docker://.../ggribo:latest               # 容器 (可选)
```

### 4.3 关键改造 (相对于 post_analysis)

| 方面 | post_analysis | pipeline 集成版 |
|------|--------------|----------------|
| ORF 选取 | 硬编码 orf_defs | `--orf-ids` 动态指定 |
| 样本选取 | 硬编码 `find_top_samples()` | `--samples` 显式指定 |
| 染色体名 | 硬编码 `"2"` | 从 metadata 自动读取 |
| 绘图注释 | 硬编码 pN/RPKM 值 | 从 expression_summary 自动读取 |
| 输出命名 | `{org}_{cat}_coverage.png` | `{orf_id}_{samples}_ggribo.png` |

### 4.4 容器方案
使用 existing 的 RiboseQC 容器或独立的 ggRibo 容器:
- ggRibo 依赖 `txdbmaker`, `GenomicFeatures`, `ggRibo`
- 已通过自定义 `Range_info` R6 class 绕过命名空间冲突
- 建议在 `containers/Singularity.ggribo.def` 中构建专用容器

---

## Phase 5 — TE Bug 修复

### 5.1 定位

在 `modules/local/deseq2_deltate/main.nf` 的 `script:` 块中使用了 `template 'deseq2_deltate.R'`。
模板变量 `$filename` 出现在 R 代码中的 Nextflow 变量引用处。

### 5.2 手工测试步骤

```bash
# 1. 提取 R 脚本模板，替换变量后手动执行
cd /tmp/te_test
cp ~/projects/Lishuqi/result/translational_efficiency/counts/merged_counts.tsv .
cp ~/projects/Lishuqi/result/translational_efficiency/counts/sample_sheet.csv .

# 2. 手动替换模板变量并运行
sed -e 's/\$task\.ext\.prefix/treatment_vs_control/g' \
    -e 's/\$counts/merged_counts.tsv/g' \
    -e 's/\$samplesheet/sample_sheet.csv/g' \
    -e 's/\$contrast_variable/group/g' \
    -e 's/\$reference/control/g' \
    -e 's/\$target_level/treatment/g' \
    -e 's/\$task\.cpus/4/g' \
    -e "s/\$filename/__BUG__/g" \   # <-- 定位 bug
    ~/riboseq/riboseq/modules/local/deseq2_deltate/templates/deseq2_deltate.R \
    > test_deseq2.R

# 3. 检查 test_deseq2.R 中的 $filename 残留
grep -n '$filename\|__BUG__' test_deseq2.R

# 4. 运行修复后的 R 脚本
Rscript test_deseq2.R
```

### 5.3 可能的 Bug 原因

1. **`$filename` 未定义**: R 模板中某处引用了 `$filename` 但 Nextflow 脚本中没有此变量
2. **复杂热图崩溃**: `ComplexHeatmap::draw()` 在无 display 的环境中可能失败 (需 `ragg` 或 `Cairo`)
3. **`ggsave` 文件名问题**: `paste0(prefix, ".fold_change.png")` 中 `$filename` 可能被误传入

### 5.4 修复方案

在 `deseq2_deltate.R` 模板中:
- 删除所有未定义的 `$filename` 引用
- 为 `ComplexHeatmap::draw()` 添加错误处理
- 确保所有 `ggsave()` 使用正确的 `paste0(prefix, ...)` 模式

---

## Phase 6 — 手工端到端测试

### 6.1 测试数据准备

```bash
export TEST_DIR=/tmp/expression_quant_test
export RESULT_DIR=~/projects/Lishuqi/result

mkdir -p $TEST_DIR/{inputs,outputs}

# 输入准备
cp $RESULT_DIR/orf_unification/unified_orfs.metadata.tsv $TEST_DIR/inputs/
cp $RESULT_DIR/orf_qc/qc_orf_confidence.tsv $TEST_DIR/inputs/
cp $RESULT_DIR/translational_efficiency/counts/merged_counts.tsv $TEST_DIR/inputs/
# 仅 Ribo samples 的 P-site bedgraph
for f in $RESULT_DIR/riboseqc/Ribo_*_P_sites_*.bedgraph; do
    ln -s "$f" $TEST_DIR/inputs/
done
```

### 6.2 运行 Phase 1: 表达量提取

```bash
cd $TEST_DIR
python3 ~/riboseq/riboseq/bin/quantify_orf_expression.py \
    --orf-meta inputs/unified_orfs.metadata.tsv \
    --orf-confidence inputs/qc_orf_confidence.tsv \
    --psites-dir inputs/ \
    --sample-pattern "Ribo_*_P_sites_plus.bedgraph" \
    --output outputs/expression_summary.tsv \
    --min-ocs 0.4 \
    --max-orfs 1000 \
    --workers 4
```

### 6.3 运行 Phase 2: RPKM/TPM

```bash
python3 ~/riboseq/riboseq/bin/calc_orf_rpkm_tpm.py \
    --expression outputs/expression_summary.tsv \
    --coverage-dir inputs/ \
    --output outputs/expression_rpkm_tpm.tsv \
    --workers 4
```

### 6.4 运行 Phase 3: 合并注释

```bash
python3 ~/riboseq/riboseq/bin/combine_orf_annotations.py \
    --orf-meta inputs/unified_orfs.metadata.tsv \
    --orf-confidence inputs/qc_orf_confidence.tsv \
    --expression outputs/expression_summary.tsv \
    --rpkm-tpm outputs/expression_rpkm_tpm.tsv \
    --te-counts inputs/merged_counts.tsv \
    --output outputs/orf_expression_combined.tsv
```

### 6.5 运行 Phase 4: ggRibo 绘图

```bash
# 从 High confidence ORFs 中取 top 5
head -6 outputs/expression_summary.tsv | tail -5 | cut -f1 > outputs/top5_orfs.txt

Rscript ~/riboseq/riboseq/bin/plot_orf_ggribo.R \
    --orf-ids-file outputs/top5_orfs.txt \
    --orf-meta inputs/unified_orfs.metadata.tsv \
    --expression outputs/expression_summary.tsv \
    --psites-dir inputs/ \
    --gtf ~/riboseq/reference_db/build_output/ensembl/Ens113_homo_sapiens/SORTED_TRANSCRIPTOME_GTF \
    --output-dir outputs/ggribo_plots/ \
    --extend 200
```

### 6.6 运行 Phase 5: TE Fix

```bash
cd $TEST_DIR/te_fix
# 复制 R 模板并手工替换变量
# 运行修复后的 deseq2 R 脚本
# 验证输出文件完整性
```

---

## 执行顺序总结

```
Phase 0: 现状审计
    ├── 0.1: 确认所有输入数据可访问 ✓ (已完成)
    ├── 0.2: 确认 TE bug 类型 → Phase 5 并行
    │
Phase 1: quantify_orf_expression.py (2-3h)
    │
Phase 2: calc_orf_rpkm_tpm.py (1h)
    │
Phase 3: combine_orf_annotations.py (1h)
    │
Phase 4: plot_orf_ggribo.R (2h, 独立)
    │
Phase 5: TE bug fix + 手工验证 (1h)
    │
Phase 6: 端到端手工测试 (1h)
```

**总计约 1 天工作量。**
