# Expression/Psite Quantification & ggRibo Integration Assessment

## 1. 现有 post_analysis 脚本分析

### 1.1 核心脚本及功能

```
post_analysis/scripts/
├── extract_amp_expression_smart.py   # ★ 核心: 从 RiboseQC bedgraph 提取每个 ORF 的 reads/pN
├── extract_amp_expression.py         #    早期版本 (全量查询, 无 gencode 预过滤)
├── calc_orf_rpkm.py                  # ★ 核心: 从 RiboseQC coverage bedgraph 计算 RPKM/TPM
├── consolidate_bedgraph.py           #    辅助: 合并零散 per-ORF bedgraph 为索引文件
├── plot_orf_ggribo.R                 # ★ 核心: ggRibo 覆盖度图 (代表性 ORF)
├── plot_sorf_strict_batch.R          #    批量绘制 sORF ggRibo 图
├── export_full_xlsx.R                #    导出 Excel 综合报告
└── match_classification.R            #    辅助: 将分类结果匹配到表达数据
```

### 1.2 数据依赖链

```
RiboseQC 输出 (*_P_sites_*.bedgraph, *_coverage_*.bedgraph)
    │
    ├──→ extract_amp_expression.py  ──→ expression_summary.tsv (per-ORF reads, pN)
    │                                       │
    │                                       ├──→ amp_expression_report.qmd (分布分析)
    │                                       │
    │                                       └──→ calc_orf_rpkm.py ──→ rpkm_tpm.tsv
    │
    └──→ plot_orf_ggribo.R (直接从 bedgraph + GTF)
            └──→ ggRibo 覆盖度 PNG
```

### 1.3 关键技术要点

| 方面 | 现有实现 | 备注 |
|------|---------|------|
| ORF 输入 | AMP 预测结果 (`*result.tsv`) | 需改为 unified ORF metadata |
| 样本发现 | gencode 结果列名 或 扫描 bedgraph 文件 | 可从 samplesheet 或 RiboseQC meta 获取 |
| Bedgraph 查询 | `awk` 按染色体+坐标区间过滤 | 高效, 依赖 tabix 排序 |
| pN 计算 | `max(coverage) / mean(coverage)` | 简单有效 |
| RPKM | `sum(RPM per nt) / length_kb` | coverage bedgraph 已经是 RPM 归一化 |
| ggRibo 绘图 | 自定义 `Range_info` R6 class + 临时 GTF | 绕过了 ggRibo 命名空间冲突 |

---

## 2. 现有 Pipeline 中已有的基础设施

### 2.1 已有的输入数据

| 数据 | 来源 | Channel |
|------|------|---------|
| P-site bedgraph (± strand) | `RIBOSEQC_ANALYSIS.out.psites_bedgraph` | `[meta, path(*_P_sites_*.bedgraph)]` |
| Coverage bedgraph (± strand) | `RIBOSEQC_ANALYSIS.out.coverage` | `[meta, path(*_coverage_*.bedgraph)]` |
| Unified ORF metadata | `UNIFY_ORF_PREDICTIONS.out.metadata` | `path(*.metadata.tsv)` |
| Unified ORF BED | `UNIFY_ORF_PREDICTIONS.out.bed` | `path(*.bed.gz)` |
| ORF confidence | `ORF_QC.out.confidence` | `path(*_orf_confidence.tsv)` |
| ORF classification | `CLASSIFY_ORFS_*` | `gencode_results.orfs.out`, etc. |

### 2.2 RiboseQC Bedgraph 文件格式

RiboseQC 为**每个样本**生成 4 个 bedgraph 文件:
- `{sample}_P_sites_plus.bedgraph` — P-site 位置 (正链)
- `{sample}_P_sites_minus.bedgraph` — P-site 位置 (负链)
- `{sample}_coverage_plus.bedgraph` — read coverage (正链, RPM 归一化)
- `{sample}_coverage_minus.bedgraph` — read coverage (负链, RPM 归一化)

格式: `chrom \t start \t end \t value` (0-based, 半开区间)

### 2.3 现有 ORF QC 模块已提供的功能

OCS = 0.30·S_translation + 0.30·S_agreement + 0.20·S_coverage + 0.15·S_periodicity + 0.05·S_readlevel

ORF QC 模块已经在 **per-ORF** 层面计算了:
- `S_translation`: 翻译得分
- `S_coverage`: P-site 密度 (unique_psites / AA length)
- `S_periodicity`: 3-nt 周期性
- `n_detecting_tools`, `detecting_tools`

**但是没有做的是**: 回到原始 bedgraph 去提取每个样本的定量表达值 (reads、pN、RPKM)。

---

## 3. 集成方案

### 3.1 方案概述

建议分 **两个新模块** 集成，利用 pipeline 已有的 bedgraph 和 unified ORF 数据:

```
Module A: EXPRESSION_QUANTIFICATION  (核心, 必须集成)
Module B: ORF_GGRIBO_PLOTS           (可选, 按需启用)
```

### 3.2 Module A: ORF 表达量定量 (EXPRESSION_QUANTIFICATION)

#### 定位
在 ORF_QC 之后、MULTIQC 之前运行。输入是 unified ORF + RiboseQC bedgraphs。

#### 输入
```
ch_unified_meta          // unified ORF metadata (含 orf_id, chrom, start, end, strand)
ch_unified_bed           // unified ORF BED
ch_psites_bedgraph       // 所有样本的 P_sites bedgraph (来自 RiboseQC)
ch_coverage_bedgraph     // 所有样本的 coverage bedgraph (来自 RiboseQC)
ch_classification        // 可选: 分类结果 (gencode + orfquant + orftype)
ch_ocs                   // 可选: ORF confidence scores (来自 ORF_QC)
```

#### 处理流程
```
Phase 1: Bedgraph 索引
  - 对所有 P-site / coverage bedgraph 做 bgzip + tabix 索引
  - 或使用 awk 按坐标区间查询 (当前 post_analysis 做法)

Phase 2: Per-ORF 表达量提取
  - 对每个统一 ORF, 查询所有样本的 P-site bedgraph
  - 提取: total_reads, mean_coverage, max_coverage, pN
  - 染色体名标准化 (处理 P1_ 前缀等)

Phase 3: RPKM/TPM 计算
  - 从 coverage bedgraph 提取 per-ORF RPM 总和
  - RPKM = sum(RPM) / (orf_length / 1000)
  - TPM = RPKM / sum(all_RPKM) * 1e6

Phase 4: 输出
  - expression_summary.tsv : ORF × Sample 矩阵 (reads + pN)
  - expression_rpkm_tpm.tsv : ORF × Sample 矩阵 (RPKM + TPM)
  - expression_combined.tsv : 合并 OCS + 分类 + 表达量
```

#### 输出
| 文件 | 内容 | 用途 |
|------|------|------|
| `{prefix}_expression_summary.tsv` | per-ORF × per-sample reads + pN | 下游分析 |
| `{prefix}_expression_rpkm_tpm.tsv` | per-ORF × per-sample RPKM + TPM | 定量比较 |
| `{prefix}_orf_expression_combined.tsv` | 合并 OCS+tier+分类+表达量 | 完整 ORF 注释表 |

#### 参数
```
--skip_orf_expression        (default: false)
--orf_expression_min_ocs     (default: 0.0, 过滤低置信度 ORF)
--orf_expression_workers     (default: 4, 并行查询线程数)
```

#### 容器需求
- Python 3.9+ with pandas, numpy
- bgzip, tabix (可选优化) 或 纯 Python awk 调用

### 3.3 Module B: ggRibo 覆盖度可视化 (ORF_GGRIBO_PLOTS)

#### 定位
可选模块，在 EXPRESSION_QUANTIFICATION 之后运行。需要 R + ggRibo。

#### 输入
```
ch_expression_summary     // 来自 Module A
ch_psites_bedgraph        // P-site bedgraph
ch_gtf                    // 参考 GTF (或 unified ORF GTF)
ch_unified_bed
```

#### 处理流程
```
Phase 1: ORF 选择
  - 按 OCS 排序取 top N (默认 N=20)
  - 或按用户指定的 ORF ID 列表
  - 对每个选中的 ORF，找表达量最高的 top K 样本

Phase 2: ggRibo 绘图
  - 为每个 ORF 创建临时 GTF (单 exon)
  - 调用 ggRibo 绘制多样本覆盖图
  - 使用自定义 Range_info + gtf_import_custom 绕过命名空间冲突

Phase 3: 输出
  - per-ORF PNG 文件
  - 汇总 HTML (可选)
```

#### 输出
| 文件 | 内容 |
|------|------|
| `{prefix}_ggribo/*.png` | 每个代表性 ORF 的多样本 ggRibo 覆盖图 |
| `{prefix}_ggribo_summary.html` | 可翻页的汇总图集 (可选) |

#### 参数
```
--skip_orf_ggribo          (default: true, 默认不运行)
--orf_ggribo_top_n         (default: 20)
--orf_ggribo_samples_per   (default: 3)
--orf_ggribo_extend        (default: 200 bp)
--orf_ggribo_orf_ids       (指定 ORF ID 列表文件)
--orf_ggribo_container      (默认使用含 ggRibo 的容器)
```

#### 容器需求
- R 4.x with ggRibo, ggplot2, txdbmaker, GenomicFeatures
- 需要专门的容器 (ggRibo Bioconductor 安装较复杂)

---

## 4. 实现优先级与复杂度评估

| 模块 | 优先级 | 复杂度 | 理由 |
|------|--------|--------|------|
| Module A: 表达量定量 | **高** | 中 | 核心功能，数据已存在 pipeline 中，主要是查询和聚合逻辑 |
| Module B: ggRibo 可视化 | **低** | 高 | 需要专门 R 容器，安装复杂；适合后处理/交互式使用 |

### 4.1 Module A 实现细节

**核心改动:**

1. 新建 `modules/local/expression_quant/main.nf`
2. 新建 `bin/quantify_orf_expression.py` (基于 extract_amp_expression_smart.py 改造)
3. 新建 `bin/calc_orf_rpkm_tpm.py` (基于 calc_orf_rpkm.py 改造)

**关键改造点 (与 post_analysis 脚本的区别):**

| 方面 | post_analysis | pipeline 集成版 |
|------|--------------|----------------|
| ORF 输入 | AMP result.tsv (特定格式) | unified ORF metadata (标准格式) |
| 样本发现 | 扫描 file system 或 gencode 列名 | 从 Nextflow channel meta 获取 sample list |
| 输入路径 | 硬编码 `/home/25119231r/riboseq/run/...` | 通过 Nextflow staging 传入 |
| 染色体名 | 手动 strip `P{version}_` 前缀 | 自动匹配 (统一 ORF 已标准化) |
| 输出格式 | organism-specific 独立文件 | 统一前缀，与 QC/MultiQC 集成 |
| 并行策略 | Python ProcessPoolExecutor 或 OMP | 利用 Nextflow 的 per-sample 并行 + Python 内部并行 |

**与现有 ORF QC 的协同:**

ORF QC 已经有 `unique_psites` (唯一 P-site 位置数) 和 `length_aa` 这些字段，但没有 **per-sample 定量值**。Module A 补充的是: 样本级别的 reads count、pN、RPKM。

可以将 OCS + 表达量 合并为一个完整的 ORF 注释表:
```
orf_id | chrom | start | end | strand | orf_type | OCS | tier |
S_translation | S_agreement | S_coverage | S_periodicity | S_readlevel |
detecting_tools | n_detecting |
sample1_reads | sample1_pN | sample1_rpkm | sample1_tpm |
sample2_reads | sample2_pN | sample2_rpkm | sample2_tpm |
...
gencode_biotype | orfquant_category | orftype_class |
```

### 4.2 Module B 实现细节

**核心挑战:**

1. **ggRibo 依赖复杂**: R 包依赖 GenomicFeatures、txdbmaker、ggplot2、ggRibo (Bioconductor)，容器构建需时
2. **容器冲突**: ggRibo 内部依赖与 ORFquant/RiboseQC 容器可能有 BiocGenerics 命名空间冲突 (已知问题)
3. **临时 GTF 生成**: 需要为每个 ORF 创建最小 GTF，这是纯文本操作，可以放在 Python 预处理中

**建议**: Module B 保留为 **post-hoc 独立脚本**，通过 `bin/plot_orf_ggribo.R` 提供，不在 pipeline 中自动运行。用户可在 pipeline 完成后调用独立脚本。

---

## 5. 具体 Pipeline 改动清单

### 5.1 新增文件

```
modules/local/expression_quant/
├── main.nf                    # Nextflow process 定义
├── meta.yml                   # 模块元数据
└── environment.yml            # Conda 环境

bin/
├── quantify_orf_expression.py # ORF 表达量提取 (基于 post_analysis 改造)
└── calc_orf_rpkm_tpm.py       # RPKM/TPM 计算 (基于 post_analysis 改造)

bin/
└── plot_orf_ggribo.R          # (可选) ggRibo 独立脚本

containers/
└── Singularity.ggribo.def     # (可选) ggRibo 容器定义
```

### 5.2 修改文件

```
workflows/riboseq/main.nf      # 新增 EXPRESSION_QUANTIFICATION 调用
nextflow.config                # 新增参数声明
conf/modules.config            # 新增模块资源配置
nextflow_schema.json           # 新增参数 schema
```

### 5.3 参数新增

```
// Expression quantification
--skip_orf_expression          (boolean, default: false)
--orf_expression_min_ocs       (number,  default: 0.0)
--orf_expression_workers       (integer, default: 4)

// ggRibo plots (optional)
--skip_orf_ggribo              (boolean, default: true)
--orf_ggribo_top_n             (integer, default: 20)
--orf_ggribo_samples_per       (integer, default: 3)
--orf_ggribo_orf_ids           (string,  path to ID list)
--orf_ggribo_container          (string)
```

### 5.4 工作流改动 (workflows/riboseq/main.nf)

在 ORF_QC 和 MULTIQC 之间插入:

```nextflow
//
// ORF Expression Quantification
// Extracts per-sample reads, pN, RPKM, TPM from RiboseQC bedgraphs for each unified ORF.
//
if (!params.skip_orf_expression) {
    // Collect per-sample bedgraph files
    ch_psites_bedgraph = RIBOSEQC_POSTFILTER.out.psites_bedgraph
        .map { meta, f -> f }
        .collect()

    ch_coverage_bedgraph = RIBOSEQC_POSTFILTER.out.coverage
        .map { meta, f -> f }
        .collect()

    EXPRESSION_QUANTIFICATION(
        ch_unify_metadata,
        ch_unify_bed,
        ch_psites_bedgraph,
        ch_coverage_bedgraph,
        ch_ocs,           // optional: from ORF_QC
        ch_classification  // optional: from CLASSIFY_ORFS
    )
}
```

---

## 6. 风险评估

| 风险 | 影响 | 缓解措施 |
|------|------|---------|
| Bedgraph 文件数量大 (n_samples × 4) | I/O 压力 | bgzip + tabix 索引; 流式处理 |
| awk 查询可能因染色体名不匹配返回空 | 漏数据 | 统一染色体名标准化逻辑 |
| coverage bedgraph 已经是 RPM，但假设不同 | RPKM 偏倚 | 在文档中注明假设; 可选保留 raw count |
| ggRibo 容器 BiocGenerics 冲突 | 运行时错误 | 自定义 Range_info 类绕过 (已在 post_analysis 中验证) |
| 大规模 ORF 集 (>500k) — 查询慢 | 运行时间过长 | 默认按 OCS 过滤; 支持并行 worker |

---

## 7. 建议实施步骤

1. **Phase 1 (1-2 天)**: 创建 `quantify_orf_expression.py`，基于 post_analysis 脚本改造，统一输入输出格式
2. **Phase 2 (0.5 天)**: 创建 `modules/local/expression_quant/main.nf`
3. **Phase 3 (0.5 天)**: 修改 `workflows/riboseq/main.nf` 和 `nextflow.config`
4. **Phase 4 (0.5 天)**: 测试 + nf-test
5. **Phase 5 (可选)**: ggRibo 独立脚本 + 容器
