# Ribo-seq ORF 检测工具统计模型汇总

本文档汇总 Ribo-TISH、Ribotricer 和 ORFquant 三个工具的统计检验模型、评分指标和整合方案。

---

## 1. Ribo-TISH (Translation Initiation Site Hallmark)

### 统计检验模型

**核心方法**：**Chi-square (χ²) 检验** + **Binomial 检验**

#### 1.1 TisPvalue (翻译起始位点显著性)
- **检验类型**：Chi-square goodness-of-fit test
- **零假设**：读段在起始密码子附近的分布是随机的
- **备择假设**：读段在起始密码子位置有显著富集
- **检验对象**：P-site 在 TIS 上下游窗口的分布模式
- **公式**：
  ```
  χ² = Σ[(Observed - Expected)² / Expected]
  TisPvalue = P(χ² > observed_χ²)
  ```

#### 1.2 RiboPvalue (核糖体密度显著性)
- **检验类型**：Binomial test
- **零假设**：ORF 区域的核糖体覆盖度与背景一致
- **备择假设**：ORF 区域有显著更高的核糖体占用
- **公式**：
  ```
  P(X ≥ k) where X ~ Binomial(n, p_background)
  ```

#### 1.3 RiboPval_adj (FDR 校正)
- **校正方法**：Benjamini-Hochberg FDR
- **作用**：多重检验校正，控制假发现率

### 输出统计量
| 字段 | 含义 | 类型 | 范围 |
|------|------|------|------|
| **TisPvalue** | TIS 位点显著性 | p-value | 0-1 |
| **RiboPvalue** | 核糖体密度显著性 | p-value | 0-1 |
| **RiboPval_adj** | FDR 校正后的 p-value | adjusted p | 0-1 |
| **TisCount** | TIS 位置的读段数 | integer | ≥0 |
| **TisLen** | ORF 长度（nt） | integer | ≥0 |

### 参考文献
- Zhang et al., *Nucleic Acids Research*, 2017
- DOI: 10.1093/nar/gkx1206

---

## 2. Ribotricer

### 统计检验模型

**核心方法**：**Phase Score（基于信息熵）** + **无参数检验**

#### 2.1 Phase Score (三周期性强度)
- **检验类型**：基于 Shannon entropy 的归一化指标
- **原理**：真实翻译的 ORF 应该有强三周期性（reads 主要落在 frame 0）
- **公式**：
  ```
  H = -Σ[p_i * log2(p_i)]  # Shannon entropy
  H_max = log2(3) = 1.585  # 最大熵（均匀分布）
  Phase_score = 1 - (H / H_max)
  
  其中 p_i = frame_i_count / total_count (i=0,1,2)
  ```
- **取值范围**：0-1
  - 1.0：完美三周期性（所有读段在 frame 0）
  - 0.33：随机分布（三个 frame 均等）
  - 0.0：反三周期性（frame 0 无读段）

#### 2.2 Read Count & Density
- **read_count**：ORF 内的总读段数（所有 frame）
- **read_density**：reads per nucleotide
- **无统计检验**：仅提供描述性统计，无 p-value

#### 2.3 Valid Codons Ratio
- **valid_codons_ratio**：有读段覆盖的密码子比例
- **作用**：评估 ORF 的完整性

### 输出统计量
| 字段 | 含义 | 类型 | 范围 |
|------|------|------|------|
| **phase_score** | 三周期性得分 | float | 0-1 |
| **read_count** | ORF 总读段数 | integer | ≥0 |
| **read_density** | 读段密度 | float | ≥0 |
| **valid_codons** | 有效密码子数 | integer | ≥0 |
| **valid_codons_ratio** | 有效密码子比例 | float | 0-1 |
| **length** | ORF 长度（nt） | integer | ≥0 |

### 参考文献
- Choudhary et al., *Nature Communications*, 2020
- DOI: 10.1038/s41467-020-17634-8

---

## 3. ORFquant

### 统计检验模型

**核心方法**：**基于 P-site 计数的描述性统计** + **标准化表达量**

#### 3.1 P-sites (核糖体 P 位点计数)
- **检验类型**：无统计检验，直接计数
- **定义**：ORF 内检测到的 P-site 数量（由 RiboseQC 提供）
- **生物学意义**：直接反映核糖体占用量
- **依赖**：RiboseQC 的 P-site 偏移校准

#### 3.2 ORF_pct_P_sites (ORF 在转录本中的 P-site 占比)
- **公式**：
  ```
  ORF_pct_P_sites = (ORF_P_sites / Transcript_total_P_sites) × 100
  ```
- **作用**：识别主要翻译的 ORF（高占比 → 主 ORF）

#### 3.3 ORF_pct_P_sites_pN (长度标准化的占比)
- **公式**：
  ```
  ORF_pct_P_sites_pN = (ORF_P_sites / ORF_length) / (Transcript_total_P_sites / Transcript_length) × 100
  ```
- **作用**：消除长度偏差，比较不同长度 ORF 的翻译效率

#### 3.4 ORFs_pM (每百万映射读段的 ORF 表达量)
- **公式**：
  ```
  ORFs_pM = (ORF_P_sites / Total_mapped_reads) × 10^6
  ```
- **作用**：跨样本标准化表达量（类似 RPM）

### 输出统计量
| 字段 | 含义 | 类型 | 范围 |
|------|------|------|------|
| **P_sites** | ORF 内 P-site 数量 | integer | ≥0 |
| **ORF_pct_P_sites** | P-site 占转录本比例 | float | 0-100 |
| **ORF_pct_P_sites_pN** | 长度标准化占比 | float | 0-100+ |
| **ORFs_pM** | 每百万读段表达量 | float | ≥0 |

### 参考文献
- Calviello et al., *Nucleic Acids Research*, 2016
- DOI: 10.1093/nar/gkv1402

---

## 4. RiboseQC 数据整合

### RiboseQC 提供的统计信息

#### 4.1 Bedgraph 文件（P-site 位置计数）
- **文件类型**：
  - `*_P_sites_plus.bedgraph` / `*_P_sites_minus.bedgraph`：所有读段
  - `*_P_sites_uniq_plus.bedgraph` / `*_P_sites_uniq_minus.bedgraph`：唯一映射读段
  - `*_coverage_*.bedgraph`：覆盖度（非 P-site）

#### 4.2 可提取的统计量
从 bedgraph 文件可计算：
1. **Total P-sites**：ORF 内所有 P-site 计数总和
2. **Unique P-sites**：仅唯一映射的 P-site 计数
3. **pN (P-sites per Nucleotide)**：
   ```
   pN = Total_P_sites / ORF_length_nt
   ```
4. **Unique pN**：
   ```
   unique_pN = Unique_P_sites / ORF_length_nt
   ```

#### 4.3 与 ORFquant 的关系
- ORFquant 的 `P_sites` 字段 = RiboseQC 统计的 P-site 总数
- ORFquant 依赖 RiboseQC 的 `*_for_ORFquant` 文件

---

## 5. 统一统计框架设计

### 5.1 需要整合的统计量

#### 核心指标（所有工具都提供或可计算）
| 统计量 | Ribo-TISH | Ribotricer | ORFquant | RiboseQC |
|--------|-----------|------------|----------|----------|
| **ORF 长度 (nt)** | ✅ TisLen | ✅ length | ✅ (从 GTF) | ✅ |
| **总读段数** | ✅ TisCount | ✅ read_count | ❌ | ✅ (coverage) |
| **唯一读段数** | ❌ | ❌ | ❌ | ✅ (uniq coverage) |
| **P-site 数** | ❌ | ❌ | ✅ P_sites | ✅ (bedgraph) |
| **唯一 P-site 数** | ❌ | ❌ | ❌ | ✅ (uniq bedgraph) |
| **pN 值** | ❌ | ✅ read_density | ✅ (可计算) | ✅ (可计算) |
| **统计显著性** | ✅ TisPvalue | ❌ | ❌ | ❌ |

### 5.2 整合策略

#### 策略 A：直接从工具输出提取（当前实现）
- **优点**：简单快速
- **缺点**：不同工具统计口径不一致
  - Ribo-TISH 的 TisCount 可能包含多映射读段
  - Ribotricer 的 read_count 可能与实际 P-site 不对应
  - ORFquant 仅提供 P-site，无总读段数

#### 策略 B：统一从 RiboseQC bedgraph 重新计算（推荐）⭐
- **优点**：
  - 统一数据源，确保一致性
  - 可同时提供 total 和 unique 计数
  - 支持 P-site 和 coverage 两种模式
- **缺点**：需要额外计算步骤

#### 策略 C：混合策略
- **P-site 相关**：从 RiboseQC 统一计算
- **统计显著性**：从各工具原始输出提取
- **工具特异性指标**：保留（如 phase_score、ORF_pct_P_sites）

---

## 6. 实施方案（推荐）

### 6.1 扩展 `unify_orf_predictions.py` 输出字段

#### 新增统计列（metadata.tsv）
```
orf_id, chrom, strand, start, end, length_aa, exon_blocks, gene_id, transcript_id, 
tools, samples, tool_scores, tool_pvalues,  # 新增 tool_pvalues
total_reads, unique_reads,                    # 新增读段计数
total_psites, unique_psites,                  # 新增 P-site 计数
pN, unique_pN,                                # 新增密度
sequence
```

### 6.2 从 RiboseQC bedgraph 计算统计量

#### 输入
- RiboseQC P-site bedgraph 文件（每个样本）
  - `{sample}_P_sites_plus.bedgraph`
  - `{sample}_P_sites_minus.bedgraph`
  - `{sample}_P_sites_uniq_plus.bedgraph`
  - `{sample}_P_sites_uniq_minus.bedgraph`

#### 计算流程
```python
def calculate_orf_statistics(orf_candidate, bedgraph_dir, samples):
    """
    为每个 ORF 计算统计量
    
    Returns:
        {
            'total_psites': int,      # 所有样本总和
            'unique_psites': int,
            'total_reads': int,       # 从 coverage bedgraph
            'unique_reads': int,
            'pN': float,              # total_psites / length_nt
            'unique_pN': float
        }
    """
    total_psites = 0
    unique_psites = 0
    
    for sample in samples:
        # 读取对应 strand 的 bedgraph
        bedgraph_file = f"{bedgraph_dir}/{sample}_P_sites_{strand}.bedgraph"
        uniq_bedgraph = f"{bedgraph_dir}/{sample}_P_sites_uniq_{strand}.bedgraph"
        
        # 计算 ORF 区域内的 P-site 计数（需处理多 exon）
        for block_start, block_end in orf_candidate.blocks:
            total_psites += count_psites_in_region(bedgraph_file, chrom, block_start, block_end)
            unique_psites += count_psites_in_region(uniq_bedgraph, chrom, block_start, block_end)
    
    pN = total_psites / orf_candidate.length_nt
    unique_pN = unique_psites / orf_candidate.length_nt
    
    return {...}
```

### 6.3 保留工具原始 p-value

#### 修改 `tool_scores` 并行存储 `tool_pvalues`
```python
# 在 ORFCandidate 类中
self.tool_scores = {tool: score}
self.tool_pvalues = {tool: pvalue}  # 新增

# metadata.tsv 输出示例
tool_scores: ribotish:9.92,orfquant:23,ribotricer:0.95
tool_pvalues: ribotish:1.2e-10,ribotish_ribo:3.5e-08  # 保留原始 p 值
```

### 6.4 新增命令行参数

```python
parser.add_argument("--bedgraph-dir", nargs='+', 
                    help="RiboseQC bedgraph directories (one per sample)")
parser.add_argument("--calculate-stats", action='store_true',
                    help="Calculate unified statistics from bedgraphs")
```

---

## 7. 输出示例

### metadata.tsv (扩展版)
```tsv
orf_id	chrom	strand	start	end	length_aa	exon_blocks	gene_id	transcript_id	tools	samples	tool_scores	tool_pvalues	total_reads	unique_reads	total_psites	unique_psites	pN	unique_pN	sequence
ORF_1_ENSG001	chr1	+	100000	100333	111	100000-100333	ENSG001	ENST001	ribotish,orfquant	sample1,sample2	ribotish:9.92,orfquant:45	ribotish:1.2e-10,ribotish_ribo:3.5e-08	150	142	45	42	0.135	0.126	ATGAAA...
ORF_2_ENSG002	chr1	-	200000	200150	50	200000-200150	ENSG002	ENST002	ribotricer	sample1	ribotricer:0.95	NA	80	75	38	35	0.253	0.233	ATGCCC...
```

### 字段说明
- **total_reads**：从 coverage bedgraph 计算（可选）
- **unique_reads**：唯一映射读段数
- **total_psites**：P-site 总数（推荐主要指标）
- **unique_psites**：唯一 P-site 数
- **pN**：P-sites per nucleotide（翻译强度）
- **unique_pN**：基于唯一映射的 pN（更可靠）
- **tool_pvalues**：各工具原始 p 值（Ribo-TISH 提供）

---

## 8. 使用建议

### 8.1 质量筛选阈值

#### 基于统计显著性
```python
# 高置信度 ORF
(tool_pvalues.ribotish < 1e-5) OR (num_tools >= 2)
```

#### 基于 P-site 密度
```python
# 高翻译活性 ORF
unique_pN > 0.1 AND unique_psites >= 10
```

#### 基于工具一致性
```python
# 多工具验证
num_tools >= 2 AND (
    ribotricer.phase_score > 0.8 OR 
    ribotish.pvalue < 0.001
)
```

### 8.2 下游分析

#### 差异翻译分析
使用 `unique_psites` 构建计数矩阵，输入 DESeq2/edgeR

#### ORF 分类
结合 `pN` 值和 `ORF_pct_P_sites` 区分主 ORF、uORF、dORF

#### 翻译效率计算
```
TE = unique_pN_ORF / mRNA_RPKM
```

---

## 9. 实施优先级

### Phase 1（立即实施）✅
1. 保留原始 p-value（不转换）
2. 添加 `tool_pvalues` 列
3. 文档更新（本文件）

### Phase 2（推荐）⭐
4. 实现从 RiboseQC bedgraph 计算统计量
5. 添加 `total_psites`, `unique_psites`, `pN`, `unique_pN` 列
6. 提供 `--bedgraph-dir` 参数

### Phase 3（可选）
7. 计算 `total_reads`, `unique_reads`（从 coverage bedgraph）
8. 添加样本级别的统计量（per-sample breakdown）
9. 整合到 Nextflow pipeline

---

## 10. 参考资料

### 工具文档
- Ribo-TISH: https://github.com/zhpn1024/ribotish
- Ribotricer: https://github.com/smithlabcode/ribotricer
- ORFquant: https://github.com/lcalviell/ORFquant
- RiboseQC: https://github.com/ohlerlab/RiboseQC

### 相关文献
1. Zhang et al. (2017) Ribo-TISH. *NAR* 45(11):e103
2. Choudhary et al. (2020) Ribotricer. *Nat Commun* 11:3826
3. Calviello et al. (2016) ORFquant. *NAR* 44(4):e29
4. Calviello et al. (2021) RiboseQC. *BMC Bioinformatics* 22:342

---

*创建日期：2026-02-01*
*作者：AI Assistant*
