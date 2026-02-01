# ORF 检测工具评分机制说明

本文档说明 ORFquant、Ribo-TISH 和 Ribotricer 三个工具如何对检测到的 ORF 进行质量评估。

---

## 1. Ribo-TISH (Translation Initiation Site Hallmark)

### 评分指标：**TisPvalue**
- **含义**：翻译起始位点（TIS）的统计显著性 p 值
- **范围**：0-1（通常为科学计数法，如 1.2e-10）
- **解读**：
  - **越小越好**（更显著）
  - p < 0.05 通常认为显著
  - p < 0.001 表示高置信度
- **统计检验**：基于 Ribo-seq 信号在 TIS 位置的三周期性（3-nt periodicity）

### unify_orf_predictions.py 中的处理
```python
score = -log10(TisPvalue)  # 转换为 -log10(p)，值越大越好
```
- **转换原因**：统一评分方向（所有工具都是分数越高越好）
- **示例**：
  - TisPvalue = 1.2e-10 → score ≈ 9.92
  - TisPvalue = 0.001 → score = 3.0
  - TisPvalue = 0.05 → score ≈ 1.30

### 其他输出信息
- **TisType**：起始密码子类型（ATG、CTG、GTG 等）
- **TisGroup**：ORF 分类（annotated、uORF、dORF、extension 等）
- **RiboPvalue**：核糖体密度的显著性（用于进一步验证）

---

## 2. Ribotricer

### 评分指标：**phase_score**
- **含义**：Ribo-seq 读段相位得分，衡量三周期性的强度
- **范围**：0-1
- **解读**：
  - **越大越好**
  - 0.95+：极高置信度（强三周期性）
  - 0.80-0.95：高置信度
  - 0.60-0.80：中等置信度
  - < 0.60：低置信度
- **计算方法**：基于 ORF 内读段的阅读框分布熵

### unify_orf_predictions.py 中的处理
```python
score = phase_score  # 直接使用（0-1，越大越好）
```

### 其他输出信息
- **read_count**：该 ORF 的总读段数
- **read_density**：读段密度（reads per nucleotide）
- **valid_codons**：有效密码子数量
- **ORF_type**：ORF 类型（annotated、uORF、dORF、internal 等）
- **status**：翻译状态（translating、non_translating）

---

## 3. ORFquant

### 评分指标：**P_sites**
- **含义**：ORF 内检测到的 P-site（核糖体 P 位点）数量
- **范围**：整数（通常 10-1000+）
- **解读**：
  - **越大越好**（更多覆盖度）
  - 23+：典型的高表达 ORF
  - 10-23：中等表达
  - < 10：低表达或假阳性风险
- **生物学意义**：直接反映核糖体占用量（ribosome occupancy）

### unify_orf_predictions.py 中的处理
```python
score = P_sites  # 直接使用（整数，越大越好）
```

### 其他输出信息（GTF 属性）
- **ORF_pct_P_sites**：该 ORF 占转录本总 P-site 的百分比
- **ORF_pct_P_sites_pN**：标准化后的百分比（考虑长度）
- **ORFs_pM**：每百万映射读段的 ORF 数（标准化表达量）

---

## 工具对比总结

| 工具 | 评分指标 | 数值范围 | 方向 | 生物学基础 | 推荐阈值 |
|------|---------|---------|------|------------|----------|
| **Ribo-TISH** | -log10(TisPvalue) | 0-∞ (通常 0-15) | 越大越好 | TIS 统计显著性 | > 3.0 (p<0.001) |
| **Ribotricer** | phase_score | 0-1 | 越大越好 | 三周期性强度 | > 0.80 |
| **ORFquant** | P_sites | 整数 (10-1000+) | 越大越好 | 核糖体占用量 | > 10 |

---

## 综合评估建议

### 高置信度 ORF（推荐用于下游分析）
满足以下**任一**条件：
1. **多工具一致**：被 2+ 工具检测到
2. **Ribo-TISH**：-log10(TisPvalue) > 5.0 (p < 1e-5)
3. **Ribotricer**：phase_score > 0.90
4. **ORFquant**：P_sites > 50

### 中等置信度 ORF（需进一步验证）
- 单工具检测，但评分较高
- 多工具检测，但评分中等

### 低置信度 ORF（可能假阳性）
- 仅单工具检测且评分低
- Ribotricer phase_score < 0.60
- ORFquant P_sites < 5

---

## 使用示例

### 在 metadata.tsv 中筛选高质量 ORF

```bash
# 筛选多工具检测的 ORF
awk -F'\t' '$10 ~ /,/' unified.metadata.tsv > high_confidence.tsv

# 解析 tool_scores 列（示例：ribotish:5.2,orfquant:45）
# Ribo-TISH 高分 ORF
awk -F'\t' '$12 ~ /ribotish:/ && $12 ~ /:[5-9]\./ || $12 ~ /:[0-9]{2}/' unified.metadata.tsv

# ORFquant 高 P-sites ORF
awk -F'\t' '$12 ~ /orfquant:[5-9][0-9]/ || $12 ~ /orfquant:[0-9]{3}/' unified.metadata.tsv
```

### 在 GTF 中筛选高质量 ORF

```bash
# 多工具检测（num_tools >= 2）
awk '$0 ~ /num_tools "2"/ || $0 ~ /num_tools "3"/' unified.gtf

# 包含特定工具
grep 'sources ".*ribotish.*"' unified.gtf
```

---

## 参考文献

1. **Ribo-TISH**: Zhang et al., *Nucleic Acids Research*, 2017
   - DOI: 10.1093/nar/gkx1206
   
2. **Ribotricer**: Choudhary et al., *Nature Communications*, 2020
   - DOI: 10.1038/s41467-020-17634-8

3. **ORFquant**: Calviello et al., *Nucleic Acids Research*, 2016
   - DOI: 10.1093/nar/gkv1402

---

*最后更新：2026-02-01*
