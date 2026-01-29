# GENCODE 转换器脚本

## 概述

此目录包含将 Ribo-seq ORF 预测工具输出转换为 gencode-riboseqORFs 兼容格式的 Python 脚本。

## 可用转换器

### ✅ ribotish_to_gencode.py

**状态**: 完成并测试（v1.2 - 已修复序列提取）

将 Ribo-TISH 预测结果转换为 GENCODE 格式。

**使用**:
```bash
python3 ribotish_to_gencode.py \
    --predict ribotish_predict.txt \
    --fasta genome.fa \
    --study_id SAMPLE1 \
    --output_prefix output \
    --min_length 16
```

**最近修复** (2026-01-20 v1.2):
- ✅ **关键修复**: 根据 ORF 长度提取序列，不再使用基因组坐标全区间
- ✅ 正确处理多个 stop codon（截断到第一个）
- ✅ BED 文件使用唯一 ORF 名称（`GENE_START_LENaa`）
- ✅ 自动调整序列长度确保一致性
- ✅ BED 和 FASTA 格式完全匹配

### ✅ ribotricer_to_gencode.py

**状态**: 完成并测试

将 Ribotricer 预测结果转换为 GENCODE 格式。

**使用**:
```bash
python3 ribotricer_to_gencode.py \
    --tsv sample_translating_ORFs.tsv \
    --fasta genome.fa \
    --study_id SAMPLE1 \
    --output_prefix output \
    --min_length 16 \
    --min_phase_score 0.5
```

### 🚧 未来计划

- ribocode_to_gencode.py - RiboCode 转换器
- rpbp_to_gencode.py - rp-bp 转换器  
- orfquant_to_gencode.py - ORFquant 转换器

## 输出格式

所有转换器生成两个文件：

### 1. FASTA 文件 (`{prefix}.gencode.fa`)

```
>{GENE}_{START}_{LENGTH}aa--{STUDY_ID}
SEQUENCE*
```

- Header 格式: `{ORF_NAME}--{STUDY_ID}`
- 序列必须以 stop codon `*` 结尾
- 长度（不含 `*`）必须与 header 声明一致

**示例**:
```
>ENSMUST00000156816_4846686_3071aa--S1
MXXXXXXXXXXXXXXXXXXXXXXX...XXXXXXXX*
```

### 2. BED 文件 (`{prefix}.gencode.bed`)

```
chr  start  end  ORF_NAME  STUDY_ID  strand
```

- **重要**: 使用 1-based 坐标（gencode-riboseqORFs 要求）
- Tab 分隔，6 列
- 无 header 行

**示例**:
```
1  4846686  4855900  ENSMUST00000156816_4846686_3071aa  S1  -
1  4846686  4846968  ENSMUST00000156816_4846686_94aa   S1  -
```

## 验证工具

### validate_output.py - 详细验证

完整的格式验证和一致性检查。

**使用**:
```bash
python3 validate_output.py \
    --bed output.gencode.bed \
    --fasta output.gencode.fa
```

**检查项**:
- ✅ BED 格式（6 列）
- ✅ ORF 名称格式和唯一性
- ✅ 坐标有效性
- ✅ FASTA header 格式
- ✅ 序列长度一致性
- ✅ Stop codon 存在
- ✅ BED 和 FASTA 交叉验证

### quick_check.sh - 快速诊断

快速的 shell 脚本诊断。

**使用**:
```bash
bash quick_check.sh \
    --bed output.gencode.bed \
    --fasta output.gencode.fa
```

**输出**:
```
==================================================
  Quick Check: GENCODE Converter Output
==================================================

📄 BED File: output.gencode.bed
   Total lines: 456
   ✅ All lines have 6 columns
   ✅ All ORF names follow GENE_START_LENaa format
   ✅ No duplicate ORF names

📄 FASTA File: output.gencode.fa
   Total ORFs: 456
   ✅ All headers follow GENE_START_LENaa--STUDY format
   ✅ All sequences end with stop codon (*)
   ✅ All sequences match their declared lengths

🔗 Cross-Validation
   ✅ BED and FASTA have matching ORF sets

==================================================
✅ Quick check PASSED
```

## 使用工作流

### 1. 运行 ORF 预测工具

```bash
# Ribo-TISH
ribotish predict -b sample.bam -g annotation.gtf -f genome.fa \
    -p sample.para.py -o sample_pred.txt

# Ribotricer
ribotricer detect-orfs -b sample.bam -r ribotricer_index \
    -o sample_translating_ORFs.tsv
```

### 2. 转换为 GENCODE 格式

```bash
# Ribo-TISH
python3 ribotish_to_gencode.py \
    --predict sample_pred.txt \
    --fasta genome.fa \
    --study_id S1 \
    --output_prefix S1

# Ribotricer
python3 ribotricer_to_gencode.py \
    --tsv sample_translating_ORFs.tsv \
    --fasta genome.fa \
    --study_id S1 \
    --output_prefix S1
```

### 3. 验证输出

```bash
# 快速检查
bash quick_check.sh --bed S1.gencode.bed --fasta S1.gencode.fa

# 详细验证
python3 validate_output.py --bed S1.gencode.bed --fasta S1.gencode.fa
```

### 4. 合并多个样本/工具

```bash
# 合并 FASTA
cat S1_ribotish.gencode.fa S1_ribotricer.gencode.fa \
    S2_ribotish.gencode.fa S2_ribotricer.gencode.fa \
    > all_orfs.fa

# 合并 BED
cat S1_ribotish.gencode.bed S1_ribotricer.gencode.bed \
    S2_ribotish.gencode.bed S2_ribotricer.gencode.bed \
    > all_orfs.bed
```

### 5. 运行 gencode-riboseqORFs

```bash
# 准备 Ensembl 注释（仅需一次）
cd gencode-riboseqORFs
bash scripts/retrieve_ensembl_data.sh 110 GRCm39

# 运行 ORF mapper
python3 ORF_mapper_to_GENCODE_v1.1.py \
    -d Ens110_GRCm39/ \
    -f all_orfs.fa \
    -b all_orfs.bed \
    -o project_name \
    -l 16 \
    -c 0.9 \
    -m longest_string
```

## 依赖项

### 必需
- Python 3.9+
- Biopython 1.79+

### 可选（用于序列提取）
- pyfaidx 0.7.2+ （推荐）
- bedtools 2.30.0+ （备用）

**安装**:
```bash
pip install biopython pyfaidx

# 或使用 conda
conda install -c bioconda biopython pyfaidx
```

## 坐标系统说明

**关键**: 不同工具使用不同的坐标系统！

| 工具 | 输入坐标 | 输出 BED | 转换 |
|------|---------|----------|------|
| Ribo-TISH | 1-based | 1-based | ✅ 无需转换 |
| Ribotricer | 1-based | 1-based | ✅ 无需转换 |
| RiboCode | 0-based | 1-based | ⚠️ 需要 +1 |
| rp-bp | 0-based | 1-based | ⚠️ 需要 +1 |
| ORFquant | GTF (1-based) | 1-based | ✅ 无需转换 |

## 常见问题

### Q1: 序列全是 'M'
**A**: 序列提取失败，检查：
- 基因组 FASTA 是否存在且完整
- 染色体名称是否匹配（chr1 vs 1）
- pyfaidx 或 bedtools 是否正确安装

### Q2: "Adjusted ORF length" 警告
**A**: 正常现象。工具报告的长度有时不精确，脚本会自动调整为实际翻译长度。

### Q3: BED 和 FASTA 数量不一致
**A**: 可能原因：
- 某些 ORF 翻译失败（检查 stderr）
- 序列提取失败
- 运行 `validate_output.py` 查看详情

### Q4: 重复的 ORF 名称
**A**: 不应该发生。检查：
- 上游预测工具输出是否有重复
- 可能需要去重

## 文档

- [测试指南](TESTING_GUIDE.md) - 详细的测试步骤
- [修复说明](FIX_NOTES.md) - v1.1 修复详情
- [转换器总览](../../bin/README_CONVERTERS.md) - 独立脚本文档

## 开发状态

| 组件 | 状态 | 优先级 |
|------|------|--------|
| ribotish_to_gencode.py | ✅ v1.2 完成 | 高 |
| ribotricer_to_gencode.py | ✅ v1.0 完成 | 高 |
| validate_output.py | ✅ 完成 | 高 |
| quick_check.sh | ✅ 完成 | 中 |
| 测试数据 | ✅ 完成 | 高 |
| ribocode_to_gencode.py | ⏳ 计划中 | 中 |
| rpbp_to_gencode.py | ⏳ 计划中 | 中 |
| orfquant_to_gencode.py | ⏳ 计划中 | 低 |

## 相关资源

- [gencode-riboseqORFs GitHub](https://github.com/jorruior/gencode-riboseqORFs)
- [Ribo-TISH](https://github.com/zhpn1024/ribotish)
- [Ribotricer](https://github.com/smithlabcode/ribotricer)
- [nf-core/riboseq](https://nf-co.re/riboseq)

## 贡献

发现 bug 或有改进建议？

1. 创建 issue 描述问题
2. 提供示例输入/输出
3. 包含错误日志

## 许可

MIT License - 详见项目根目录 LICENSE 文件

---

**最后更新**: 2026-01-20  
**维护者**: nf-core/riboseq 团队  
**版本**: v1.2
