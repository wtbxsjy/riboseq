# 单工具 Singularity 验证脚本（从 BAM 开始）

这组脚本用于在跑整条 Nextflow 流程前，先逐个验证每个工具在你的环境里能用（镜像版本与 pipeline module 保持一致）。

## 通用要求

- 需要 `singularity` 可用（或 Apptainer 兼容命令）。
- 建议在一个"干净的工作目录"里运行脚本，输出会写到当前目录下的子目录。
- 大多数脚本要求输入 BAM **已排序并已建立索引**（`.bai` 或 `.csi`）。脚本里对常见情况会自动补 `samtools index` / `samtools faidx`。
- 多线程：
  - `--cpus` 会用于 `samtools index`（以及 BAM 过滤脚本中的 `samtools view`），这类步骤通常能明显加速。
  - 其它工具是否真正多线程取决于工具本身；如果脚本没有显式线程参数，可用 `--args` 透传（例如 ribotricer/ribocode 的额外参数）。

## WSL / OneDrive 路径（可选）

如果你的数据在 `/mnt/c/...`，有时需要显式 bind：

```bash
export BIND_EXTRA="/mnt:/mnt"
```

所有脚本都支持 `BIND_EXTRA`（逗号分隔多个 bind 路径）。

## 容器缓存

脚本默认把镜像 pull 到 `scripts/singularity_single_tool_tests/containers/`，避免重复下载。

## 脚本列表

### ORF 预测工具（01-13）

- `00_run_order_example.sh`：把 01..11 的示例按顺序串联（带大量注释；默认 `DRY_RUN=1` 只打印不执行）。
- `01_sorf_bam_filter.sh`：sORF 预测前 BAM 过滤（unique / contig regex / read length / flags / 去 duplicate）。
- `02_riboseqc_prepareannotation.sh`：RiboseQC 准备 annotation（`*_Rannot`）。
- `03_riboseqc_analysis.sh`：RiboseQC 分析（产出 `*_for_ORFquant` 等）。
- `04_orfquant_run.sh`：ORFquant（读取 `*_for_ORFquant` + `*_Rannot`）。
- `05_ribotish_quality.sh`：RiboTISH quality（产 offset 参数 `*.para.py`）。
- `06_ribotish_predict.sh`：RiboTISH predict（per-sample）。
- `07_ribotricer_prepareorfs.sh`：Ribotricer prepare-orfs（候选 ORFs index）。
- `08_ribotricer_detectorfs.sh`：Ribotricer detect-orfs（per-sample）。
- `09_rpbp_prepare_genome.sh`：rp-bp prepare-rpbp-genome。
- `10_rpbp_predict.sh`：rp-bp predict（per-sample）。
- `11_ribocode_detect.sh`：RiboCode onestep（注意：更推荐 transcriptome BAM；低深度数据可能失败属正常）。
- `12_riboseqc_orfquant_from_filtered_bam.sh`：从过滤后 BAM 运行 RiboseQC + ORFquant。
- `13_orfquant_prepareannotation.sh`：ORFquant 注释准备。

### GENCODE 注释工具（14-16）✨ 新增

- **`14_ribotish_to_gencode.sh`**：将 Ribo-TISH 预测结果转换为 gencode-riboseqORFs 兼容格式
- **`15_ribotricer_to_gencode.sh`**：将 Ribotricer 预测结果转换为 gencode-riboseqORFs 兼容格式
- **`16_gencode_orf_mapper.sh`**：运行 gencode-riboseqORFs 统一 ORF 注释和去冗余

## 最小运行示例

### 基础 ORF 预测

```bash
cd scripts/singularity_single_tool_tests

# 1) 过滤
./01_sorf_bam_filter.sh \
  --sample S1 \
  --bam /path/to/S1.bam \
  --fai /path/to/genome.fa.fai

# 2) RiboTISH quality + predict
./05_ribotish_quality.sh --sample S1 --bam /path/to/S1.filtered.bam --gtf /path/to/annot.gtf
./06_ribotish_predict.sh --sample S1 --bam /path/to/S1.filtered.bam --gtf /path/to/annot.gtf --fasta /path/to/genome.fa --ribopara S1.para.py
```

### GENCODE 注释流程 ✨ 新增

```bash
cd scripts/singularity_single_tool_tests

# 步骤 1: 运行 ORF 预测工具（以 Ribo-TISH 和 Ribotricer 为例）

# Ribo-TISH
./05_ribotish_quality.sh \
  --sample S1 \
  --bam S1.filtered.bam \
  --gtf annotation.gtf

./06_ribotish_predict.sh \
  --sample S1 \
  --bam S1.filtered.bam \
  --gtf annotation.gtf \
  --fasta genome.fa \
  --ribopara out_ribotish_quality/S1.para.py

# Ribotricer
./07_ribotricer_prepareorfs.sh \
  --gtf annotation.gtf \
  --fasta genome.fa \
  --prefix ribotricer_idx

./08_ribotricer_detectorfs.sh \
  --sample S1 \
  --bam S1.filtered.bam \
  --ribotricer-index out_ribotricer_prepareorfs/ribotricer_idx

# 步骤 2: 转换为 GENCODE 格式

./14_ribotish_to_gencode.sh \
  --sample S1 \
  --predict out_ribotish_predict/S1_pred.txt \
  --fasta genome.fa

./15_ribotricer_to_gencode.sh \
  --sample S1 \
  --tsv out_ribotricer_detectorfs/S1_translating_ORFs.tsv \
  --fasta genome.fa

# 步骤 3: 合并所有样本的 ORF（多样本时）

cat out_ribotish_to_gencode/S1.gencode.fa \
    out_ribotricer_to_gencode/S1.gencode.fa \
    > merged_orfs.fa

cat out_ribotish_to_gencode/S1.gencode.bed \
    out_ribotricer_to_gencode/S1.gencode.bed \
    > merged_orfs.bed

# 步骤 4: 运行 GENCODE ORF mapper（统一注释）

# 注意：需要先准备 Ensembl 注释目录
# 参考 gencode-riboseqORFs 文档：
# bash retrieve_ensembl_data.sh 110 GRCh38

./16_gencode_orf_mapper.sh \
  --project MyProject \
  --fasta merged_orfs.fa \
  --bed merged_orfs.bed \
  --ensembl-dir /path/to/Ens110
```

## GENCODE 注释流程详细说明

### 工作流程

```
ORF 预测工具输出 → 格式转换 → 合并文件 → GENCODE 映射 → 统一注释
   (Ribo-TISH,        (14-15)     (cat)      (16)       (.orfs.out)
    Ribotricer)
```

### 14_ribotish_to_gencode.sh

**用途**：将 Ribo-TISH 预测文件转换为 gencode-riboseqORFs 格式

**输入**：
- `--predict`: Ribo-TISH 的 `*_pred.txt` 文件
- `--fasta`: 基因组 FASTA（用于序列提取）
- `--sample`: 样本 ID（作为 study_id）

**输出**：
- `${SAMPLE}.gencode.fa`: GENCODE 格式 FASTA（`>ORF_NAME--STUDY_ID`）
- `${SAMPLE}.gencode.bed`: GENCODE 格式 BED（1-based 坐标）

**示例**：
```bash
./14_ribotish_to_gencode.sh \
  --sample S1 \
  --predict out_ribotish_predict/S1_pred.txt \
  --fasta genome.fa \
  --min-length 16 \
  --outdir out_ribotish_to_gencode
```

### 15_ribotricer_to_gencode.sh

**用途**：将 Ribotricer TSV 输出转换为 gencode-riboseqORFs 格式

**输入**：
- `--tsv`: Ribotricer 的 `*_translating_ORFs.tsv` 文件
- `--fasta`: 基因组 FASTA（用于序列提取）
- `--sample`: 样本 ID（作为 study_id）

**输出**：
- `${SAMPLE}.gencode.fa`: GENCODE 格式 FASTA
- `${SAMPLE}.gencode.bed`: GENCODE 格式 BED（1-based 坐标）

**质量过滤**：
- `--min-length`: 最小 ORF 长度（氨基酸，默认 16）
- `--min-phase-score`: 最小 phase score（默认 0.5）

**示例**：
```bash
./15_ribotricer_to_gencode.sh \
  --sample S1 \
  --tsv out_ribotricer_detectorfs/S1_translating_ORFs.tsv \
  --fasta genome.fa \
  --min-length 16 \
  --min-phase-score 0.5 \
  --outdir out_ribotricer_to_gencode
```

### 16_gencode_orf_mapper.sh

**用途**：统一多个工具/样本的 ORF 注释，去冗余，映射到 GENCODE

**输入**：
- `--fasta`: 合并的 GENCODE 格式 FASTA
- `--bed`: 合并的 GENCODE 格式 BED（1-based）
- `--ensembl-dir`: Ensembl 注释目录（需预先下载）
- `--project`: 项目 ID（输出文件前缀）

**Ensembl 目录要求**：
必须包含以下文件（通过 gencode-riboseqORFs 的 `retrieve_ensembl_data.sh` 生成）：
- `PROTEOME_FASTA`
- `TRANSCRIPTOME_FASTA`
- `SORTED_TRANSCRIPTOME_GTF`
- `TRANSCRIPT_SUPPORT`
- `PSITES_BED`

**输出**：
- `${PROJECT}.orfs.fa`: 统一的 ORF 序列
- `${PROJECT}.orfs.bed`: 统一的 ORF 坐标（1-based）
- `${PROJECT}.orfs.gtf`: GENCODE 格式 GTF 注释
- `${PROJECT}.orfs.out`: 详细 ORF 特征表（含分类、跨样本检测等）
- `${PROJECT}.altmapped`: 可选映射
- `${PROJECT}.unmapped`: 未映射 ORF

**去冗余参数**：
- `--collapse-threshold`: 相似度阈值（默认 0.9）
- `--collapse-method`: 去冗余方法（`longest_string` 或 `psite_overlap`）

**示例**：
```bash
./16_gencode_orf_mapper.sh \
  --project MyStudy \
  --fasta merged_orfs.fa \
  --bed merged_orfs.bed \
  --ensembl-dir /data/Ens110 \
  --min-length 16 \
  --collapse-threshold 0.9 \
  --collapse-method longest_string \
  --outdir out_gencode_orf_mapper
```

## 完整示例：从 BAM 到 GENCODE 注释

```bash
#!/bin/bash
# 完整的 GENCODE 注释流程示例

SAMPLE="S1"
BAM="S1_ribo.bam"
GENOME_FA="genome.fa"
GTF="annotation.gtf"
ENSEMBL_DIR="/data/Ens110"

# 0. 准备工作目录
cd scripts/singularity_single_tool_tests
export BIND_EXTRA="/mnt:/mnt"  # 如果需要

# 1. BAM 过滤
./01_sorf_bam_filter.sh \
  --sample $SAMPLE \
  --bam $BAM \
  --fai ${GENOME_FA}.fai \
  --outdir out_sorf_filter

# 2. Ribo-TISH
./05_ribotish_quality.sh \
  --sample $SAMPLE \
  --bam out_sorf_filter/${SAMPLE}.filtered.bam \
  --gtf $GTF

./06_ribotish_predict.sh \
  --sample $SAMPLE \
  --bam out_sorf_filter/${SAMPLE}.filtered.bam \
  --gtf $GTF \
  --fasta $GENOME_FA \
  --ribopara out_ribotish_quality/${SAMPLE}.para.py

# 3. Ribotricer
./07_ribotricer_prepareorfs.sh \
  --gtf $GTF \
  --fasta $GENOME_FA \
  --prefix genome_idx

./08_ribotricer_detectorfs.sh \
  --sample $SAMPLE \
  --bam out_sorf_filter/${SAMPLE}.filtered.bam \
  --ribotricer-index out_ribotricer_prepareorfs/genome_idx

# 4. 转换为 GENCODE 格式
./14_ribotish_to_gencode.sh \
  --sample $SAMPLE \
  --predict out_ribotish_predict/${SAMPLE}_pred.txt \
  --fasta $GENOME_FA

./15_ribotricer_to_gencode.sh \
  --sample $SAMPLE \
  --tsv out_ribotricer_detectorfs/${SAMPLE}_translating_ORFs.tsv \
  --fasta $GENOME_FA

# 5. 合并 ORF
cat out_ribotish_to_gencode/${SAMPLE}.gencode.fa \
    out_ribotricer_to_gencode/${SAMPLE}.gencode.fa \
    > merged_${SAMPLE}.fa

cat out_ribotish_to_gencode/${SAMPLE}.gencode.bed \
    out_ribotricer_to_gencode/${SAMPLE}.gencode.bed \
    > merged_${SAMPLE}.bed

# 6. GENCODE 映射
./16_gencode_orf_mapper.sh \
  --project MyProject \
  --fasta merged_${SAMPLE}.fa \
  --bed merged_${SAMPLE}.bed \
  --ensembl-dir $ENSEMBL_DIR

echo "✅ 完成！查看结果："
echo "   - ORF 分类表: out_gencode_orf_mapper/MyProject.orfs.out"
echo "   - GTF 注释: out_gencode_orf_mapper/MyProject.orfs.gtf"
```

## 注意事项

### GENCODE 注释特殊要求

1. **坐标系统**：
   - 输出的 BED 文件使用 **1-based** 坐标（gencode-riboseqORFs 要求）
   - 这与标准 BED 格式的 0-based 不同

2. **Ensembl 注释准备**：
   ```bash
   # 下载并准备 Ensembl 注释（在 gencode-riboseqORFs 仓库中）
   git clone https://github.com/jorruior/gencode-riboseqORFs.git
   cd gencode-riboseqORFs
   bash scripts/retrieve_ensembl_data.sh 110 GRCh38
   ```

3. **容器限制**：
   - `16_gencode_orf_mapper.sh` 需要容器中包含 `ORF_mapper_to_GENCODE_v1.1.py`
   - 或者需要挂载 gencode-riboseqORFs 仓库

4. **多样本合并**：
   - 可以合并来自不同样本、不同工具的 ORF
   - `study_id` 用于追踪 ORF 来源

## 输出解读

### 14-15 转换器输出

**FASTA 格式**：
```
>ENST00000456328_100000_111aa--S1
MAAGTLQSQLQNLQ*
```
- Header: `{GENE_ID}_{START}_{LENGTH}aa--{STUDY_ID}`
- 序列必须以 stop codon `*` 结尾

**BED 格式**（6 列，1-based）：
```
chr1  100000  100333  ENST00000456328_100000_111aa  S1  +
```

### 16 GENCODE mapper 输出

**主要输出文件 `*.orfs.out`**：
包含详细的 ORF 特征，例如：
- `orf_id`: 唯一 ORF 标识
- `orf_biotype`: ORF 类型（uORF, dORF, annotated, lncRNA_ORF 等）
- `gene_name`: 宿主基因名称
- `n_datasets`: 检测到该 ORF 的数据集/样本数量
- `X_S1`: 每个样本的检测标记（1/0）

**用途**：
- 识别高置信度 ORF（多样本/工具支持）
- ORF 功能分类
- 跨研究比较

## 故障排查

**问题 1**：`No ORFs passed the filters`

**原因**：ORF 太短或质量不达标

**解决**：
```bash
# 降低过滤阈值
--min-length 10 \
--min-phase-score 0.3
```

**问题 2**：`Missing required Ensembl file`

**原因**：Ensembl 目录不完整

**解决**：
```bash
# 重新下载 Ensembl 注释
cd gencode-riboseqORFs
bash scripts/retrieve_ensembl_data.sh 110 GRCh38
```

**问题 3**：容器中找不到 `ORF_mapper_to_GENCODE_v1.1.py`

**解决**：
```bash
# 临时方案：挂载 gencode-riboseqORFs 仓库
export BIND_EXTRA="/path/to/gencode-riboseqORFs:/opt/gencode-riboseqORFs"
```

## 相关资源

- [gencode-riboseqORFs GitHub](https://github.com/jorruior/gencode-riboseqORFs)
- [独立转换脚本文档](../../bin/README_CONVERTERS.md)
- [测试数据](../../test_data/)
- [GENCODE 集成总结](../../GENCODE_SESSION_UPDATE_20260117.md)
