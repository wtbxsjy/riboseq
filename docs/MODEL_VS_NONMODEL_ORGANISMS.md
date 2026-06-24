# 模式生物 vs 非模式生物处理差异分析

> 适用版本：`dev` 分支，2026-06-24  
> 分析范围：人/鼠（模式生物）vs 水稻/玉米（非模式生物）

## 一、BSgenome / 基因组包 — 最核心的架构差异

这是技术层面最根本的区别。

### 背景

| 生物类型 | BSgenome 包来源 | ORFquant 原生行为 |
|----------|----------------|-------------------|
| 模式生物（人/鼠） | Bioconductor 预编译包（如 `BSgenome.Hsapiens.UCSC.hg38`） | `library(genome_package)` 直接加载 |
| 非模式生物（水稻/玉米） | 不存在预编译包 | `genome_package = NULL` → `library(NULL)` 崩溃 |

### 流水线的三层防护

| 模块 | 策略 | 代码位置 |
|------|------|----------|
| **RiboseQC** | 始终 `forge_BSgenome=TRUE`，从用户 FASTA 临时构建 BSgenome。模式生物虽冗余但无害；非模式生物必须如此 | `modules/local/riboseqc/prepareannotation/main.nf:41` |
| **ORFquant** | 运行时 monkey-patch `load_annotation()`，增加 NULL 守卫 | `modules/local/orfquant/main.nf:102-121` |
| **PREPARE_FOR_ORFQUANT** | 同样的 patch，确保 P-site 校正后重新量化也能处理 NULL | `modules/local/prepare_for_orfquant/main.nf:58-79` |

### ORFquant NULL 守卫补丁

```r
# modules/local/orfquant/main.nf:102-121
unlockBinding("load_annotation", asNamespace("ORFquant"))
patched_load_annotation <- function(path) {
    GTF_annotation <- get(load(path))
    if (is(GTF_annotation$genome, "FaFile")) {
        genome_sequence <- GTF_annotation$genome
    } else {
        genome_pkg <- GTF_annotation$genome_package
        if (!is.null(genome_pkg) && nchar(genome_pkg) > 0) {
            library(genome_pkg, character.only = TRUE)
            genome_sequence <- get(genome_pkg)
        } else {
            genome_sequence <- NULL        # ← 非模式生物走这条路径
        }
    }
    GTF_annotation <<- GTF_annotation
    genome_seq <<- genome_sequence
}
assign("load_annotation", patched_load_annotation, envir = asNamespace("ORFquant"))
lockBinding("load_annotation", asNamespace("ORFquant"))
```

> ⚠️ **关键影响**：当 `genome_sequence <- NULL` 时，ORFquant 无法提取基因组序列来验证 ORF 编码潜能，会导致绝大多数基因因缺少序列证据而无法预测 ORF。因此非模式生物必须使用 `forge_BSgenome=TRUE` 的 RiboseQC 输出（即 `for_ORFquant` 文件中的 `genome` 字段为 `FaFile` 对象，走 `is(GTF_annotation$genome, "FaFile")` 分支）。

> 相关 Bug（CLAUDE.md 第 16 条）：`--additional_fasta` 自动生成的 GTF 只有 `exon` 特征，导致 RiboseQC TxDb 构建失败 → `for_ORFquant` 不可用 → ORFquant 无结果。

---

## 二、GENCODE vs Ensembl 注释体系

通过 `--gencode` 参数控制（`nextflow.config:25`，默认 `false`）。

| 维度 | `--gencode true`（人/鼠 GENCODE） | `--gencode false`（水稻/玉米 Ensembl Plants） |
|------|-----------------------------------|-----------------------------------------------|
| GTF 基因类型属性 | `gene_type` | `gene_biotype` |
| Salmon 索引 | 传 `--gencode` 标志，适配管道符分隔的 FASTA header | 不传额外标志 |
| 转录本 FASTA 预处理 | `cut -d "\|" -f1` 截掉管道符后的字段 | 不需要预处理 |
| 典型参考来源 | GENCODE（人 v47, 鼠 vM36） | Ensembl Plants（水稻 IRGSP-1.0, 玉米 Zm-B73-REFERENCE-NAM-5.0） |

### 关键代码

**`subworkflows/local/prepare_genome/main.nf:116`**
```groovy
def biotype = gencode ? "gene_type" : "gene_biotype"
```

**`subworkflows/local/prepare_genome/main.nf:145-149`**
```groovy
if (gencode) {
    PREPROCESS_TRANSCRIPTS_FASTA_GENCODE(ch_transcript_fasta)
    ch_transcript_fasta = PREPROCESS_TRANSCRIPTS_FASTA_GENCODE.out.fasta
}
```

**`conf/modules.config:61`**
```groovy
params.gencode ? '--gencode' : ''
```

---

## 三、lncRNA 生物型分类

流水线中**唯一硬编码了物种名**的地方。两个独立的分类器各维护一份 `LNCRNA_BIOTYPES` 集合，必须同步更新。

| 分类器 | 文件 | 行号 |
|--------|------|------|
| ORFquant 分类 (R) | `scripts/class_orf/orfquant_orf_classify.R` | 221-230 |
| ORF-type 分类 (Python) | `scripts/class_orf/class_ORFtype.py` | 26-35 |

### 完整的 biotype 列表

**模式生物（人/鼠）专有**：
```
lncRNA, lincRNA, antisense, processed_transcript,
sense_intronic, sense_overlapping, non_coding,
3prime_overlapping_ncrna, bidirectional_promoter_lncrna
```

**非模式生物额外 biotype**（代码中明确标注了物种）：
```
ncRNA            # Oryza sativa (水稻): 通用长链非编码 RNA
antisense_RNA    # Oryza sativa: 反义长链非编码 RNA
misc_non_coding  # Zea mays (玉米): 非编码基因汇总类别
```

### 分类逻辑

```r
# orfquant_orf_classify.R:234
lnc_mask <- gene_feat$gene_biotype %in% LNCRNA_BIOTYPES
lncrna_genes <- unique(gene_feat$gene_id[lnc_mask])
```

如果 ORF 的 `gene_id` 属于 lncRNA 基因，则分类为 `"lncRNA"`，不再进行位置性分类。

> ⚠️ **关键依赖**：如果某物种的 lncRNA biotype 不在集合中，该物种的 lncRNA 基因上的 ORF 将被错误分类为位置性类别（如 `"novel_Upstream"`），造成系统性偏差。添加新非模式物种时，必须检查其在 Ensembl Plants 中的 `gene_biotype` 取值。

---

## 四、Contig 过滤正则

`nextflow.config:110`

```
sorf_exclude_contigs_regex =
  ^(chr)?(M|MT|Mt|chrM|chrMT|chrMt|ChrM|ChrMT|ChrMt)$  ← 线粒体（人/鼠/植物）
  |^(chr)?(C|CP|Pt|chrC|chrCP|chrPt|ChrC|ChrCP|ChrPt)$ ← 叶绿体/质体（仅植物）
  |^chrUn_.*|.*_random$|.*_alt$|.*_fix$                  ← 未定位 scaffold（人/鼠）
```

| 生物类型 | 需要排除的 contig | 正则覆盖情况 |
|----------|-------------------|-------------|
| 人/鼠 | `chrM`, `chrUn_*`, `*_random`, `*_alt`, `*_fix` | ✅ 默认即可 |
| 水稻 (IRGSP-1.0) | `Mt`, `Pt` (无 `chr` 前缀) | ✅ 默认即可 |
| 玉米 (B73) | `Mt`, `Pt`, `B73V4_ctg*` | ⚠️ 需自行添加 contig 模式 |

---

## 五、性能参数与数据特征

### Ribotricer 周期性检测灵敏度

`nextflow.config:101`
```
ribotricer_phase_score_cutoff = null  // 默认 3/7 ≈ 0.429
```

| 数据类型 | 推荐 cutoff | 原因 |
|----------|------------|------|
| 哺乳动物高深度 Ribo-seq | null (0.429) | 3-nt 周期性信号强 |
| 植物/低深度样本 | **0.1** | metagene 周期性较弱 |

代码注释：*"Lower this (e.g. 0.1) for plant/low-depth samples where metagene periodicity is weak."*

### 基因组大小

| 物种 | 参考基因组 | 大小 | STAR 索引内存 |
|------|-----------|------|--------------|
| 人 | GRCh38 | ~3.1 GB | 30 GB+ |
| 鼠 | GRCm39 | ~2.7 GB | ~28 GB |
| 水稻 | IRGSP-1.0 | ~400 MB | ~8 GB |
| 玉米 | B73-NAM-5.0 | ~2.3 GB | ~20 GB |

### 分类资源

`conf/modules.config:959-969` 中 `CLASSIFY_ORFS_ORFQUANT` 和 `CLASSIFY_ORFS_GENCODE` 按 "Mouse-scale" 分配 16-32GB / 24h。对于基因组较小的非模式生物可能过度分配，但不影响正确性。

---

## 六、设计哲学

流水线**没有** `if model_organism then X else Y` 的显式分支，采用**参数驱动 + 通用补丁**策略：

```
┌─────────────────────────────────────────────────────────┐
│                 参数驱动（用户配置）                       │
│  --gencode        → 切换 GENCODE / Ensembl 行为          │
│  --sorf_exclude*  → 自定义 contig 过滤                   │
│  --ribotricer_phase_score_cutoff → 调节周期性灵敏度      │
├─────────────────────────────────────────────────────────┤
│                 通用补丁（代码内建）                       │
│  forge_BSgenome=TRUE  → RiboseQC 始终可用               │
│  NULL genome_pkg guard → ORFquant 适应所有物种            │
│  lncRNA biotypes 全集  → 分类器覆盖人/鼠/稻/玉米          │
│  版本号剥离            → riboWaltz 不受 BAM/GTF 差异影响  │
└─────────────────────────────────────────────────────────┘
```

## 七、新物种适配检查清单

添加新的非模式生物时，按以下清单逐项检查：

1. **lncRNA biotype** — 查看该物种在 Ensembl/Ensembl Plants 中的 `gene_biotype` 取值，若有不在现有 `LNCRNA_BIOTYPES` 中的非编码类别，同步更新：
   - `scripts/class_orf/orfquant_orf_classify.R` 中的 `LNCRNA_BIOTYPES`
   - `scripts/class_orf/class_ORFtype.py` 中的 `LNCRNA_BIOTYPES`

2. **Contig 命名** — 检查线粒体/叶绿体 contig 名称是否符合 `sorf_exclude_contigs_regex` 的默认值，不符合则传自定义正则

3. **GENCODE 标志** — 确认参考 GTF 来源（GENCODE vs Ensembl），正确设置 `--gencode`

4. **Ribotricer 阈值** — 若数据深度低或周期性弱，降低 `--ribotricer_phase_score_cutoff`（如 `0.1`）

5. **BSgenome 不可用** — 确保 RiboseQC 正常运行（`forge_BSgenome=TRUE`），产生可用的 `for_ORFquant` 文件

6. **参考染色体命名** — 检查 BAM 中的染色体名与 GTF/FASTA 是否一致，必要时在参考准备阶段处理
