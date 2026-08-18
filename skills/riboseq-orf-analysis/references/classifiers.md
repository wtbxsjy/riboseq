# 三个分类器：输入输出与已知问题

模块：`modules/local/classify_orfs/main.nf`，三个 process 并行：
`CLASSIFY_ORFS_GENCODE` / `CLASSIFY_ORFS_ORFQUANT` / `CLASSIFY_ORFS_ORF_TYPE`。
输出目录：`result/orf_classification/{gencode,orfquant,orf_type}/`。
实测速度：ORF_TYPE 最快（几分钟）、GENCODE 次之（17 万 ORF 约 1 分钟，
配置 cpus 16 / 64GB / 48h 时实测 5m8s、3.6GB）。

## 1. GENCODE 分类器（默认启用）

**脚本链**：`classify_orfs_wrapper.py --mode gencode` →
`gencode-riboseqORFs/ORF_mapper_to_GENCODE_v1.1.py`（或 fast/indexed_fast 实现）。

**硬要求**：
- `--orf_classify_ensembl_dir`（缺了直接报错），目录含 5 个标准文件
  （TRANSCRIPTOME_FASTA / SORTED_TRANSCRIPTOME_GTF / PROTEOME_FASTA /
  TRANSCRIPT_SUPPORT / PSITES_BED），准备方法见 riboseq-data-prep skill
- 容器 `gencode_orf_mapper.sif` 必须含 bedtools + BioPython

**wrapper 自动转换**（CLI 传 `--input <prefix>`）：
- BED12 → BED6（col[4] 填 study_id = metadata `samples` 列的第一个样本）
- 蛋白 FASTA key：`{orf_id}--{study_id}`（mapper 硬性格式）

**输出**（rice 实测 `.out.gz` 压缩；human 会话为未压缩 `.out` 41MB）：

```
gencode_results.orfs.out(.gz)    # 主结果
gencode_results.orfs.{bed,gtf.gz,fa,allframes.bed,frames.bed}
gencode_results.logs
```

`.out` 列：`orf_id version chrm starts ends strand trans gene gene_name
orf_biotype gene_biotype pep orf_length ...`（第 7 列 = orf_biotype）。

**orf_biotype 取值**：CDS / dORF / uORF / doORF / uoORF / intORF / lncRNA。

**验证**：`zcat gencode_results.orfs.out.gz | cut -f7 | sort | uniq -c`；
lncRNA >90% = 蛋白 header 不匹配坑（见 riboseq-data-prep）。
未与任何注释转录本重叠的 ORF 归 "intergenic"（下游分析时与非 CDS 一起处理，
回收约 53% 否则被静默丢弃的 ORFs）。

**实现选择**：`--gencode_impl original|fast|indexed_fast`（mouse 用 indexed_fast；
982K ORFs 时 mapper 崩溃 → 分批 + `csv.field_size_limit(sys.maxsize)`）。

## 2. ORFquant 分类器（当前默认 skip）

**状态**：`nextflow.config` 里 `skip_orf_classify_orfquant = true`
（"Temporarily disabled; pending data.table rewrite"）。
参数不在 nextflow_schema.json → 显式传参会有 `WARN: invalid input values`（无害）。

**设计**（`scripts/class_orf/orfquant_orf_classify.R`）：
- 输出 `orfquant_classification.tsv`：`ORF_category_Gen`（基因组级）、
  `ORF_category_Tx`（转录本级）、`ORF_category_Tx_compatible`（最佳 isoform 级）
- `project_to_tx_coords()` 把 ORF 基因组块沿转录本外显子链投影到转录本坐标
  （修复了旧基因组近似对 ~18% ORF 的错误）
- `normalize_annotation()` 读 exon + CDS feature，返回 exon_txs 与
  cds_txs_tx_coords

重新启用前确认 data.table 重写完成。

## 3. ORF-type 分类器

**脚本**：`class_orf/class_ORFtype.py`（wrapper `--mode orf_type`，需要 `--gtf`）。

**输出**：

```
orftype_classification.tsv            # 主结果
orftype_results.{orfs.bed,orfs.fa,orfs.gtf,orfs.out,orfs.pep.fa,logs}
```

tsv 列：`orf_id orf_name orf_type chrom strand start end source gene_id gene_name
transcript_id orf_length score source_count samples tools sequence start_codon
aa_sequence cds_overlap exon_blocks orf_type_category`。

- **行数 ≥ unified ORF 数**（含 frame 变体；human 1,171,859 行 vs 956,971 ORFs）
- orf_type 取值：canonical_CDS / uORF / dORF / overlap_uORF 等基因级分类
- 不依赖 Ensembl 目录，速度最快

## 与 ORF_QC 的关系

ORF_QC 用分类信息计算置信度（OCS 五维加权：translation/agreement/coverage/
periodicity/readlevel）；历史教训：ORF_QC 与 CLASSIFY 的 channel 单消费者死锁
（已修：ORF_QC 改为磁盘 glob 读取 unified 文件 + optional 输出占位）。
