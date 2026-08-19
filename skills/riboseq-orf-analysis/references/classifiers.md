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

`.out` 列（**有 header 行**）：`orf_id version chrm starts ends strand trans gene
gene_name orf_biotype gene_biotype pep orf_length ...`（**orf_biotype = 第 10 列**，
第 7 列是 trans；统计时先 `tail -n +2`）。

**orf_biotype 取值**：CDS / dORF / uORF / doORF / uoORF / intORF / lncRNA。

**验证**：`zcat gencode_results.orfs.out.gz | tail -n +2 | cut -f10 | sort | uniq -c`。
健康参考（本机 rice / maize，7 类齐全）：rice 151,848 ORF → CDS 107K / dORF 19K /
doORF 13K / intORF 7.5K / uORF 3.2K / uoORF 1.4K / lncRNA 562；maize 232,594 ORF →
CDS 195K / doORF 11K / intORF 8.5K / uORF 7K / uoORF 4.1K / dORF 3.9K / lncRNA 2.1K。
未与任何注释转录本重叠的 ORF 归 "intergenic"（下游分析时与非 CDS 一起处理，
回收约 53% 否则被静默丢弃的 ORFs）。

### ⚠️ 2 类塌缩坑（2026-08-19 修复，commit 7281ecd）

**症状**：biotype 分布只剩 lncRNA/CDS 两类（或 lncRNA >90%）；或分类直接崩溃
报 duplicate key / ValueError。Arabidopsis 与 Lishuqi 生产均实测塌缩。

**根因**（`functions.py::load_fasta()` 原版直接 `SeqIO.index`，两种失败模式）：
1. **版本后缀不匹配**：GTF `protein_id` 无版本（`ENSP00000493376`）而 FASTA
   header 带版本（`ENSP00000493376.2`）→ 查找 KeyError → 绝大多数 ORF 映射失败
   → 塌缩成 2 类。点号版本的物种（人/拟南芥 Gencode/Ensembl）必踩。
2. **FASTA 重复 key**：`retrieve_ensembl_data.sh` 旧版 `strip_fasta_versions`
   用 `FS="."` 取 $1 切掉版本 → `AT1G01020.1` 塌成基因级 `AT1G01020`，
   isoform 撞 key → SeqIO.index 直接抛 ValueError。

**修复**：
- `_dual_key_index()`：每个记录同时注册精确 key 与剥版本 key
  （`_strip_version()` = 去掉 `\.\d+$`，first-wins）；SeqIO.index 抛 ValueError
  时回退为逐条 parse + `setdefault` 去重（不再崩）。
- `retrieve_ensembl_data.sh`：`strip_fasta_versions` 改为保留第一个空白分隔
  token（`>AT1G01020.1 cdna...` 完整保留版本后缀）。

**影响面与重跑安全性**：
- 无点号版本的物种（rice `Os01t0100100-01`、maize `Zm00001eb000010_P001`）不受
  影响：剥版本是 no-op，双 key 索引退化为精确索引，**重跑结果不变**（本机
  2026-08 验证：maize proteome 72,539 条构建 1.3s、0 条剥版本额外键）。
- 生效路径：pipeline 的 CLASSIFY_ORFS_GENCODE stage 的是仓库 `scripts/` 下的
  mapper + functions.py（wrapper 从自身目录解析脚本），16 单工具脚本也 `cp`
  仓库副本 → 两条路径修复都生效；容器内 `/opt/gencode-riboseqORFs` 的 clone
  不会被用到。
- **仅 original 实现**：`--gencode_impl fast|indexed_fast`（`scripts/class_orf/`
  下的 run_gencode_classify_fast/indexed.py）有各自加载器（load_fasta_by_orf_id /
  load_fasta_records），未经过此修复，用这两种实现时另行验证。
- 内存注意：`_dual_key_index` 把全部记录物化进内存（比 SeqIO.index 高，maize
  转录组+蛋白组约几百 MB 级），process_medium 标签够用；超大转录组（人 GENCODE
  ~250K 条）也实测无碍。

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
