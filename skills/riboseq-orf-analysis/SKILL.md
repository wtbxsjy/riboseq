---
name: riboseq-orf-analysis
description: >
  ORF 预测、unification（多工具合并）与 classification（分类）的分析流程。
  涵盖 pipeline 内 per-sample 预测参数（Ribo-TISH/Ribotricer/ORFquant/PRICE/RiboCode/rp-bp）、
  单工具独立运行（scripts/singularity_single_tool_tests/，含各工具前置依赖如
  ORFquant 需 RiboseQC _for_ORFquant）、unified_orfs 输出数据模型（metadata.tsv 列、
  工具名大小写约定、合并规则）、三个分类器（GENCODE orf_biotype / ORFquant ORF_category /
  ORF-type）、以及跳过 pipeline 的手动 bypass 路径（scripts/run_orf.py + singularity 容器执行）。
  当用户提到 ORF 预测/calling、unify/unification、分类/classify、uORF/dORF 分析、
  unified_orfs、orf_biotype、Ribo-TISH/Ribotricer/ORFquant 结果处理时使用本 skill。
---

# riboseq ORF 分析（calling → unification → classification）

## 决策树

- **Pipeline 已跑/要跑 ORF 预测+统一+分类** → 第 1-2 节（pipeline 内参数）
- **只想单独跑某个预测工具（不跑 pipeline）** → 第 3 节 → `references/single_tool.md`
- **Pipeline 跑完预测但 unify/分类失败或想单独重跑** → 第 4 节（手动 bypass）
- **要理解 unified_orfs 输出文件/列** → 第 5 节 → `references/unify_reference.md`
- **要理解分类器行为/输出** → 第 6 节 → `references/classifiers.md`
- **遇到 ORFquant 报错（BSgenome/library NULL）** → 第 7 节 → `references/orfquant_saga.md`

## 1. 全链路概览

```
5 个预测工具（per-sample，无 pooled 模式）
  Ribo-TISH → {sample}_pred.txt
  Ribotricer → {sample}_translating_ORFs.tsv
  ORFquant  → {sample}_Detected_ORFs.gtf.gz   (依赖 RiboseQC _for_ORFquant)
  PRICE     → {sample}.orfs.tsv
  RiboCode  → {sample}_collapsed.gtf.gz 等    (可选 --run_ribocode，需转录组 BAM)
  rp-bp     → bayes-factors.bed.gz            (可选 --run_rpbp)
      ↓ UNIFY_ORF_PREDICTIONS（scripts/unify_orf_predictions.py，容器 unify_orf.sif）
  unified_orfs.{bed,gtf,metadata.tsv,stats.txt} + expression_{summary,rpkm_tpm}.tsv
      ↓ 三个分类器并行 + ORF_QC + EXPRESSION_QUANT
  orf_classification/gencode   → gencode_results.orfs.out(.gz)  [orf_biotype]
  orf_classification/orfquant  → orfquant_classification.tsv   [默认 skip]
  orf_classification/orf_type  → orftype_classification.tsv     [orf_type]
  orf/unified_orfs_orf_confidence.tsv                            [OCS 置信度]
```

所有预测工具吃 sORF 过滤后的 BAM（`{sample}.sorf.filtered.bam`）；跳过 RiboseQC 会自动跳过 ORFquant。

## 2. Pipeline 内运行参数速查

```bash
# 关闭不需要的工具
--skip_rpbp true --skip_orfquant true --skip_orf_classify_orfquant true
# 启用 RiboCode / rp-bp
--run_ribocode --run_rpbp
# 关键调参（植物低深度）
--ribotricer_phase_score_cutoff 0.1
# unify 调参（wheat/soybean 用 24，rice 用 6）
--unify_orf_min_len 6 --unify_orf_frame_merge_min_overlap 0.9
# 分类模式（只跑 orf_type）
--orf_classify_mode orf_type
# 自定义容器
--unify_orf_container ... --gencode_orf_mapper_container ... --rpbp_container ...
```

ORF_QC 默认参数（一般不用动）：`orf_qc_min_consensus_tools=2`、
`orf_qc_confidence_weights='0.30,0.30,0.20,0.15,0.05'`、
`orf_qc_periodicity_min_f0=0.6`、`orf_qc_offset_max_delta=1`。

## 3. 单工具独立运行（不跑 pipeline）

全套单工具脚本：`scripts/singularity_single_tool_tests/01..20_*.sh`——每个都 mirror 对应
pipeline module（容器版本/参数一致），`00_run_order_example.sh` 是按依赖顺序串联的模板
（`DRY_RUN=0` 后真跑，需先填 BAM/FASTA/GTF/FAI/RRNA 路径）。

依赖链（⚠️ = 硬前置，缺了跑不了）：

```
01 sorf_bam_filter (BAM+.fai) ──→ {sample}.sorf.filtered.bam   ← 所有预测工具的共同输入
02 riboseqc_prepareannotation (GTF+FASTA) ──→ *_Rannot         （一次性）
03 riboseqc_analysis (filtered BAM+Rannot) ──→ *_for_ORFquant + *_P_sites_calcs
    ⚠️ 04 orfquant_run 必需 02 的 *_Rannot + 03 的 *_for_ORFquant（ORFquant 不吃 BAM）
05 ribotish_quality ──→ {sample}.para.py
    ⚠️ 06 ribotish_predict (+FASTA) 必需 05 的 .para.py
07 ribotricer_prepareorfs ──→ 候选 ORF index（一次性，与样本无关）
    08 ribotricer_detectorfs (+filtered BAM, --stranded) ──→ *_translating_ORFs.tsv
09 rpbp_prepare_genome（⚠️ 必需 rRNA FASTA）──→ *.orfs-genomic/exons.bed.gz
    10 rpbp_predict (+filtered BAM) ──→ bayes-factors.bed.gz
11 ribocode_detect（⚠️ 推荐 transcriptome BAM；低深度失败属正常）
12 = 03+04 一条龙；13 = orfquant_prepareannotation
PRICE 无单工具脚本（只有 pipeline module：-reads/-prefix/-genomic，OML 要绝对路径）
各工具输出 ──→ 17 unify_predictions ──→ 18 classify_orfs   ← 手动统一/分类（见下）
```

通用：BAM 已排序建索引；干净工作目录；镜像缓存到该目录 `containers/`；
WSL 路径 `export BIND_EXTRA="/mnt:/mnt"`。

**手动统一/分类**（单工具跑完后，不需要 pipeline）：
`17_unify_predictions.sh`（只支持 Ribo-TISH/Ribotricer/ORFquant 三种输入；
历史坑：旧版传的 `--merge-tolerance` flag 已被重构移除会报错，2026-08-18 已修）→
`18_classify_orfs.sh --mode gencode|orfquant|orf_type`（gencode 需 `--ensembl-dir`，
orfquant/orf_type 需 `--gtf`）。样本多或要 RiboCode/PRICE → 第 4 节 run_orf.py。

各工具详细命令/输入输出/容器解析 + 17/18 完整用法 → `references/single_tool.md`；
GENCODE 注释链（14/15/16 脚本）也在其中。

## 4. 手动 bypass（跳过 pipeline）

统一入口 `scripts/run_orf.py`（自动检测 orfont 加速包，没有则退回原脚本）：

```bash
# unify：合并各工具结果
python3 scripts/run_orf.py unify \
    --gtf ref.gtf --fasta ref.fa --output unified_orfs \
    --min-len 6 --threads 2 --frame-merge-min-overlap 0.9 \
    --ribotish <postfilter>/*_pred.txt \
    --ribotricer <postfilter>/*_translating_ORFs.tsv \
    --ribocode <ribocode>/*_collapsed.gtf.gz \
    --orfquant <orfquant>/*_Detected_ORFs.gtf.gz \
    --price <price>/*.orfs.tsv

# 分类（wrapper 实际 CLI 是 --input/--output_dir，不是 --bed/--metadata！）
python3 scripts/run_orf.py classify-gencode \
    --input unified_orfs --output_dir out_gencode \
    --ensembl_dir Ens58_oryza_sativa --cpus 16
python3 scripts/run_orf.py classify-orftype \
    --input unified_orfs --output_dir out_orftype --gtf ref.gtf
```

⚠️ **run_orf.py 的 docstring 已过时**（写的 `--bed u.bed --metadata u.tsv`），
实际 wrapper 参数是：`--mode gencode|orfquant|orf_type`、`--input <unify 输出前缀或完整路径>`、
`--output_dir`、`--gtf`（orfquant/orf_type 需要）、`--fasta`、`--ensembl_dir`（gencode 需要）、
`--gencode_impl original|fast|indexed_fast`、`--cpus`。`--input` 传 unify 的 prefix
（如 `unified_orfs`），wrapper 自动找 `{prefix}.metadata.tsv / .orfs.fa / .bed(.gz)`。

**容器内执行**（unify_orf.sif 缺 duckdb，需先 pip 装；用自带包装脚本，
自动处理 `--no-home --pid`、HOME/PYTHONUSERBASE、pip install duckdb；
`--dry-run` 可预览命令，`--no-duckdb` 跳过安装，`--log FILE` tee 日志）：

```bash
# run/ 在本机是仓库的兄弟目录（~/riboseq/run/），不是仓库内
bash skills/riboseq-orf-analysis/scripts/run_orf_in_container.sh \
  ~/riboseq/run/rice/containers/unify_orf.sif \
  -- unify --gtf ... --fasta ... --output unified_orfs \
      --ribotish <postfilter>/*_pred.txt ...
```

完整手动链路（singularity exec 细节、duckdb 流式分批、orfont 加速实测）→
`references/manual_bypass.md`。

## 5. unify 输出数据模型

输出 8 类文件：`unified_orfs.bed.gz`（BED12）、`.gtf.gz`、`.metadata.tsv`、
`.stats.txt`、`_expression_summary.tsv`、`_expression_rpkm_tpm.tsv`、versions.yml。

- metadata.tsv 26 列（orf_id/chrom/strand/start/end/length_aa/exon_blocks/gene_id/
  transcript_id/tools/samples/tool_scores/tool_pvalues/unique_psites/pN/sequence/
  aa_sequence/is_cds_overlap/...），`orf_id = ORF_{i}_{gid}`
- **工具名大写约定**（解析 tools/sources 列必须用）：`'Ribo-TISH'`、`'Ribotricer'`、
  `'ORFquant'`、`'PRICE'`、`'RiboCode'`
- 合并三步：exact-match → frame-aware（同框合并）→ overlap grouping（选代表）
- 实测规模：rice 23 样本 369,938 ORFs / 35min / 65.7GB；human 62 样本 956,971 /
  18.3% frame-merge reduction

完整列清单与 stats.txt 解读 → `references/unify_reference.md`。

## 6. 三个分类器

| 分类器 | 输出 | 要点 |
|---|---|---|
| GENCODE（默认） | `gencode_results.orfs.out(.gz)`：orf_biotype ∈ CDS/dORF/uORF/doORF/uoORF/intORF/lncRNA | 必需 `--ensembl_dir`（5 个标准文件）；容器需 bedtools+BioPython；wrapper 自动做 BED12→BED6 和蛋白 FASTA（key `{orf_id}--{sample_id}`） |
| ORFquant | `orfquant_classification.tsv`：ORF_category_Gen/Tx/Tx_compatible | **当前默认 skip**（`--skip_orf_classify_orfquant` 不在 schema，传了会 WARN）；数据表重写中 |
| ORF-type | `orftype_classification.tsv` | 最快（几分钟）；行数 ≥ unified ORF 数（含 frame 变体）；canonical_CDS/uORF/dORF/overlap_uORF 等 |

细节（实现选择、坑）→ `references/classifiers.md`。

## 7. 坑速查

1. **RiboCode `.gtf.gz` 解析失败**：老代码 `endswith('.gtf')` 匹配不到 `.gtf.gz` → 62 个文件全部 "missing required columns"（已修，`_open()` fallback）。
2. **ORFquant `library(NULL)` / `GTF_annotation not found`**：非模式生物无 BSgenome → 5 轮 monkey-patch 迭代史与最终 FaFile 方案 → `references/orfquant_saga.md`。
3. **channel 单消费者死锁**：`.into{}` 在 resume 的缓存任务上不可用；ORF_QC 最终改为磁盘 glob + `optional: true` 占位文件（已修进代码，遇到时按 riboseq-pipeline-run 的 troubleshooting §4 处理）。
4. **PRICE/GEDI 参数**：v1.0.5 只认 `-reads/-prefix/-genomic`，`-genomic` 要 OML 绝对路径。
5. **GENCODE 分类结果 lncRNA >90% / 只剩 2 类 biotype**：蛋白 FASTA header 版本后缀与 GTF protein_id 不匹配（`ENSP00000493376.2` vs `ENSP00000493376`）→ 2026-08-19 已修（commit 7281ecd，`_dual_key_index` 双 key 索引；FASTA 重复 key 崩溃同源）。无点号版本的物种（rice/maize）不受影响，重跑结果不变。诊断与细节 → `references/classifiers.md` §1。
6. **unify 输入路径与后缀**：各工具按后缀推断 sample id（`infer_sample_id_from_prediction_path`），文件名不符合 `{sample}_xxx` 约定时 sample 归属错乱。
7. **contig 命名不统一（PRICE 缺 `chr` 前缀）——整条 PRICE 支线被静默清零**（2026-09-10 实测 PRJEB26593）：PRICE（GEDI）的 `.orfs.tsv` / `_Detected_ORFs.gtf` 用 Ensembl 风格（`1`…`X`,`Y`,`MT`），其余 4 工具与参考 FASTA/BAM 都是 `chr` 前缀。后果：PRICE ORF 取序列得到全 `N`（`extract_sequence()` 直接 `fasta[cand.chrom]`）、featureCounts/psite 定量恒 0、与同名坐标的其他工具 ORF 永远合不上（`id_key` 含 chrom）——327,119/1,215,479 = **26.9% 的 ORF 被清零**（metadata `total_psites>0` 占 0.00% vs 其他工具 99.98%；补 `chr` 前缀后单样本 55.1% 有 ≥1 P-site）。**已修**（未提交）：`_make_chrom_normalizer(gtf_index)` 在解析完所有文件后、合并/查表/取序列前按参考命名归一化（幂等、基于 `_chrom_aliases()`，GL/KI 脚手架不动），并重建 `cand.id_key`。⚠️ 修复后必须重跑 UNIFY，且**所有 ORF ID 会整体重排**，post_analysis 集合需按坐标重映射；新增预测工具时先核对 contig 约定。

## 完成后

- ORF 定量（RPKM/TPM）、P-site 过滤、ggRibo 图 → **riboseq-quant-te** skill
- 结果体检 → riboseq-pipeline-run 的 `check_pipeline_outputs.sh`
