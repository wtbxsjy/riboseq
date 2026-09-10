# GSE120762（Jackson 2018 Nature）小鼠粘膜免疫重分析计划

> 创建：2026-09-01。论文：*The Translation of Non-Canonical Open Reading Frames Controls Mucosal Immunity*（Jackson et al., Nature 2018）。
> 数据目录：`~/riboseq/data/mouse_GSE120762/`（16 runs，28 fastq.gz，HiSeq 2000）。
> 项目目录（待建）：`~/riboseq/run/mouse_GSE120762/`。

## 1. 数据集构成（ENA 元数据 + 论文 Methods）

| 组 | runs | 布局 | 说明 | 论文用途 |
|---|---|---|---|---|
| RibosomeProfiling | SRR7956050-53 | SE 75nt | TruSeq Ribo Profile (mammalian) kit，CHX 50μg/ml 2min，Rep1/2 × CHX_NT/CHX_LPS | ORF 发现、TE、RRS |
| RNA-seq | SRR7956038-43 | PE 75bp | poly(A)+，BMDM WT×3 / LPS×3 | TE 分母、差异表达 |
| RiboTag | SRR7956044-49 | PE 75bp | HA-RPL22 IP（核糖体结合 mRNA），NT/6hr/24hr ×2 rep | lncRNA "expressed" 判定（FPKM≥1） |

模型：BMDM（骨髓来源巨噬细胞），LPS 刺激。核心生物学问题：lncRNA 中非经典 ORF 的翻译（代表性基因 Aw112010）。

## 2. 论文分析要点（Methods 摘要）

- 原始流程：TopHat2 → GRCm38；过滤 AS 0..-2、比对数 ≤2；Cufflinks + DESeq2。
- lncRNA ORF 发现：ORFfinder / RiboCode 1.2.10，起始密码子 **ATG/CTG/TTG/GTG**，ORF > 30nt；RibORF PME>0.6；BLASTX→CHESS2 + PhyloCSF 评估编码潜力。
- TE = 核糖体读段覆盖外显子数 / RNA 读段覆盖外显子数（长度归一化）；RRS = TE ÷（3'UTR 的 ribo/RNA 比）；阈值 **RRS≥7、TE≥0.0001** 视为蛋白编码基因组水平。
- RiboTag：lncRNA "expressed" = 三种 RiboTag 处理中至少一种 FPKM≥1；RNA-seq "detectable" = FPKM≥0.1。

## 3. 重分析策略

### 阶段 1：主流程（10 样本，一个 pipeline run）

samplesheet：4 RP（type=riboseq）+ 6 RNA-seq（type=rnaseq）。

- **分组**：RP 的 CHX_NT 与 RNA-seq 的 WT 均为"未受 LPS 刺激"→ 统一 group=`NT`；LPS 刺激 → group=`LPS`。
  - NT: SRR7956050, 7956052, 7956038, 7956039, 7956040（ribo 2 + rna 3）
  - LPS: SRR7956051, 7956053, 7956041, 7956042, 7956043（ribo 2 + rna 3）
- **deltaTE**：`scripts/contrasts.csv` → `lps_vs_nt,group,NT,LPS`。每 seq_type×group 格 ≥2 样本（2/2 与 3/3），满足交互模型守卫。
- **ORF 发现**：Ribo-TISH + Ribotricer + ORFquant + RiboCode（当前默认开启；RiboCode 默认起始密码子含 ATG/CTG/TTG/GTG，与论文一致）→ unify（min_len 6aa，frame-merge overlap 0.9）→ 三个分类器（GENCODE / ORFquant / ORF-type）。lncRNA 编码 ORF 候选 = `orf_biotype=lncRNA`。
- **TE 定量**：QUANTIFY_ORFS（featureCounts，unified ORF BED）→ MERGE_COUNTS → DESEQ2 deltaTE（lps_vs_nt）。

### 阶段 2：RiboTag 补充分析（流程外，run 1 完成后）

RiboTag 不进阶段 1 samplesheet 的原因：deseq2_deltate.R 将所有 type=rnaseq 行并入 "rna" 臂 —— RiboTag（核糖体结合 mRNA IP）与 poly(A) RNA-seq 是不同 assay，混入会污染 deltaTE 的 rna 臂解释。

完成后手工执行（复用 run 1 的 STAR index 与容器）：
1. STAR 比对 6 个 PE RiboTag 样本（同 run 1 参数）。
2. featureCounts：对 unified_orfs.bed + GENCODE 基因注释分别计数。
3. 按论文阈值：lncRNA "expressed"（任一 RiboTag 处理 FPKM≥1）+ RNA-seq "detectable"（FPKM≥0.1）过滤；计算 TE（RRS）等指标。

### 阶段 3：下游解读

- 对 unified ORFs 施加论文阈值（RRS≥7 且 TE≥0.0001）筛选翻译活跃的非经典 ORF。
- 重点核对 Aw112010（论文主角 lncRNA）在 GRCm39 上的 ORF 预测与分类结果。

## 4. 参考与容器配置

- **参考（GRCm39，与人类项目"现代 GENCODE"原则一致；论文用 GRCm38，重分析统一 GRCm39）**：
  - 复用 `~/riboseq/run/mouse_Mucosal_Immunity/reference/` 下已解压文件（源头 MyDrive）：`mouse.genome.fa`（GENCODE M35 / Ensembl 112）、`mouse.gtf`、`mouse.transcripts.fa`、`mouse_final_contamination.fasta`。
- **分类目录**：`~/riboseq/run/mouse_tissue/reference/ensembl/Ens110_mus_musculus`（5 symlink 齐全：TRANSCRIPTOME_FASTA / SORTED_TRANSCRIPTOME_GTF / PROTEOME_FASTA / TRANSCRIPT_SUPPORT / PSITES_BED）。
  - ⚠️ 已知版本差：run GTF 为 Ensembl 112（GENCODE M35），分类目录为 Ensembl 110。两者同 GRCm39、同 ENSMUST ID 空间，旧 mouse_Mucosal_Immunity 项目即用此组合。若分类遗漏明显，后续可从 MyDrive 重建 M35 匹配目录。
- **容器**：从 `~/riboseq/run/human_GSE208041/containers/` 复制全套（orfquant_patched.sif、unify_orf.sif、gencode_orf_mapper.sif、ribowaltz.sif、price.sif、deseq2_deltate.sif）。prepare_workflow.py 无 --deseq2-container 参数，deseq2 容器手动复制并在 run_pipeline.sh 中加 `--deseq2_container`。
- **strandedness**：全部 `auto`（ENA 无链信息；与 PRJEB26593/GSE208041 惯例一致）。PE 支持已在 main.nf 确认（fastq_2 分支存在）。

## 5. run_pipeline.sh 参数（镜像 GSE208041 + PRJEB26593 对比部分）

```
--aligner star  --max_memory 150.GB  --max_cpus 16  --max_time 48.h
--save_reference true
--fasta/--gtf/--transcript_fasta/--contaminant_fasta   (mouse 参考)
--orf_classify_ensembl_dir ~/riboseq/run/mouse_tissue/reference/ensembl/Ens110_mus_musculus
--deseq2_container .../containers/deseq2_deltate.sif
--skip_rpbp true
--orfquant_psite_correction true
--ribotish_fail_on_empty/--ribotricer_fail_on_empty/--orfquant_fail_on_empty/
--rpbp_fail_on_empty/--price_fail_on_empty/--ribocode_fail_on_empty true
--unify_orf_min_len 6  --unify_orf_frame_merge_min_overlap 0.9
--contrasts scripts/contrasts.csv
```

- sORF 读长过滤：**保持默认 28-30nt**（TruSeq Ribo Profile CHX 常规足迹；旧 MM 项目设 0/0 属其特殊处理，不沿用）。若 RiboseQC 显示 31nt 主导再调 `--sorf_read_len_max 32` + ribowaltz 读长。
- `--gencode_classify_impl`：默认 original（与两个人类项目一致）。
- 无 pathogen 参数；有 contrasts + rnaseq → deltaTE 完整模式。
- 资源：机器 56 核/188G 已跑 PRJEB26593(42) + GSE208041(32)，本项目 --max_cpus 16。

## 6. 执行步骤（待下载完成）

1. `prepare_workflow.py -w ~/riboseq/run/mouse_GSE120762/ -d ~/riboseq/data/mouse_GSE120762 --species mouse -r ~/riboseq/run/mouse_Mucosal_Immunity/reference --contaminant-dir ~/MyDrive/sORF_Discovery_Project/contamination_indices --orf-classify-ensembl-dir ~/riboseq/run/mouse_tissue/reference/ensembl/Ens110_mus_musculus --orfquant-container ... --unify-orf-container ... --gencode-orf-mapper-container ... --ribowaltz-container ... --price-container ... --fail-on-empty --unify-orf-min-len 6 --unify-orf-frame-merge-min-overlap 0.9 --orfquant-psite-correction --max-cpus 16 --max-memory 150.GB --max-time 48.h`
   - 注意：prepare 需在下载完成后跑（它会在当时建立 data/ 符号链接）。
   - prepare 生成的 samplesheet 是单一 type → 手工覆盖为 10 样本混合 type 版本。
2. 手工写 `scripts/contrasts.csv`（`lps_vs_nt,group,NT,LPS`）。
3. 复制 deseq2_deltate.sif 进项目 containers/；编辑 run_pipeline.sh（加 deseq2_container、contrasts、6 个 fail_on_empty、unify 参数；核对与第 5 节一致）。
4. 清理 prepare 可能引入的跨物种参考符号链接污染（历史问题）。
5. `bash run/human_GSE208041/process/run_pipeline.sh` 同款方式启动 + Monitor 盯关键节点（STAR 完成、RiboseQC、ORF 预测、UNIFY、CLASSIFY×3、QUANTIFY/MERGE/DESEQ、MULTIQC）。
