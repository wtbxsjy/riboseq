# 单工具独立运行（从 BAM 开始，不跑 pipeline）

脚本位置：`scripts/singularity_single_tool_tests/`（01..20 编号脚本）。每个脚本注释里
都标了 `Mirrors: modules/local/xxx`，容器版本/参数与 pipeline module 一致。
`00_run_order_example.sh` 是把 01..11 按依赖顺序串联的模板（默认 `DRY_RUN=1` 只打印）。

## 依赖关系图

```
BAM + .fai ──→ 01 sorf_bam_filter ──→ {sample}.sorf.filtered.bam   ← 所有预测工具的共同输入
GTF + FASTA ──→ 02 riboseqc_prepareannotation ──→ *_Rannot          （一次性）
filtered BAM + Rannot ──→ 03 riboseqc_analysis ──→ *_for_ORFquant + *_P_sites_calcs
                                  │
                                  └─→ 04 orfquant_run   ⚠️ 必需 02 的 *_Rannot + 03 的 *_for_ORFquant
filtered BAM + GTF ──→ 05 ribotish_quality ──→ {sample}.para.py
                                  │
                                  └─→ 06 ribotish_predict (+FASTA)  ⚠️ 必需 05 的 .para.py
GTF + FASTA ──→ 07 ribotricer_prepareorfs ──→ {prefix}_candidate_orfs.tsv  （一次性，与样本无关）
                                  │
filtered BAM + index ──→ 08 ribotricer_detectorfs ──→ *_translating_ORFs.tsv
FASTA + GTF + rRNA ──→ 09 rpbp_prepare_genome ──→ *.orfs-genomic/exons.bed.gz  （一次性）
                                  │
filtered BAM + bed.gz ──→ 10 rpbp_predict ──→ bayes-factors.bed.gz
BAM + GTF + FASTA ──→ 11 ribocode_detect   ⚠️ 推荐 transcriptome BAM
12 = 03+04 组合（filtered BAM 一条龙出 ORFquant 结果）
13 orfquant_prepareannotation（Rannot 生成的 ORFquant 侧注释准备）
```

## 前置条件速查表

| 工具 | 脚本 | 硬前置 | 注意 |
|---|---|---|---|
| sORF BAM 过滤 | 01 | BAM + .fai | 过滤规则与 pipeline sorf_* 参数一致（unique/contig regex/28-30nt/flags） |
| ORFquant | 04 | **RiboseQC 输出 `*_for_ORFquant` + `*_Rannot`** | ORFquant 不吃 BAM；没有 RiboseQC 前置产物就跑不了 |
| Ribo-TISH | 05→06 | predict 必需 quality 的 `*.para.py` | 两步串联，缺一不可 |
| Ribotricer | 07→08 | detect 必需 prepare 的候选 index | index 一次性、可跨样本复用 |
| rp-bp | 09→10 | 必需 rRNA FASTA（prepare-genome 用） | index 一次性；predict 需要 orfs-genomic + orfs-exons 两个 bed.gz |
| RiboCode | 11 | GTF + FASTA + 转录组 BAM（推荐） | 低深度数据 periodicity 不足失败属正常 |
| PRICE | 无单工具脚本 | GTF/OML | 只能走 pipeline module 或手写 singularity 命令（见下） |

通用要求：BAM 已排序已建索引（脚本会自动补 `samtools index`/`faidx`）；在干净工作目录
运行；镜像缓存到 `scripts/singularity_single_tool_tests/containers/`；WSL/OneDrive 路径
用 `export BIND_EXTRA="/mnt:/mnt"`。

## 各工具详解

### ORFquant（04/12/13）⚠️ 前置最重

```bash
./04_orfquant_run.sh \
  --sample S1 \
  --for-orfquant out_riboseqc_analysis/S1_for_ORFquant \
  --annotation out_riboseqc_annot/REPLACE_ME_Rannot \
  --fasta genome.fa --cpus 8
```

- 输入：RiboseQC 的 `*_for_ORFquant`（03 产出）+ `*_Rannot`（02 产出）+ genome FASTA
- 容器解析顺序：`--container` 指定 → 当前目录 `orfquant.sif` → 自动 pull
  `orfquant:1.1.0--r40_1`（R 4.0）；也可 `--orfquant-pkg ORFquant-1.02.tar.gz` 避免
  运行时从 GitHub 下载源码
- 输出：`{SAMPLE}_final_ORFquant_results*`
- 脚本内置：namespace 冲突修复（conflicted/Position/combine）+ 并行模式失败自动
  n_cores=1 重试 + 最终化错误抑制（详情见 `orfquant_saga.md`）
- `12_riboseqc_orfquant_from_filtered_bam.sh` = 03+04 一条龙；`13` 只做注释准备
- P-site 偏移校正链：03 → `19_extract_rl_cutoff.sh` → `20_prepare_for_orfquant_corrected.sh`
  → 04（用校正后的 for_ORFquant 跑，提高准确性）

### Ribo-TISH（05→06）

```bash
./05_ribotish_quality.sh --sample S1 --bam S1.sorf.filtered.bam --gtf annot.gtf
./06_ribotish_predict.sh --sample S1 --bam S1.sorf.filtered.bam \
    --gtf annot.gtf --fasta genome.fa --ribopara out_ribotish_quality/S1.para.py
```

- quality 产出 offset 参数 `S1.para.py`；predict 产出 `S1_pred.txt`
- 对应 pipeline：Ribo-TISH 每样本 `{sample}_pred.txt`

### Ribotricer（07→08）

```bash
./07_ribotricer_prepareorfs.sh --gtf annot.gtf --fasta genome.fa --prefix genome
./08_ribotricer_detectorfs.sh --sample S1 --bam S1.sorf.filtered.bam \
    --index out_ribotricer_prepareorfs/genome_candidate_orfs.tsv \
    --stranded forward
```

- prepare-orfs 的 prefix 与样本无关（索引名）；产出 `${prefix}_candidate_orfs.tsv`
- detect 产出 `S1_translating_ORFs.tsv`
- ⚠️ 真实 flag 是 `--index`；该目录 README.md 示例里写的 `--ribotricer-index` 是过时写法
- unstranded 模式在 ribotricer 里有已知问题，脚本默认不传

### rp-bp（09→10）

```bash
./09_rpbp_prepare_genome.sh --fasta genome.fa --gtf annot.gtf --rrna rrna.fa
./10_rpbp_predict.sh --sample S1 --bam S1.sorf.filtered.bam \
    --orfs-genomic out_10/transcript-index/genome.orfs-genomic.bed.gz \
    --orfs-exons  out_10/transcript-index/genome.orfs-exons.bed.gz
```

- rRNA FASTA 是 prepare-genome 的必需输入
- predict 产出 `bayes-factors.bed.gz`（pipeline 中对应 `--run_rpbp` 的输出）

### RiboCode（11）

```bash
./11_ribocode_detect.sh --sample S1 --bam S1.toTranscriptome.bam \
    --gtf annot.gtf --fasta genome.fa --stranded forward
```

- **更推荐 transcriptome BAM**（pipeline FASTQ 模式会额外生成）；genome BAM 也能跑但效果差
- 容器：ribocode:1.2.11 + samtools 双镜像；低深度数据失败常见
- `--args` 可透传 RiboCode_onestep 额外参数

### PRICE（无单工具脚本）

pipeline module `modules/local/price/` 的关键点（手写 singularity 命令时照抄）：
v1.0.5 只认 `-reads/-prefix/-genomic`；`-genomic` 要 **OML 绝对路径**
（`$(realpath ${genome_name}.oml)`，GEDI 按名字查找）；容器 gedi_price。

## 完整串联模板

```bash
cd scripts/singularity_single_tool_tests
# 填好 BAM/FASTA/GTF/FAI/RRNA 后：
DRY_RUN=0 CPUS=16 ./00_run_order_example.sh
# RANNOT / FOR_ORFQUANT 变量要指向 02/03 实际产出的 *_Rannot / *_for_ORFquant
```

## 手动 unify 与分类（17/18 脚本）

单工具跑完后，不需要 pipeline，用这两个脚本即可完成统一与分类：

### 17_unify_predictions.sh（wraps `scripts/unify_orf_predictions.py`）

```bash
./17_unify_predictions.sh \
  --gtf annot.gtf --fasta genome.fa \
  --ribotish "s1_pred.txt s2_pred.txt" \
  --ribotricer "s1_orfs.tsv" --orfquant "s1_orfquant.gtf" \
  --output ./results/unified
```

- 只支持 Ribo-TISH / Ribotricer / ORFquant 三种输入（脚本无 `--ribocode`/`--price`；
  需要时用 `run_orf.py unify`，见 `manual_bypass.md`）
- P-site 统计：`--bedgraph-dir`（RiboseQC bedgraph 目录）+ `--sample-list "s1,s2"`
- 合并参数：`--frame-merge-min-overlap`（默认 0.9）、`--min-overlap`（默认 0.5）、
  `--min-len`（默认 10）、`--no-frame-merge`
- ⚠️ 历史坑（2026-08-18 修复）：旧版脚本传 `--merge-tolerance`（bp）和
  `--no-overlap-group`，重构后的 unify_orf_predictions.py 已无这两个 flag，
  argparse 严格解析会直接报错；现已改为 `--frame-merge-min-overlap`，
  `--no-overlap-group` 接受但忽略并打 WARN
- 容器：python:3.9 镜像 + 运行时 `pip install biopython pyfaidx`

### 18_classify_orfs.sh（wraps `scripts/classify_orfs_wrapper.py`）

```bash
./18_classify_orfs.sh --mode gencode --input ./results/unified \
  --ensembl-dir ./refs/Ens58 --output-dir ./class_res
./18_classify_orfs.sh --mode orf_type --input ./results/unified \
  --gtf annot.gtf --output-dir ./class_res
```

- `--mode gencode|orfquant|orf_type`；`--input` 传 17 的输出 prefix
  （自动找 `.bed`/`.gtf`/`.metadata.tsv`）
- gencode 需要 `--ensembl-dir`（5 文件目录，见 riboseq-data-prep）+ `--fasta`
  （input 为 prefix 时）；orfquant/orf_type 需要 `--gtf`
- 脚本内部把 `--output-dir` 正确转成 wrapper 的 `--output_dir`，直接照用即可

### 与 run_orf.py 的关系

17/18 走原版脚本（无 DuckDB 加速、无 RiboCode/PRICE 输入）；样本多/数据大或
需要 RiboCode/PRICE 时用 SKILL.md §4 的 `run_orf.py` + `run_orf_in_container.sh`。

## 结果进入 GENCODE 注释链

另一个路径：`14/15_*_to_gencode.sh` 格式转换 → `16_gencode_orf_mapper.sh`
（需要 5 文件 Ensembl 目录，见 riboseq-data-prep）。
