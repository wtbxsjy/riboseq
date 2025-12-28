# 单工具 Singularity 验证脚本（从 BAM 开始）

这组脚本用于在跑整条 Nextflow 流程前，先逐个验证每个工具在你的环境里能用（镜像版本与 pipeline module 保持一致）。

## 通用要求

- 需要 `singularity` 可用（或 Apptainer 兼容命令）。
- 建议在一个“干净的工作目录”里运行脚本，输出会写到当前目录下的子目录。
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

## 最小运行示例

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
