# 各物种运行命令与参数差异（真实运行记录）

磁盘原件：`run/{project}/process/run_pipeline.sh`。以下按项目记录与 rice 模板的差异，
完整命令直接读对应脚本。机器规格：56 核 / 188GB，local executor，无 SLURM。

## rice（23 样本，模板基准）

- aligner star、350.GB / 64 cpu / 96.h
- `--skip_rpbp --skip_orfquant --skip_orf_classify_orfquant`
- `--ribotricer_phase_score_cutoff 0.1`（植物低深度适配）
- `--unify_orf_min_len 6 --unify_orf_frame_merge_min_overlap 0.9`
- `--sorf_read_len_min 0 --sorf_read_len_max 0`（关闭读长过滤）
- `-c scripts/gencode_fast_24.config` + `-resume -bg` + 三个 report
- 最终统计：succeeded=6, cached=743, peakMemory=136GB, 8h15m

## maize（97 HQ 样本）

- `--input samplesheet_hq.csv`（salmon 比对率 >20% 过滤后）
- 额外 `-c process/orfquant_override.config`：

```groovy
process {
    withName: 'ORFQUANT_RUN' { cpus = 16; errorStrategy = 'ignore' }
    withName: 'ORF_QC'       { errorStrategy = 'ignore' }
}
```

- 后期换 `orfquant_mirai.sif`（v1.3.2 mirai 并行后端）；无 `-resume/-bg`，用 nohup
- 教训：运行中途改 errorStrategy → hash 全失效，3257 个 work dir 报废

## wheat（hisat2 唯一实跑）

- `--aligner hisat2 --bam_csi_index`、150.GB / 42 cpu / 48.h
- `--unify_orf_min_len 24`、`--orf_classify_mode orf_type`
- `--unify_orf_merge_tolerance 3 --unify_orf_min_overlap 0.5`
- `-resume -process.maxForks=8`
- `--orf_classify_ensembl_dir ~/riboseq/run/wheat/reference/ensembl_lib/`

## soybean / mouse_tissue / mouse_tissue_old

- soybean：star、150GB/42cpu/48h、min_len 24、orf_classify_mode orf_type、-resume
- mouse_tissue：`--gencode_classify_impl indexed_fast`、`-process.maxForks=8`
- mouse_tissue_old：+ `--extra_unify_orf_predictions_args "--stats-mode preload"`、
  `--merge_replicates`、`-process.maxForks=6`

## human_GSE158930（双基因组 SARS-CoV-2）

- `--pathogen_fasta SARS-CoV-2.genome.fa --pathogen_gtf <pathogen.gtf>`
- `--pathogen_contig_pattern '^NC_045512'`（CLAUDE.md gotcha 16：双基因组建议预拼接后走 --fasta/--gtf）
- `--sorf_read_len_min 28 --sorf_read_len_max 33`
- `--ribowaltz_read_lengths [28,29,30,31,32,33]`
- `--ribotish_fail_on_empty true`（+ ribotricer/orfquant/rpbp/price/ribocode 同系列）
- `--rpbp_container rpbp.sif`

## mouse_Mucosal_Immunity

- `--igenomes_ignore true --igenomes_base ''`、`--skip_rpbp true --skip_price true`
- `nohup nextflow run ... -bg`

## 进程级 config 覆盖（`-c` 追加）

`gencode_fast_24.config`（rice/maize 通用）：

```groovy
process {
    withName: 'CLASSIFY_ORFS_GENCODE' { cpus = 16; memory = 64.GB; time = 48.h }
}
```

用法：大任务单独放宽资源、或 `errorStrategy='ignore'` 让非关键任务失败不中断
整条 pipeline（ORF_QC 案例）。⚠️ 运行中途改这类文件 = hash 失效（见铁律 1）。

## BAM 输入模式

- samplesheet 列：`sample,bam,bam_index,strandedness,type`
- strandedness 必须显式（forward/reverse/unstranded），不能 auto
- UMI dedup 与 RiboCode 自动跳过；sORF 过滤仍生效
