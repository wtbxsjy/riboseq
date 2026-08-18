---
name: riboseq-pipeline-run
description: >
  运行、resume、监控和排查本仓库 nf-core/riboseq pipeline 的 nextflow run 命令。
  涵盖标准启动模板（run/{project}/process/run_pipeline.sh）、各物种参数差异、
  nohup/-bg 后台运行与日志监控、会话锁/缓存/resume 失败/容器拉取失败/OOM 等
  常见报错的排查修复、容器构建与管理。
  当用户提到运行 pipeline、nextflow run、-resume、pipeline 报错/卡住/死锁/失败、
  清理 work 目录或会话缓存、构建 singularity/apptainer 容器时使用本 skill——
  即使只是说"帮我跑一下 rice"或"resume 卡住了"。
---

# riboseq pipeline 运行

运行入口是 `run/{project}/process/run_pipeline.sh`（由 riboseq-data-prep 的
prepare_workflow.py 生成）。本 skill 覆盖运行、后台监控、停止清理、排错和容器管理。

## 决策树

- **首次运行/正常重跑** → 第 1-2 节（标准模板 + 物种差异）
- **后台长期运行** → 第 3 节（nohup + 监控脚本）
- **卡住/报错/resume 失败** → 第 4 节速查 → `references/troubleshooting.md`
- **构建/更新容器** → 第 5 节 → `references/container_management.md`
- **跑完确认各阶段结果是否正常** → 第 8 节 + `scripts/check_pipeline_outputs.sh`
- **改了代码/config 想知道 resume 重算范围** → 第 6 节铁律 1 → `references/resume_impact.md`

## 1. 标准启动模板（rice 真实命令，逐段注释）

```bash
nextflow run /home/25119231r/riboseq/riboseq/main.nf \
    -profile singularity \
    -w /home/25119231r/riboseq/run/rice/process/work \
    --input /home/25119231r/riboseq/run/rice/scripts/samplesheet.csv \
    --outdir /home/25119231r/riboseq/run/rice/result \
    --aligner star \
    --max_memory '350.GB' --max_cpus 64 --max_time '96.h' \
    --save_reference true \
    --orfquant_container .../containers/orfquant_patched.sif \
    --unify_orf_container .../containers/unify_orf.sif \
    --gencode_orf_mapper_container .../containers/gencode_orf_mapper.sif \
    --ribowaltz_container .../containers/ribowaltz.sif \
    --price_container .../containers/gedi.sif \
    --fasta .../reference/rice.genome.fa \
    --gtf .../reference/rice.gtf \
    --transcript_fasta .../reference/rice.transcripts.fa \
    --contaminant_fasta .../reference/rice_final_contamination.fasta \
    --orf_classify_ensembl_dir .../reference/ensembl/Ens58_oryza_sativa \
    --skip_rpbp true \
    --skip_orfquant true \
    --skip_orf_classify_orfquant true \
    --ribotricer_phase_score_cutoff 0.1 \
    --skip_collect_qc_stats true \
    --orfquant_psite_correction true \
    --unify_orf_min_len 6 \
    --unify_orf_frame_merge_min_overlap 0.9 \
    --sorf_read_len_min 0 --sorf_read_len_max 0 \
    -c .../scripts/gencode_fast_24.config \
    -resume -bg \
    -with-report .../result/pipeline_report.html \
    -with-timeline .../result/timeline.html \
    -with-dag .../result/flowchart.html
```

要点：

- **脚本开头固定 `rm -f timeline.html flowchart.html pipeline_report.html`**——report 文件残留会导致 resume/重跑报 "already exists"，这是脚本内置的防线，手动构造命令时不要省。
- **多行反斜杠续行内禁止 `#` 注释**（会吞掉后续所有参数）；禁用参数 = 整行删除。
- **`bash run_pipeline.sh -resume` 无效**——脚本不接受参数，`-resume` 被丢弃会从头重跑。要 resume 就把 nextflow 命令复制出来手动加 `-resume`。
- `--sorf_read_len_min 0 --sorf_read_len_max 0` = 关闭 sORF 读长过滤（默认是 28-30）。
- `--skip_orf_classify_orfquant` **不在 nextflow_schema.json 里**，会打 `WARN: invalid input values`——无害，属预期（ORFquant 分类临时停用）。
- `-c` 追加的进程级 config 见 `references/run_commands_by_species.md`。

## 2. 各物种差异速查

| 项目 | 差异点 |
|---|---|
| rice | 模板本身：star、350GB/64cpu/96h、min_len 6、ribotricer cutoff 0.1、无长度过滤 |
| maize | `--input samplesheet_hq.csv`（HQ 过滤后）；`-c orfquant_override.config`（ORFQUANT_RUN cpus=16 + ORF_QC errorStrategy='ignore'）；中途换过 orfquant_mirai.sif |
| wheat | **唯一 hisat2 实跑**：`--aligner hisat2`、`--bam_csi_index`、150GB/42cpu/48h、`--unify_orf_min_len 24`、`--orf_classify_mode orf_type`、`--unify_orf_merge_tolerance 3 --unify_orf_min_overlap 0.5`、`-process.maxForks=8` |
| soybean | star、min_len 24、orf_classify_mode orf_type、resume |
| mouse_tissue | `--gencode_classify_impl indexed_fast`、`-process.maxForks=8` |
| human_GSE158930 | 双基因组：`--pathogen_fasta/--pathogen_gtf/--pathogen_contig_pattern '^NC_045512'`、`--sorf_read_len_min 28 --sorf_read_len_max 33`、`--ribowaltz_read_lengths [28,29,30,31,32,33]`、`--ribotish_fail_on_empty true` 等 fail_on_empty 系列、`--rpbp_container rpbp.sif` |
| mouse_Mucosal_Immunity | `--igenomes_ignore true --igenomes_base ''`、`--skip_rpbp true --skip_price true` |

完整命令原文 → `references/run_commands_by_species.md`。

**BAM 输入模式**：samplesheet 用 `bam`/`bam_index` 列，strandedness 不能 auto，UMI 去重和 RiboCode 自动跳过。

## 3. 后台运行与监控

```bash
# 后台启动（标准姿势）
cd ~/riboseq/run/rice/process && nohup bash run_pipeline.sh > run_rice.log 2>&1 & echo "PID: $!"

# 带时间戳的 resume 日志
nohup bash run_pipeline.sh > pipeline_resume_$(date +%Y%m%d_%H%M%S).log 2>&1 &
```

**自带监控脚本**（后台轮询直到完成或超时，增量展示进度）：

```bash
bash skills/riboseq-pipeline-run/scripts/monitor_pipeline.sh \
  ~/riboseq/run/rice/process/.nextflow.log --interval 300 --timeout-hours 96
```

**手动进度查询惯用式**：

```bash
grep "Submitted\|Cached\|ERROR" .nextflow.log | grep -v DEBUG | tail -5
grep -c "Submitted" .nextflow.log      # 累计提交数
```

**判死锁前先查任务进程**：`"no new submissions in many hours"` 不一定死锁——
可能是某任务（如 ORFquant 多线程、EXPRESSION_QUANT 大扫描）自身卡住，
`ps aux | grep` 对应任务进程确认后再处置。

## 4. 停止与清理

```bash
# 停止
pkill -f "nextflow.*rice"          # 按项目名
kill -9 $(ps aux | grep "nextflow.*one.jar" | grep -v grep | awk '{print $2}')

# 清理会话锁（force-kill 后 resume 前必做）
rm -rf .nextflow/cache/<session-uuid>/ && rm -f .nextflow.pid
rm -f .nextflow/cache/*/db/LOCK
```

**自带清理脚本**（默认干跑预览，`--apply` 才真删）：

```bash
bash skills/riboseq-pipeline-run/scripts/cleanup_nextflow.sh \
  ~/riboseq/run/rice/process --result-dir ~/riboseq/run/rice/result \
  --clear-lock --failed-only --task-pattern "UNIFY_ORF|ORF_QC" --apply
```

按任务名删 work dir 强制重算的惯用式：

```bash
find work/ -name ".command.sh" -exec grep -l "UNIFY_ORF\|ORF_QC" {} \; | while read f; do rm -rf "$(dirname "$f")"; done
find work/ -name ".exitcode" -exec grep -l "^[^0]" {} \; | while read f; do rm -rf "$(dirname "$f")"; done
```

## 5. 排错速查

| 症状 | 处置 | 详情 |
|---|---|---|
| resume 被锁/卡住不动 | 删会话缓存 + LOCK + .nextflow.pid | troubleshooting §1 |
| report 文件 "already exists" | 删 timeline/flowchart/pipeline_report.html | troubleshooting §2 |
| 改了代码 resume 还是旧结果 | 按任务名删对应 work dir 强制重算 | troubleshooting §3 |
| 换 samplesheet 后 resume 失败 | 清理被移除样本的 work dir | troubleshooting §3 |
| PoisonPill / DataflowBroadcast 报错 | 多为瞬时，清理缓存重跑 | troubleshooting §4 |
| WARN: undefined parameter | config 里初始化默认值 | troubleshooting §5 |
| PRICE 容器 Wave 拉取失败 (status 255) | 改本地 gedi.sif：`--price_container /path/gedi.sif` | troubleshooting §6 |
| 任务 exit 137（OOM） | `dmesg | grep -i "oom\|killed"` 确认后降内存/并行度 | troubleshooting §7 |
| STAR exit 143 | 多为受牵连（如 riboWaltz 注释失败），删 work dir 重跑 | troubleshooting §7 |
| EXPRESSION_QUANT futex_wait 数小时 | 大样本正常慢，勿误判死锁 | troubleshooting §8 |

## 6. 铁律（每条都有惨痛教训支撑）

1. **运行中途不要改任何 config/代码**——maize 案例改 errorStrategy 导致 3257 个 work dir 全部 hash 失效重跑数天。用 git commit hash 锁定版本。**改了什么会重算多少、怎么只重算受影响任务** → `references/resume_impact.md`（影响矩阵 + 修改后 resume 剧本）。
2. **不要动 `.nextflow/cache` 里的 LevelDB**——删了 work dir 就变孤儿数据，`cache rebuild` 也救不回来（1.4TB 教训）。
3. **清理前先干跑预览**（cleanup_nextflow.sh 默认 dry-run）。
4. **resume 前先清锁**：`rm -f .nextflow/cache/*/db/LOCK`。

## 7. 容器管理

- 构建：`apptainer build --fakeroot -F out.sif containers/Singularity.xxx.def`（卡住就去掉 `--fakeroot` 重试；旧源码被缓存时 `apptainer cache clean --force`）
- 免重建技巧：`NXF_SINGULARITY_BINDPATH="/host/src:/opt/pkg"` 运行时 bind 开发源码进容器
- conda 双模式：`-profile test,conda`（ORFquant 用 conda 绕开 Apptainer 问题）
- 完整命令与 sandbox/OOM 限制 → `references/container_management.md`

## 8. 运行后结果检查

每个子步骤都有明确的正常输出形态。一键体检：

```bash
bash skills/riboseq-pipeline-run/scripts/check_pipeline_outputs.sh \
  run/rice/result --samples-file run/rice/scripts/samplesheet.csv --expected-orfs 370000
# 退出码：0 通过 / 1 有警告 / 2 有失败
```

检查项：阶段目录齐全 → riboseqc P_sites bedgraph 数 == 样本数 → 各 ORF 工具
per-sample 输出数与空文件 → metadata 行数（±50% 期望值容差）→ expression
summary 行数 == metadata 行数 → ORF_QC confidence 行数一致性 → GENCODE biotype
分布（**lncRNA >90% = 蛋白 header 不匹配坑**）。

基准健康数字（rice 23 样本：369,939 ORFs / metadata 271MB；maize 97 样本 ~660K；
human 62 样本 956,971）与逐阶段"异常信号解读"→ `references/results_checklist.md`。

**注意**：缺文件时先判断是"被 skip 的参数"还是"失败被 errorStrategy='ignore'
静默吞掉"——后者才是问题。

## 运行完成后

- 结果在 `result/`（orf_unification / orf_classification / riboseqc / expression ...）
- ORF 定量与过滤 → **riboseq-quant-te** skill；手动 unify/分类 → **riboseq-orf-analysis** skill
