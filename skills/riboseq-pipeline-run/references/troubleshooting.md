# Pipeline 运行排错决策树

每条含：症状 → 原因 → 修复。命令中的路径按实际项目替换。

## 1. 会话锁：force-kill 后 resume 卡住/报锁错误

**症状**：`-resume` 后 pipeline 不启动或报 session cache 锁冲突。

**原因**：Nextflow 被 `kill -9` 后 `.nextflow/cache/<uuid>/db/LOCK` 残留。

**修复**（resume 前的标准前置动作）：

```bash
cd run/{project}/process
rm -rf .nextflow/cache/<session-uuid>/     # 明确知道是哪个坏 session 时
rm -f .nextflow.pid
rm -f .nextflow/cache/*/db/LOCK
```

- 更粗暴：`rm -rf .nextflow`（放弃全部缓存，全新跑）
- resume 失败有时是**新建了空 session**：删除新的 session 目录、用原始 session
  缓存目录 resume 才能命中缓存
- ⚠️ 绝不要删 LevelDB 后又想靠 `cache rebuild` 恢复——LevelDB 文件删了
  work dir 就是孤儿数据，1.4TB 教训，不可恢复

## 2. report 文件 "already exists"

**症状**：resume/重跑时报 timeline/flowchart/pipeline_report.html 已存在。

**修复**：`rm -f result/{timeline,flowchart,pipeline_report}.html`
（run_pipeline.sh 开头已内置，手动构造命令时不要省）。

## 3. resume 复用旧结果 / 换样本表后失败

**症状 A**：修了代码后 resume，任务仍用旧 `.command.sh` 跑（结果没变）。

**修复**：按任务名删除对应 work dir 强制重算：

```bash
cd run/{project}/process
find work/ -name ".command.sh" -exec grep -l "UNIFY_ORF\|ORF_QC" {} \; \
  | while read f; do rm -rf "$(dirname "$f")"; done
# 或按退出码删所有失败任务
find work/ -name ".exitcode" -exec grep -l "^[^0]" {} \; \
  | while read f; do rm -rf "$(dirname "$f")"; done
```

**症状 B**：换 samplesheet（如 samplesheet_hq.csv）后 resume 报 DAG 路径解析错误。

**修复**：先 kill，清理被移除样本的残留 work dir
（`find work -maxdepth 2 -name '.command.log' | xargs grep -l "$sample"` 定位后 rm），
再 resume。

**原理**：改 workflow 代码后模块 hash 变化，`-resume` 全失效是预期行为——要么
接受重跑，要么只删受影响任务（`-c` 进程级 config 不改变模块 hash，可安全追加）。

## 4. PoisonPill / DataflowBroadcast / channel 死锁

**症状**：`PoisonPill` 传入进程崩溃、`.into{}` 在缓存任务上报
`DataflowBroadcast` 不可用、ORF_QC 无限等待空 channel。

**处置**：
- 多为瞬时问题：清会话缓存重跑一次
- 历史根因已修（ORF_QC 改为磁盘 glob + `optional: true` 占位文件），
  若再遇到 `.into{}`/`.ifEmpty()` 报错，优先怀疑是 resume 的缓存任务与新
  workflow 代码不匹配 → 按 §3 删对应任务 work dir

## 5. WARN: undefined parameter

**症状**：`WARN: Access to undefined parameter 'xxx' -- Initialise it to a default value`

**修复**：在 `nextflow.config` 的 `params {}` 里加默认值。
注意 Nextflow 26.x 布尔参数 CLI 传入是字符串（schema 需
`"type": ["boolean","string"]`，见 CLAUDE.md gotcha 15）。

## 6. PRICE 容器拉取失败

**症状**：

```
Failed to pull singularity image ... oras://community.wave.seqera.io/library/gedi_price:latest
status : 255  hint: Try and increase singularity.pullTimeout
```

**修复**：不走 Wave，用本地镜像：

```bash
--price_container /home/25119231r/MyDrive/sORF_Discovery_Project/containers/gedi.sif
```

（GEDI v1.0.6 另有参数格式坑：只认 `-reads/-prefix/-genomic`，`-genomic` 要
OML 绝对路径；修复后 18,091 ORFs/sample。）

## 7. 任务 exit 137（OOM）与 exit 143（SIGTERM）

**症状**：任务 exit code 137 / 143。

**判定 OOM**：

```bash
dmesg | grep -i "oom\|killed process"
journalctl -k | grep -i oom
cat /sys/fs/cgroup/memory.max    # sandbox 内存上限可能远小于宿主
```

**已知案例**：
- ORFquant sandbox 内 visible memory 仅 ~54GB（宿主 188GB），16 mirai daemon ×
  2.5GB 全灭 → 最终方案：sandbox + mclapply n_cores=1 + 每样本串行 +
  MAX_PARALLEL=8 跨样本并行（见 container_management.md）
- STAR exit 143：多为被 riboWaltz 注释构建失败牵连，删 work dir 重跑
- exit 137 修复后记得确认 `errorStrategy` 是否把失败静默了

## 8. "Pipeline seems deadlocked" 误判

**症状**：`no new submissions in many hours`。

**处置**：先查具体任务进程（`ps aux | grep`）：
- EXPRESSION_QUANT 主进程 futex_wait 等 worker 5 小时属大样本正常慢（已优化为
  27h → <10s 的 streaming 单次扫描，旧版本才会遇到）
- ORFquant 多线程卡住是任务自身问题，不是 pipeline 死锁
- 确认真死锁再考虑 §1/§4 的清理

## 9. 其他速查

| 症状 | 修复 |
|---|---|
| resume 时 report 残留 | §2 |
| `-resume` 从头重跑 | 检查是否误用 `bash run_pipeline.sh -resume`（参数被丢弃）；手动复制命令加 `-resume` |
| 结果目录中间产物脏了 | 清 `result/orf_unification orf orf_classification expression_quant` 强制重跑下游 |
| 全量重跑 | `rm -rf process/work process/.nextflow`（谨慎：丢全部缓存） |
