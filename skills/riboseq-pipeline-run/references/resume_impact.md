# 流程修改对 -resume 的影响矩阵

Nextflow 的 `-resume` 按**任务签名（task hash）**命中缓存。任务签名包含：
进程名、模块脚本（main.nf）、任务命令引用的所有输入文件和参数、以及
**进程级指令**（cpus/memory/time/errorStrategy/maxForks 等在 Nextflow 22+
版本都参与签名）。改了什么决定重算多少。

## 影响矩阵

| 改动 | 重算范围 | 说明 |
|---|---|---|
| 改某模块 main.nf / bin 脚本 | 该模块全部任务 + 所有下游（DAG 级联） | 最常见的"改 bug 重算"场景 |
| 改 `-c` 配置里该进程的指令（cpus/memory/time/errorStrategy/maxForks） | 该进程的全部任务 | **maize 血泪教训**：中途改 errorStrategy → 3257 个 work dir 全部 hash 失效。改指令=改签名 |
| 改样本表（增删样本/改路径/改列值） | 受影响样本的整条下游链 | 被移除样本的 work dir 变成孤儿（占磁盘），需手动清理 |
| 改参考文件（--fasta/--gtf/--transcript_fasta/--contaminant_fasta） | 几乎全链 | 参考是比对起的源头输入 |
| 改参与任务命令的 CLI 参数（如 --ribotricer_phase_score_cutoff） | 引用该参数的模块 + 下游 | 参数进入任务命令 → 进签名 |
| 改不参与任何任务命令的参数（如 --max_memory 只在调度层面） | 不影响任务 hash | 但 --max_* 若作为进程指令处理仍会影响（见上一行） |
| 删除 result/ 输出文件 | **不影响缓存** | ⚠️ 常见误区：删了结果想重算，resume 不会重算也不会重新 publish（缓存任务不重新执行）——必须删对应 work dir |
| 删除 .nextflow/cache 的 LevelDB | work dir 全部变孤儿 | 不可恢复（`cache rebuild` 救不了，1.4TB 教训） |

## 修改后 resume 的标准剧本

1. **先想清楚改动属于矩阵哪一行**，评估重算面。
2. **小修 bug**：只删受影响任务的 work dir，resume 后只重算这些任务 + 下游：

   ```bash
   cd run/{project}/process
   find work/ -name ".command.sh" -exec grep -l "任务名正则" {} \; \
     | while read f; do rm -rf "$(dirname "$f")"; done
   # 再用 -resume 启动
   ```

3. **改动进了任务签名但你以为没有**：确认改的脚本真的是该任务执行的脚本
   （有些模块调 bin/ 下脚本，改 scripts/ 下同名文件无效）；确认后按上一步
   强制重算。
4. **改 config 指令**：没有"只重算一部分"的余地（该进程全部任务失效），
   要么接受，要么运行前改好。
5. **换样本表**：先 kill → 清理被移除样本的 work dir
   （`find work -maxdepth 2 -name '.command.log' | xargs grep -l "$sample"`）→ resume。
6. **对比 DAG**：改 workflow 后重新生成 `-with-dag flowchart.html` 与旧版对比，
   确认拓扑变化符合预期。
7. **删除 result 后想重算**：删 result 文件 + 删对应 work dir（否则 resume 静默
   跳过，结果文件也不会回来）。

## 铁律重述

运行期间不修改任何 config/代码/样本表——用 git commit hash 锁定版本；
非改不可时先 kill pipeline，改完按矩阵评估重算面再 resume。
