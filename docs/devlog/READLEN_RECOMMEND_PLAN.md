# 计划：在流程中增加"原始 BAM 读长评估与推荐"步骤

> 2026-09-11 制定。背景分析见 `READLENGTH_THRESHOLD_ANALYSIS.md`。

## 0. 现状与关键发现

**已存在但未启用**：`workflows/riboseq/main.nf:877-884` 已有 `RIBOSEQC_PREFILTER`，跑在
**未过滤 BAM**（`ch_bams_for_prefilter`，line 425）上，由 `params.run_prefilter_qc`
（`nextflow.config:193`，默认 `false`）控制，输出发布到 `riboseq_qc/riboseqc/prefilter`
但**不接任何下游消费者**。

**读长假设存在于两处**（改动必须同时覆盖）：

| 位置 | 当前值 | 影响 |
|---|---|---|
| `nextflow.config:108-109` `sorf_read_len_min/max` | 28 / 30 | `SORF_BAM_FILTER` 抛弃读段 |
| `nextflow.config:169` `ribowaltz_read_lengths` | `[28,29,30]` | 逐读长 P-site offset 只在窗内计算 |

**为什么必须修**（实测，见分析文档 §3）：

1. 28-30 丢弃小鼠 63%、人类 57% 的已比对读段，而 RiboseQC/riboWaltz 跑在过滤**之后**，
   结构上无法评估被丢弃的部分（两者的逐读长表都只有三行）。
2. 小鼠 27 nt：20-27 万 CDS 内读段、框偏好 55% —— **高于 RiboseQC 自身选读长用的 50% 门槛**，
   却被硬阈值整体丢弃。
3. RiboseQC 在**原始 BAM** 上直接算会**误导**：小鼠 31 nt 在原始 BAM 上框偏好 94.3%，
   加 NH:i:1 过滤后跌到 55% —— 94% 是多比对假象（31 nt 唯一比对率仅 34%）。
   **因此推荐必须建立在"唯一比对+去污染 contig"之后、读长过滤之前的数据上。**

## 1. 目标与验收标准

**目标**：流程自动产出每样本的读长质量表与一个队列级推荐窗口，可选择性回灌到过滤器。

**验收标准**：
- [ ] 新项目跑完即得 `readlength_recommendation.tsv`（逐样本 + 队列级），无需手工分析
- [ ] 推荐数字与独立实现（`post_analysis/scripts/frame_preference_by_length.sh`）在 28/29/30
      上偏差 < 5 个百分点
- [ ] `--sorf_read_len_auto false`（默认）时，过滤结果与当前流程**逐字节一致**
- [ ] MultiQC 报告中可见读长质量表
- [ ] 全流程 `nf-test --profile debug,test,docker` 通过

## 2. 设计

### 2.1 数据流（新增部分）

```
align (未过滤 genome BAM)
  │
  ├─► SORF_BAM_FILTER_STAGE1 ──┬─► stage1.bam (唯一比对 + 去 contig，**不切读长**)
  │   (唯一+contig)             │        │
  │                             │        └─► RIBOSEQC_ANALYSIS (prefilter, 已有)
  │                             │                 │
  │                             │                 └─► P_sites_calcs（全读长）
  │                             │                          │
  │                             │                          └─► RECOMMEND_SORF_READLENGTH
  │                             │                                   │
  │                             │                       队列级读长窗口 (tsv)
  │                             │                                   │
  │                             └─► SORF_BAM_FILTER_STAGE2 ◄─────────┘
  │                                 (读长切分)  ← 关闭时用 params 常量（行为不变）
  │                                        │
  └────────────────────────────────────────┴─► 其余流程不变
```

要点：**stage1 的输出同时喂给 prefilter RiboseQC 和 stage2**，因此读长判据来自
"干净但未切读长"的 BAM —— 正是实测验证过的口径。

### 2.2 为什么拆成两个进程

Nextflow 中进程不能消费自己的输出。stage1 写一次 BAM，stage2 再读一次做长度切分，
比"先跑一次过滤、再跑一次不切长度的"省一半 I/O。

### 2.3 决策点（需要你拍板）

| # | 决策 | 选项 | 建议 |
|---|---|---|---|
| D1 | 是否复用已有 `RIBOSEQC_PREFILTER` | 复用 / 新建模块 | **复用**（零重复实现） |
| D2 | 推荐的应用方式 | 仅报告 / 自动回灌 | **默认仅报告**（`--sorf_read_len_auto` 打开回灌）。理由：自动改 ORF 集会让已发表分析不可复现 |
| D3 | 窗口粒度 | 逐样本 / 队列级 | **队列级**。逐样本窗口会让不同样本覆盖不同的读长集合，破坏跨样本可比性 |
| D4 | 合格判据 | — | `frame_preference >= 50`（RiboseQC 自己的门槛）**且** `pct_map >= 1%`（排除长尾噪声读长） |
| D5 | 窗口选取 | — | 覆盖 ≥90% 合格读段的**最小连续区间**；若合格读长不连续，取最大连续段并在报告中告警 |
| D6 | 队列聚合规则 | — | 队列窗口 = 各样本推荐窗口的**并集**（保证无样本被过度剥夺），并在表中给出每样本的独立推荐值供人工复核 |

### 2.4 新增/修改文件清单

**新增**

| 文件 | 内容 |
|---|---|
| `modules/local/recommend_sorf_readlength/main.nf` | 进程定义：收集全部 `*_P_sites_calcs` → 输出推荐表 |
| `bin/recommend_sorf_readlength.py` | 逻辑（纯标准库，便于 nf-test） |
| `modules/local/recommend_sorf_readlength/environment.yml` | python 环境 |
| `modules/local/recommend_sorf_readlength/tests/` | nf-test（含"低质量样本不产生空窗口"边界） |

**修改**

| 文件 | 改动 |
|---|---|
| `modules/local/sorf_bam_filter/main.nf` | 增加 `val stage`（`'unique'` / `'length'`）；stage1 跳过读长判断，stage2 只做读长判断。输出名区分 |
| `workflows/riboseq/main.nf:433-455` | 拆成两段调用 + 接入推荐通道 |
| `workflows/riboseq/main.nf:875-884` | prefilter 的输入由 `ch_bams_for_prefilter` 改为 stage1 BAM；`run_prefilter_qc` 与新的 `run_readlen_recommendation` 解耦 |
| `nextflow.config` | 新增 4 个参数（见下） |
| `nextflow_schema.json` | 同步 4 个参数（**注意 gotcha 15**：布尔必须 `["boolean","string"]` + enum） |
| `conf/modules.config` | `RECOMMEND_SORF_READLENGTH` 的 publishDir / ext.args |
| `assets/multiqc_config.yml:152` `custom_data:` | 新增 readlength 段落 |
| `CLAUDE.md` | 新增 gotcha：读长假设的两处耦合 + 原始 BAM 评估的多比对陷阱 |

**新参数**（默认值均保持现有行为）：

```
run_readlen_recommendation = false        // 是否执行推荐步骤
sorf_read_len_auto         = false        // 推荐值是否回灌到 stage2
readlen_frame_pref_min     = 50           // 合格判据：框偏好下限 (%)
readlen_pct_map_min        = 1            // 合格判据：读长占比下限 (%)
```

### 2.5 回灌的实现要点（避免哈希失效扩大化）

`SORF_BAM_FILTER` 的读长入参当前是 `params.sorf_read_len_min/max`（`val`）。
实现时应**替换为 channel 入参**（`Channel.value(...)` vs 推荐通道），而**不是新增输入槽**：

- 新增输入槽 → 输入元组形状改变 → 所有任务的 hash 改变 → 已跑项目整条 DAG 重算。
- 替换为同值通道 → 值相同时 hash 理论上不变，可保住既有缓存。

⚠️ **此点必须在实施时用小规模数据实测确认**（跑两次 `-resume`，比对
`.nextflow.log` 的 `Cached process` 数量），不可默认成立。

### 2.6 成本

| 项 | 实测/估算 |
|---|---|
| 现有 `RIBOSEQC_ANALYSIS`（已过滤 BAM） | 2m36s – 8m31s / 样本（小鼠） |
| 新增 prefilter RiboseQC（stage1 BAM，读段数约 4×） | 估 10–30 min / 样本，`label 'process_medium'`（6 cpu / 36 GB） |
| 新增 stage1 进程 | 与现有 `SORF_BAM_FILTER` 相当（小鼠实测 26–48 min / 样本） |
| `RECOMMEND_SORF_READLENGTH` | 秒级 |

净增约 **1.5–2× 单样本前处理时间**。对 12 样本队列 = 数小时（可并行）。

### 2.7 风险与对策

| 风险 | 对策 |
|---|---|
| **已有项目的 `-resume` 缓存全部失效** | 明确定性为 breaking change；已完成的 GSE120762 / PRJEB26593 **不重跑流程**，改用独立脚本（已完成）；本改动只惠及新项目 |
| 推荐给不出窗口（样本质量普遍差，如 PRJEB26593 最佳仅 52%） | 判据失败时回退到 `params.sorf_read_len_min/max` 并在 MultiQC 高亮告警，**不得**产生空窗口导致 BAM 为空（`SORF_BAM_FILTER` 零读段会导致下游 ORF 工具报错，见 gotcha 20 同类问题） |
| 少数样本拉偏队列窗口（逐样本并集可能过宽） | 输出逐样本推荐值；并集外若某样本贡献 >20% 的读长 则告警 |
| `-stub-run` 无法端到端验证（gotcha 24） | 用 1–2 个真实小样本端到端验证，不用 stub-run 做最终确认 |

## 3. 实施步骤

1. **独立验证先行**（不碰流程代码）：对 mouse_GSE120762 的一个样本，手工跑
   stage1 过滤 + prefilter RiboseQC，确认 `P_sites_calcs` 真的列出**全部**读长
   （预期 20–76），并与 `frame_preference_by_length.sh` 的输出对齐。
   *这是整个计划的前提假设，必须先证实。*
2. 实现 `bin/recommend_sorf_readlength.py` + 单元测试（喂合成的 `P_sites_calcs`）。
3. 实现 nf 模块与 wiring（D1–D6 按第 2.3 节建议取值）。
4. schema / config / MultiQC / CLAUDE.md 同步。
5. 哈希稳定性实测（第 2.5 节）。
6. 小规模端到端：2 样本、`--run_readlen_recommendation true`，核对
   (a) 推荐表数值、(b) `--sorf_read_len_auto false` 时输出与改动前逐字节一致。
7. `nf-core pipelines lint` + `nf-test`；分支 `dev` 提交并在 CHANGELOG.md 记录。

## 4. 回滚

- 默认参数下行为与现状一致，无需回滚路径。
- 若第 5 步实测发现哈希不稳定：把新增输入槽方案改为"独立进程分支"
  （新 `SORF_BAM_FILTER_LEN` 只在 `run_readlen_recommendation=true` 时走），
  彻底隔离新旧路径，接受两条路径的代码重复。

## 5. 不在本次范围

- 用推荐窗口**重跑**已完成的 GSE120762 / PRJEB26593（涉及数天机时，需单独决策）
- `ribowaltz_read_lengths` 的自动联动：本期只做**告警**（推荐窗口超出该参数时在
  MultiQC 提示），自动改写在下一期，因为它同样影响 hash 与 P-site 回退链
