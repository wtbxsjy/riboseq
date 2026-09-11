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

### 2.3 决策点

| # | 决策 | 结论 |
|---|---|---|
| D1 | 是否复用已有 `RIBOSEQC_PREFILTER` | **复用**（零重复实现） |
| D2 | 推荐的应用方式 | **默认仅报告**（`--sorf_read_len_auto` 打开回灌）。自动改 ORF 集会让已发表分析不可复现 |
| D3 | 窗口粒度 | **队列级**。逐样本窗口会让不同样本覆盖不同的读长集合，破坏跨样本可比性 |
| D4 | 合格判据 | **四层**（见 2.3.1），**不是**固定 50% |
| D5 | 窗口选取 | 候选**连续区间按读段质量**排序取最大者；带洞集合需读长列表参数（见 2.3.3） |
| D6 | 队列聚合 | **逐读长多数票**（≥半数样本），不是并集 |

#### 2.3.1 为什么不能写死 50%（2026-09-11 实测）

固定阈值是"选择题"的判据（挑一个读长做周期分析锚点），被挪用成"删除题"的判据
（决定丢掉哪些读段）就是范畴错误：**库差时合格集为空**，而硬编码 28-30 之所以存在
正是因为它保证非空。改用四层判据：

```
L1 效应量   E(L) = (p_max − 1/3) / (2/3)      0=随机 1=完美，跨库可比
L2 相对门   E(L) >= max(E_FLOOR, β · E_best)  随库自适应，不用改数字
L3 显著性   z = (p_max−1/3)/sqrt((1/3)(2/3)/n) >= 4.9   给小 n 兜底
L4 丰度地板 n(L)/Σn >= 1%
```

**非空回退阶梯**（消除"库差到没窗口"）：四层全空 → 取 E_best 单个读长 → 仍空 →
回退 `params.sorf_read_len_min/max` 并告警。硬编码从主逻辑降级为最后保险。

**再加一道库质量门 `E_BEST_MIN`**：队列层 `E_best` 低于 0.40 时判定"该库读长谱是平的，
不支持数据驱动选择"，**不发窗口**，回退默认值 + 高亮告警。理由见 2.3.2 —— 没有这道门，
噪声会被相对门捡进来。

#### 2.3.2 β 敏感性实测（8 样本，4 小鼠 + 4 人类）

跑法：`post_analysis/readlen_analysis/beta_sweep.py`。

| | 队列层 E_best | β=0.3 | β=0.4 | β=0.5 | β=0.6 | β=0.7 |
|---|---|---|---|---|---|---|
| 小鼠 | **0.884** | [27,28,29] | [28,29] | [28,29] | [28,29] | [28,29] |
| 人类 | **0.223** | — | — | — | — | — |

- **小鼠**：窗口在 β∈[0.4, 0.7] 上完全稳定 —— β 这个自由常数其实不敏感，
  这大幅降低了"常数拍脑袋"的风险。
- **人类**：E_best=0.223 < 0.40 → 判定不可用，回退默认 + 告警。该库整个 E 谱在
  0.03–0.28 之间平铺（这就是它 12 个样本全部没过 RiboseQC `max_coverage` 门的原因）。

**实测暴露并修掉的两个设计缺陷**（都会被数据触发，不是理论担忧）：

1. **连续区间不能按"长度个数"选**：人类 β=0.5 时选出 {26,28,29,30,32,33,34,35}，
   最长连续段是 [32..35]，算法会把窗口推荐成 32–35 而丢掉真正的 28-30。
   改为按**读段质量**排序候选区间。
2. **相对门在平谱上会捡噪声**：人类 E_best 仅 0.18–0.28，β·E_best 是极低门槛，
   噪声读长纷纷入选。故加 `E_BEST_MIN` 库质量门（D4 的补充）。

#### 2.3.3 带洞集合需要读长列表参数

小鼠在 β=0.5 时的真实答案是 **{25, 28, 29}** —— 25 nt 中位 E=0.576，四个样本一致
（0.592/0.515/0.775/0.561），是个真实信号，但与主峰隔离。而
`SORF_BAM_FILTER` 只能表达 `seqlen>=rlmin && seqlen<=rlmax` 一个连续区间。
取 [25,29] 会把该区间最差的 26（E=0.127）和 27（E=0.326）一起纳入。
→ 建议把参数升级为**读长列表** `--sorf_read_lengths 25,28,29,30`（min/max 保留兼容）。
这同时解决了 25 nt 的归属问题。

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

1. ~~**独立验证先行**（不碰流程代码）~~ **✅ 已完成 2026-09-11，8 样本（4 小鼠 + 4 人类）**
   - **前提证实**：RiboseQC 跑在 stage-1 BAM（唯一比对+去 contig、未切读长）上，
     `P_sites_calcs` 列出 **12–20 个读长**（流程内那版固定 3 行）。
   - **交叉验证**：与独立实现（`frame_preference_by_length.sh`）的逐读长排序
     **中位相关 0.905**（0.867–0.991），top5 重合多为 4/5。
   - 副产物：发现被窗口整体丢弃的真实信号 **25 nt**（4/4 小鼠样本两法 top-3 内）。
   - 详见 `READLENGTH_THRESHOLD_ANALYSIS.md` §7–8。

   **踩坑记录（实施时不要重犯）**：
   - RiboseQC 会**按 basename 相对路径**在后期重新打开基因组 FASTA
     （`Building aggregate P-sites profiles` → `scanFaIndex` → `.io_check_exists`），
     独立运行时必须把 FASTA 及其 `.fai` **链进工作目录**，否则跑到 ~13 分钟才失败。
   - RiboseQC 启动时会**删掉同目录下其它 `tmp_RiboseQC_*`**：并发跑多个样本而共用一个
     cwd 时，后启动的会删掉先启动者的 scratch 目录，后者在
     `Calculating P-sites positions and junctions` 阶段以
     `gzfile: cannot open .../P_sites_stats` 失败。**必须每样本一个工作目录**。
     注意 `P_sites_calcs` 是在该失败点**之前**写出的，所以崩掉的运行仍有可用的 calcs
     —— 但不要依赖这一点。
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
