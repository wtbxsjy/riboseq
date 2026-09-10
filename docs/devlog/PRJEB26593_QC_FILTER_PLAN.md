# PRJEB26593 Ribo-seq 表达量质量控制方案（对照 rice/maize post_analysis 过滤流程）

- 日期: 2026-09-07
- 状态: 已确认（2026-09-07），执行中 —— Step 3 P-site 纯度已后台启动；Step 1/2 + pN 分布先行
- 目标: 对 PRJEB26593 统一 ORF 集（1,215,479 个）套用 rice/maize post_analysis 的过滤级联
  （psite reads > 9、psite 比例 ≥ 50%、psite 位置 > 2 nt、长度 ≤ 50 aa、top 10% per biotype），
  产出高置信 ORF 列表 + 过滤计数报告。

## 1. 输入文件核对（已验证存在）

| 文件 | 路径 | 说明 |
|---|---|---|
| `unified_orfs.bed.gz` (25MB) | `run/human_PRJEB26593/result/orf_unification/` | BED12, 1,215,479 ORF |
| `unified_orfs_expression_summary.tsv` (230MB) | 同上 | 列: orf_id, chrom, start, end, strand, 12×`{sample}_reads`, 12×`{sample}_pN`, total_reads, n_expressed_samples |
| `unified_orfs_expression_rpkm_tpm.tsv` (388MB) | 同上 | RPKM/TPM（下游翻译效率用，过滤本身不用） |
| `unified_orfs.metadata.tsv` (3.3GB) | 同上 | 全注释 |
| riboseqc bedgraphs (5.5GB) | `run/human_PRJEB26593/result/riboseqc/` | 12 样本 × `{sample}_P_sites_{plus,minus}.bedgraph` + `{sample}_coverage_{plus,minus}.bedgraph`（非 uniq 版，27–46MB/文件，ERR2603026 最小 8.6MB） |
| `gencode_results.orfs.out.gz` | `result/orf_classification/gencode/` | **473,069 行**（只覆盖统一集的 ~39%），列含 `orf_biotype`；统一 ORF id 在 `all_orf_names`/`phaseI_id` 列 |
| `orftype_classification.tsv` (3.3GB) | `result/orf_classification/orf_type/` | **1,215,479 行 = 全覆盖**。列: orf_id(统一 id), length_aa, tools, samples, total_reads, total_psites, pN, is_cds_overlap, overlapping_genes, orf_type_category |
| ORFquant 分类 | — | result/ 下未找到 `orfquant_classification.tsv`，本 QC 不需要它 |

## 2. 过滤级联（与 rice slide-rice-filter.js 对齐）

rice 级联（参考基准）: 394,053 → 排除 CDS 285,921 → Stage1 263,749 → Stage2 250,190 → ≤50aa 102,582 → top10% 10,255。

PRJEB26593 对应步骤:

### Step 1 — CDS 排除（分类注释合并）
- 主依据 `orftype_classification.tsv`（全覆盖）: 排除 `is_cds_overlap == 1`（canonical CDS 重叠）。
- 补充 `gencode_results.orfs.out.gz`: `orf_biotype == 'CDS'` 的 ORF 也排除；该文件同时提供
  `dORF/uORF/doORF/uoORF/intORF/lncRNA` 等细分 biotype，供 Step 5 分组使用。
- 未在 gencode 结果中的 ORF 按 orf_type_category（novel/uORF/dORF/…）归类，视为 novel 保留。
- 输出: `prjeb26593_step1_nonCDS.tsv` + 排除计数。

### Step 2 — Stage 1 表达量预过滤（reads + pN）
- 输入 `unified_orfs_expression_summary.tsv`。
- 条件: **存在 ≥1 个样本满足 `{sample}_reads > 9 AND {sample}_pN > 0.5`**（≥1 样本即保留，
  时间序列 0/2/5/10h 下保留条件性翻译 ORF）。
- pN 语义（unify 脚本）: `max_bg_value × n_intervals / sum_bg_values`，峰度×覆盖的周期性代理指标。
- **阈值决策（用户确认 2026-09-07）**: 不直接沿用 rice 的 0.5 —— 先绘制 `{sample}_pN` 分布
  （本数据集 12 个样本 frame_preference < 50%，弱周期性，pN 整体偏低），再依据分布与候选阈值计数表定档。
- 输出: `prjeb26593_stage1_passed.tsv` + `pn_distribution.png` + `stage1_threshold_table.tsv`。

### Step 3 — Stage 2 P-site 纯度过滤（compute_psite_purity.py 复用，已修复 1 处 bug）
- 命令: `python3 post_analysis/scripts/compute_psite_purity.py --bed unified_orfs.bed --riboseqc-dir result/riboseqc --output prjeb26593_psite_purity.tsv --workers 8`
- 机制（已验证源码）: 按 `{sample}_P_sites_plus.bedgraph` 发现样本 → 逐样本读 P_sites{plus,minus} +
  coverage{plus,minus} → 按 chrom bisect 索引与 1.2M ORF 求交 → 链特异位置（+ 从 orf_start，− 从 orf_end）。
- 输出列: 每样本 `{sample}_p_site_GSE/_reads_GSE/_p_site_pct/_p_site_pos/_not_p_site_GSE` +
  全局 `total_psites/total_reads/global_p_site_pct/global_p_site_pos/global_p_site_pos_sd/n_samples_with_psites`。
- **执行中发现并修复的两个问题（2026-09-08，详见 §8）**:
  1. 每样本 `p_site_pos` 列恒为 1.0（脚本用错 key：值加权和当位置加权和除）→ 已修复脚本
     （本运行未重跑 9.5h 任务；位置判据改用 ORF 级 `global_p_site_pos`，本就正确）。
  2. `p_site_pct` 量纲无效：RiboseQC coverage bedgraph 是 RPM 归一化的（全基因组积分 ≈1e6），
     而 P_sites bedgraph 是原始计数 → 原始脚本的 pct 无意义。
- 过滤条件（修正后）: **存在 ≥1 个样本满足 `p_site_GSE > 9 AND pct ≥ 0.5`，
  其中 `pct = {sample}_p_site_GSE / {sample}_reads`（表达矩阵原始 reads，与 P-site 同量纲）**，
  且 ORF 级 `global_p_site_pos > 2`（加权平均 P-site 位置，链特异）。
  实测 pct 分布 q50=0.49（区间 0.35–0.625），区分力正常；pos>2 移除 0 个。
- 输出: `prjeb26593_stage2_passed.tsv`。

### Step 4 — ORF 长度 ≤ 50 aa
- 直接用 `orftype_classification.tsv` 的 `length_aa` 列（全覆盖，无需重算）。
- 保留 `length_aa ≤ 50`（rice 的 sORF 取向；如需全部长度可在执行时换向或出全集）。
- 输出: `prjeb26593_filtered_len50aa.tsv`。

### Step 5 — Top 10% per biotype（ggRibo 作图子集）
- 按 biotype 分组（gencode `orf_biotype` 优先，缺失用 `orf_type_category`），每组按
  `total_psites` 降序取前 10%（用户确认）。
- 供 ggRibo 可视化（yuanliang 规则 max_file=10 样本；PRJEB26593 12 样本，选 10 个或按时间点取代表）。
- P-site offset 已有: `result/riboseq_qc/psite_correction/{sample}_corrected_for_ORFquant` × 12。
- 输出: `prjeb26593_top10_per_biotype.tsv` + ggRibo 图。

## 3. 输出目录与脚本

- 输出目录: `/home/25119231r/riboseq/post_analysis/PRJEB26593/`（与 rice/maize 并列）。
- Step 3 直接复用 `post_analysis/scripts/compute_psite_purity.py`（不改动）。
- Step 1/2/4/5 写一个过滤脚本 `post_analysis/scripts/filter_orf_expression.py`（读上述 TSV，
  按级联出每步计数 + 通过列表），避免在对话里散落一次性命令。
- 最终报告: `prjeb26593_filter_report.md`（每步通过/排除计数表）。

## 4. 与 rice 的差异与风险点

1. **pN 阈值**: 本数据集弱周期性，p_site/reads 两套值整体偏低；建议严格 0.5 / 放宽 0.3 双档先看分布。
2. **gencode 覆盖不足**（473K/1.2M）: CDS 排除主依据切换到 orftype（全覆盖），gencode 仅作 biotype 细分注释。
3. **两套 reads 来源**: Stage1 的 `{sample}_reads` 来自 unify 内部流式 coverage；Stage2 的 `reads_GSE`
   来自 RiboseQC coverage bedgraph，两者数值略有差异（rice 同样两套，不影响逻辑）。
4. **orfquant_classification.tsv 不存在是预期行为**: rice/maize 里该文件来自手动跑的 ORFquant；
   本 fork `nextflow.config:139` 默认 `skip_orf_classify_orfquant = true`（待 data.table 重写），
   PRJEB26593 管线只跑了 GENCODE + ORF-type 两个分类器。本 QC 不需要该文件。
5. **资源协调**（gotcha 27）: 共享 188GB 机器，其他租户常态占 ~95GB。Step 3 用 `--workers 8`
   （每 worker 内存 ~200-300MB，主进程 ORF 索引 ~300MB）；GSE208041 当前已暂停，若其重启需错峰。

## 5. 预计耗时与资源

| 步骤 | 耗时（估） | 峰值内存 |
|---|---|---|
| Step 1 分类合并 | ~10 min（3.3GB TSV 流式 join） | ~2GB |
| Step 2 Stage1 | ~5 min（230MB awk/python） | ~1GB |
| Step 3 P-site 纯度 | **1–3 h**（48 个 bedgraph 共 ~1.7GB 流式扫描 × 1.2M ORF 索引） | ~4–6GB（8 workers） |
| Step 4 长度 | 秒级 | — |
| Step 5 top10% + ggRibo | 10 min + 作图数小时（按需） | ~2GB |

## 6. 交付物

1. `prjeb26593_filter_report.md` — 每步计数表（对照 rice 的 394K→10,255 格式）
2. 各步通过列表 TSV（step1/2/3/4/5）
3. `prjeb26593_psite_purity.tsv` — 12 样本 × 5 指标的完整纯度矩阵
4. ggRibo 图（top10% 子集）
5. 可选: PPT 新增 QC 过滤页

## 7. 已确认决策（2026-09-07）

1. pN 阈值: 先绘制分布、给出候选阈值计数表，再定档（不预先锁 0.5）。
2. Step 4: 保留长度过滤（length_aa ≤ 50）。
3. Step 5 排名依据: total_psites。
4. Step 2/3 的"≥1 样本"逻辑: 确认保留（时间序列）。
5. orfquant_classification.tsv: 不需要；缺失属预期（见风险点 4）。
6. PPT: 已同步更新 —— 修正 slide-02 分类描述（两路: GENCODE + ORF-type），新增 QC 过滤页（slide-04）。

## 8. 执行进度（2026-09-07）

- **Step 3 P-site 纯度**: 后台运行中（8 workers，各 ~99% CPU；12 样本、1,215,479 ORF 已索引完成）。
- **Step 1 CDS 排除: 完成**。1,215,479 → 排除 is_cds_overlap=1 的 1,015,770 + gencode CDS 补充 9,665
  → **保留 190,044**（15.6%）。gencode out 的 `all_orf_names` 列已确认即统一 orf_id
  （473,068 行，CDS 映射到 404,794 个统一 id）。
- **Stage 1: 完成，pN 经分布检验后被弃用**。本 fork 的 pN = `max_bg_value × n_intervals / sum_bg_values`，
  因 max ≥ sum/n_intervals，**pN ≥ 1 恒成立**（实测 10,125,599 条 reads>9 记录 min=1.0，<1 占比 0%；
  分布 q50=2.57, q90=6.87, max=76.6）。rice 的 pN>0.5 在此定义下恒真、无区分力。
  结论: Stage 1 有效门控 = `{sample}_reads > 9`（≥1 样本），通过 **65,261**（non-CDS 的 34.3%）；
  psite 比例 ≥50% 的要求由 Stage 2 的 p_site_pct ≥ 0.5 承担（与用户原始需求一致）。
- 产物: `step1_nonCDS.tsv.gz`、`prjeb26593_stage1_passed.tsv`（含 n_samples_reads_gt9）、
  `pn_distribution.png`、`stage1_threshold_table.tsv`、`step1_stats.json`，均位于
  `/home/25119231r/riboseq/post_analysis/PRJEB26593/`。
- PPT: slide-02 已修正、slide-04 新增（Step1/Stage1 数字已填入），已重新编译。

## 9. 执行进度（2026-09-08）— 级联全部完成

- **Step 3 纯度矩阵**: 完成（888,360 ORF 有指标，耗时 ~9.5h，2026-09-07 晚 ~ 09-08 早）。
- **Stage 2（修正定义）**: 完成。`psite > 9` 任一 ≥1 样本 59,086 → 加 `pct ≥ 0.5` 后 **50,214**
  （stage1 的 76.9%；pct = purity p_site_GSE / 表达矩阵 {sample}_reads，两者同为原始计数）；
  `global_p_site_pos > 2` 移除 0 个。
- **Step 4 长度 ≤ 50 aa**: 13,552（stage2 的 27.0%）。
- **Step 5 top 10% per biotype（按 total_psites 排名）**: **1,360** ORF / 11 biotypes
  （intORF 314, isoform 383, lncRNA 221, uORF 152, novel 86, novel_upstream 79, doORF 58,
  dORF 36, uoORF 26, novel_downstream 4, truncated 1）。
- **级联总表**: 1,215,479 → CDS 排除 190,044 → Stage1 65,261 → Stage2 50,214 → ≤50aa 13,552 → top10% **1,360**。
- 产物: `prjeb26593_stage2_passed.tsv`、`prjeb26593_filtered_len50aa.tsv`、
  `prjeb26593_top10_per_biotype.tsv`、`prjeb26593_filter_stats.json`（含各级联计数与 biotype 分布）、
  `prjeb26593_filter_report.md`，均位于 `/home/25119231r/riboseq/post_analysis/PRJEB26593/`。
- 脚本修复: `compute_psite_purity.py` 每样本 pos 用错 key（`_pos_wt`→`_pos_sum`）已修复并同步修了
  全局聚合块里一处遗留的 `pos_wt` 引用（NameError）；本次运行未重跑 9.5h 纯度任务（结论不依赖每样本 pos）。
- PPT: slide-04 已填入全部最终数字（50,214 / 13,552 / 1,360），已重新编译。
