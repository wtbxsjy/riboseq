# ORF QC Module — 使用文档

## 当前状态

ORF QC 模块已实现，但由于 Nextflow `-resume` 时的 channel 毒丸（PoisonPill）问题，**暂时从 workflow DAG 中移除**。目前通过独立脚本在 pipeline 完成后运行。

> **TODO**: 修复 `-resume` 时新进程的 channel 解析问题，重新集成到 workflow。

---

## 快速开始

```bash
cd /media/home/renzhe/riboseq/riboseq

# 设置路径
RESULT_DIR=/path/to/your/result
OUTPUT_DIR=${RESULT_DIR}/orf_qc
mkdir -p ${OUTPUT_DIR}

# 一键运行
python3 bin/extract_orf_qc_metrics.py \
    --results-dir ${RESULT_DIR} \
    --output ${OUTPUT_DIR}/tool_data.json

python3 bin/harmonize_orf_qc.py \
    --input ${OUTPUT_DIR}/tool_data.json \
    --output-prefix ${OUTPUT_DIR}/qc

python3 bin/compare_orf_tools.py \
    --tool-data ${OUTPUT_DIR}/tool_data.json \
    --unified-meta ${RESULT_DIR}/orf_unification/unified_orfs.metadata.tsv \
    --periodicity ${OUTPUT_DIR}/qc_periodicity.json \
    --output-prefix ${OUTPUT_DIR}/qc

python3 bin/generate_orf_qc_report.py \
    --confidence ${OUTPUT_DIR}/qc_orf_confidence.tsv \
    --periodicity ${OUTPUT_DIR}/qc_periodicity.json \
    --flags ${OUTPUT_DIR}/qc_sample_flags.json \
    --sample-id "MyProject" \
    --output ${OUTPUT_DIR}/qc_report.html
```

---

## 四个阶段

### Phase 1: `extract_orf_qc_metrics.py` — 指标提取

自动扫描结果目录，检测哪些工具产生了输出，解析所有格式。

**输入**: 结果目录（`--results-dir`）或文件列表（`--file-list`）

**支持的 8 个工具**:

| 工具 | 识别文件 | 提取内容 |
|------|---------|---------|
| RiboCode | `*_collapsed.txt` | per-ORF p 值、frame 覆盖率、RPKM |
| RiboseQC | `*_P_sites_calcs` | P-site offset、frame preference |
| riboWaltz | `*_psite_offset.tsv` | 两步校正 offset、total_percentage |
| Ribotricer | `*_translating_ORFs.tsv` | phase_score、valid_codons_ratio |
| Ribo-TISH | `*_pred.txt`, `*.para.py` | Fisher p-value、Frame Q-value |
| PRICE | `*.orfs.tsv` | GEDI ORF 检测结果、p-value |
| rp-bp | `*bayes-factors.bed.gz` | Bayes factor |
| ORFquant | `*_Detected_ORFs.gtf.gz` | ORF 分类、P-site 统计 |

**输出**: `tool_data.json` — 结构化 JSON

### Phase 2: `harmonize_orf_qc.py` — 指标对齐

**P-site offset harmonization**: 跨工具交叉验证，以 riboWaltz 为权威参考：
```
read_length | ribowaltz | riboseqc | ribocode | consensus | delta | flag
    28      |    12     |    12    |    12    |    12     |   0   |  OK
    29      |    12     |    12    |    12    |    12     |   0   |  OK
```

**Periodicity assessment**: 归一化各工具的三周期性指标到 [0,1]：
```
Ribotricer: 0.849
RiboseQC:   0.732
RiboCode:   0.633
Aggregate:  0.738 (GOOD)
```

**Sample flags**: 自动检测质量问题：
- `NO_PERIODICITY`: 无 read length 显示三周期性
- `LOW_PERIODICITY`: <50% read length 有三周期性
- `P_SITE_DISCORDANCE`: P-site offset 跨工具不一致
- `TOOL_FAILURE`: 工具运行失败或未产生输出

**输出**:
- `qc_psite_harmonized.tsv` — 跨工具 P-site 共识表
- `qc_periodicity.json` — 周期性和 flags
- `qc_sample_flags.json` — 机器可读质量标记

### Phase 3+4: `compare_orf_tools.py` — 跨工具比较 + OCS 评分

**跨工具比较**: 成对 Jaccard 重叠系数
```
RiboCode ↔ Ribo-TISH: Jaccard 0.580 (102,998 overlaps)
PRICE    ↔ Ribo-TISH: Jaccard 0.052 (13,101 overlaps)
```

**ORF Confidence Score (OCS)**: 每个 unified ORF 的 0-1 置信度

| 分量 | 权重 | 数据来源 |
|------|------|---------|
| S_translation | 0.30 | unified metadata `tool_scores` 列 |
| S_agreement  | 0.30 | unified metadata `tools` 列（跨工具支持数） |
| S_coverage   | 0.20 | unified metadata `unique_psites / length_aa` |
| S_periodicity| 0.15 | unified metadata `pN` 列 |
| S_readlevel  | 0.05 | Phase 2 全局周期评估 |

**置信度分级**:
| 等级 | OCS 范围 | 含义 |
|------|---------|------|
| High | ≥0.7 | 多工具支持、强周期性、高覆盖率 |
| Medium | 0.4-0.7 | 中等证据、1-2 工具支持 |
| Low | 0.2-0.4 | 弱证据、单工具、边界显著性 |
| Uncertain | <0.2 | 低于阈值，可能是假阳性 |

**输出**:
- `qc_orf_confidence.tsv` — 每个 ORF 的 OCS + 分量 + 等级
- `qc_tool_agreement.tsv` — 成对 Jaccard 矩阵

### Phase 5: `generate_orf_qc_report.py` — HTML 报告

交互式 Plotly 报告，包含 5 个标签页：

1. **Dashboard**: 关键指标汇总 + quality flags
2. **Read-Level QC**: P-site offset harmonization 表 + 周期性热图
3. **ORF Confidence**: OCS 分布 + 分类分布 + 长度分布
4. **Cross-Tool**: UpSet 图 + Jaccard 热图
5. **Detail Table**: 可搜索、可过滤的 per-ORF 详情表

---

## 输出文件一览

```
${OUTPUT_DIR}/
├── tool_data.json              # Phase 1: 结构化工具数据
├── qc_psite_harmonized.tsv     # Phase 2: P-site 共识表
├── qc_periodicity.json         # Phase 2: 周期性评估
├── qc_sample_flags.json        # Phase 2: 质量标记
├── qc_orf_confidence.tsv       # Phase 3+4: ORF 置信度
├── qc_tool_agreement.tsv       # Phase 3+4: 工具一致性
└── qc_report.html              # Phase 5: 交互式报告
```

---

## 故障排查

### 工具状态为 NOT_FOUND
确认 pipeline 已运行该工具，且结果文件在 `${RESULT_DIR}/orf_predictions/<tool>/` 下。

### OCS 全为 Uncertain (0.1-0.2)
检查 unified metadata 的 `tool_scores`、`pN`、`unique_psites` 列是否有值。若统一 ORF 结果来自旧版 pipeline，可能缺少这些列。重新运行 `UNIFY_ORF_PREDICTIONS` 即可。

### HTML 报告空白
确认 `orf_confidence.tsv` 包含数据。可以用 `head` 检查前几行。

### tool_data.json 过大
对于大规模数据（如人类基因组），`tool_data.json` 可能超过 1GB。考虑：
- 使用 `--file-list` 指定单个样本
- 分样本运行后合并结果

---

## 已知限制

1. **Workflow 集成**: 当前仅在脚本模式可用，不在 Nextflow DAG 中
2. **-resume 兼容**: 新进程加入 DAG 后 `-resume` 会触发 channel PoisonPill
3. **跨工具 ORF 匹配**: 当前用统一 ORF metadata 的 `tools` 列（仅覆盖 unification 成功解析的工具）
4. **内存**: 大型项目（>50 万 ORFs）需要 >4GB RAM
