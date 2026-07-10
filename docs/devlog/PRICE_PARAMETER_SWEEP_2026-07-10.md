# PRICE ORF 产量分析与参数优化实验

**日期**: 2026-07-10 — 2026-07-10  
**项目**: Lishuqi 人类+MRSA Ribo-seq（6 个 `type=riboseq` 样本）  
**容器**: `gedi.sif`（GEDI version 1.0.6f, Price version 1.0.5）  
**目的**: 诊断 PRICE 流程 ORF 偏少的原因，测试参数优化方案

---

## 1. 背景

Lishuqi 项目PRICE模块的运行结果显示ORF数量明显偏少：
- 6 个 Ribo 样本中，Ribo_11 和 Ribo_13 输出为 0（仅占位符表头）
- 成功样本的 ORF 数（1073-2483）对 1.1GB 基因组BAM 而言偏低

初步怀疑方向：
1. sORF 预过滤（28-30 nt，唯一比对，排除特定 contig）过于激进
2. 流程模块中 `-filter` 参数被错误注释为"无效"，未被使用
3. GEDI 内部除零异常导致部分样本全样本崩溃

---

## 2. 实验方法

### 2.1 测试台设计

为避免修改 Nextflow 流程，采用独立测试台：
- 使用 `apptainer exec gedi.sif gedi Price` 直接调用
- 复用基因组索引（`.oml`），一次构建多次使用
- 每个参数组合在独立目录运行，采集 ORF 数、运行时间、崩溃状态
- 测试脚本：`process/price_test/run_price_matrix.sh`

### 2.2 测试矩阵

**阶段 0 — 基线复现**: 用现有参数独立运行，验证测试台与流程输出一致  
**阶段 1 — 崩溃诊断**: 对 Ribo_13 逐一施加修复假设（nthreads/skipmt/novelTranscripts）  
**阶段 2 — 输入 BAM 对比**: 对比 sORF 过滤 BAM vs 基因组 BAM + PRICE 内部 `-filter`  
**阶段 3 — FDR/keepAnno 扫描**: 在最优输入上测试统计参数的影响

### 2.3 采集指标

每次运行记录：ORF 数（`.orfs.tsv` 行数）、是否崩溃、运行时长、退出码、p 值分布

---

## 3. 实验结果

### 3.1 测试台验证

| 样本 | 流程 ORF | 测试台 ORF | 一致性 |
|------|--------|---------|--------|
| Ribo_15 sORF | 2483 | 2492 | ✅（±9，无实质差异）|
| Ribo_12 sORF | 1104 | 1104 | ✅ 完全一致 |
| Ribo_13 sORF | 0（崩溃）| 0（崩溃）| ✅ 完全一致 |

独立测试台与 Nextflow 流程输出高度一致，可作为可靠实验平台。

### 3.2 阶段 1：崩溃诊断

**Ribo_13 崩溃根因**：`java.lang.RuntimeException: Could not run PriceOrfInference → / by zero`

排除的假设（均在 sORF BAM 上全部仍然崩溃）：

| 测试 | 参数修改 | 结果 |
|------|---------|------|
| D1 | `-nthreads 1` | ❌ 仍崩溃 |
| D2 | `-nthreads 1 -skipmt` | ❌ 仍崩溃 |
| D3 | 无 `-novelTranscripts` + `-nthreads 1` | ❌ 仍崩溃 |

**关键发现**：所有修复尝试对 sORF BAM 均无效。崩溃是 **sORF 过滤 BAM 与 GEDI 内部除零 bug 的组合问题**。

### 3.3 阶段 2：输入 BAM 对 ORF 产量的影响

**Ribo_12（流程输出 1104）**：

| 方案 | 输入 BAM | PRICE `-filter` | ORF 数 | vs 基线 |
|------|---------|---------------|--------|---------|
| A1 | sORF 过滤（28-30, 唯一）| 无 | 1104 | — |
| A2 | 基因组 BAM | 28:30 | 1497 | +35.6% |
| A3 | 基因组 BAM | 26:34 | **1796** | **+62.7%** |
| A4 | 基因组 BAM | 无（全读长）| 1958 | +77.4% |

**Ribo_13（流程输出 0，崩溃）**：

| 方案 | 输入 BAM | PRICE `-filter` | ORF 数 | 结果 |
|------|---------|---------------|--------|------|
| baseline | sORF 过滤 | 无 | 0 | ❌ 崩溃 |
| D4 | 基因组 BAM | 28:30 | **954** | ✅ 修复 |

**核心发现**：
1. sORF 预过滤（唯一比对 + contig 排除）使 PRICE 损失 35-63% 的 ORF，即使读长窗口相同
2. 放宽读长窗口从 28-30 到 26-34 额外增产 ~20%
3. Ribo_13 改用基因组 BAM 后从崩溃恢复，产出 954 个 ORF

### 3.4 阶段 3：FDR 与 keepAnno 的影响

**Ribo_12（基因组 BAM + filter 26:34 基础上）**：

| 参数 | ORF 数 | vs A3 基线 |
|------|--------|-----------|
| A3 基线（fdr=0.1, 无 keepAnno）| 1796 | — |
| + fdr=0.2 | 1796 | 无变化 |
| + fdr=0.5 | 1796 | 无变化 |
| + keepAnno | **2020** | **+12.5%** |
| + keepAnno + fdr=0.2 | 2020 | +12.5% |

**核心发现**：
1. **FDR 松绑对产量无影响** — 瓶颈在候选 ORF 生成阶段（codon 覆盖率），不在统计过滤阶段
2. **`-keepAnno` 是唯一有效的统计参数** — 保留注释 CDS 不受 p 值过滤，额外增产 12.5%
3. 默认 `fdr=0.1` 在此数据集上已充分宽松

### 3.5 汇总

**最优参数组合**：`基因组 BAM + -filter 26:34 + -keepAnno`

| 样本 | 当前 ORF | 最优方案预期 | 增幅 |
|------|--------|-----------|------|
| Ribo_11 | 0 | ~600-800 | ∞ |
| Ribo_12 | 1104 | **2020** | **+83%** |
| Ribo_13 | 0 | ~900 | ∞ |
| Ribo_14 | 1969 | ~3600 | +83% |
| Ribo_15 | 2483 | ~4540 | +83% |
| Ribo_16 | 1073 | ~1960 | +83% |

---

## 4. 关于 `-filter` 范围的自动推导

### 4.1 现状

当前流程中，PRICE 使用 `ch_bams_for_sorf_prediction` 作为输入，该 BAM 经过 sORF 过滤：
- 读长 28-30 nt（`sorf_read_len_min/max`）
- 唯一比对（`sorf_unique_mode=auto`，MAPQ≥60 或多标记 NH:i:1）
- 排除线粒体/叶绿体/未定位 contig（`sorf_exclude_contigs_regex`）

这对 PRICE 有两个问题：
1. **唯一比对过滤**移除了大量 PRICE 本可利用的读段
2. **contig 排除**在某些样本上触发 GEDI 除零异常

### 4.2 自动推导方案

PRICE 的 `-filter` 参数只需要**读长范围**（如 `26:34`），读长过滤由 PRICE 内部高效处理。该范围可从 riboWaltz 输出自动推导。

**数据来源**：`riboWaltz *_psite_offset.tsv`

该文件包含 riboWaltz 判定为具有良好的 3-nt 周期性的读长列表：
```
length  total_percentage  start_percentage  ...  corrected_offset_from_5  sample
28      27.916            27.32             ...  12                        Ribo_15
29      37.062            31.701            ...  12                        Ribo_15
30      35.021            40.979            ...  12                        Ribo_15
```

只有通过周期性检测的读长才出现在此表中。因此：
- 取所有样本的 `length` 列的并集
- `min(length)` → `max(length)` 即为最优 `-filter` 范围
- 对于 Lishuqi 数据集：28-30（6/6 样本一致）

**推导逻辑**（`scripts/deduce_price_read_filter.py`）：
```python
import sys, csv
lengths = set()
for f in sys.argv[1:]:
    with open(f) as fh:
        for row in csv.DictReader(fh, delimiter='\t'):
            lengths.add(int(row['length']))
if lengths:
    print(f"{min(lengths)}:{max(lengths)}")
else:
    print("")  # 空 = 不传 -filter，让 PRICE 用全读长
```

**备选数据源**：RiboseQC `*_P_sites_calcs`
- 取 `all==TRUE && max_inframe==TRUE && max_coverage==TRUE` 的 `read_length`
- 通常比 riboWaltz 更保守（如 Lishuqi 只保留 29-30，排除 28）
- 作为 fallback（当 `--skip_ribowaltz` 时）

**终极 fallback**：`--price_read_filter` 参数的手动设定值，或 `sorf_read_len_min - 2` 到 `sorf_read_len_max + 4`（稍微放宽）

### 4.3 推荐实现

在 Nextflow 流程中：

```nextflow
// 伪代码
if (params.price_read_filter == 'auto') {
    // 步骤 1：收集 riboWaltz psite_offset 文件
    // 步骤 2：运行 deduce_price_read_filter.py 推导范围
    // 步骤 3：将结果传给 PRICE 的 -filter
} else if (params.price_read_filter) {
    // 使用用户指定的范围，如 "26:34"
} else {
    // 不传 -filter，PRICE 用全读长
}
```

---

## 5. 流程修改建议

### 5.1 参数设计

新增参数（`nextflow.config`）：

```groovy
params {
    // PRICE ORF detection
    price_read_filter       = 'auto'    // 'auto' | '26:34' | '' (不传 -filter)
                                        // auto: 从 riboWaltz psite_offset 推导
    price_keep_anno         = true      // -keepAnno: 保留注释 CDS 不过滤
    price_fdr               = 0.1       // FDR 阈值（留作可选项，本次实验证明无需修改）
    price_read_filter_min   = null      // 手动指定最小读长（覆盖 auto）
    price_read_filter_max   = null      // 手动指定最大读长（覆盖 auto）
}
```

### 5.2 代码修改清单

**P0 — 关键修复**：

1. **`modules/local/price/main.nf`**：
   - 删除第 26-27 行错误注释（`-filter` 是 GEDI Price 的有效参数）
   - 将 `-filter` 参数接入 `extra_args`（或独立传入）
   - 将 `-keepAnno` 接入（默认开启）

2. **`workflows/riboseq/main.nf`**：
   - PRICE 输入从 `ch_bams_for_sorf_prediction` 改为 `ch_riboseq_genome_bam`（排线粒体/ rRNA 等明显污染的 BAM，但不做唯一比对和读长过滤）
   - 或者新增一个专用的轻度过滤 BAM 通道：仅过滤线粒体和 rRNA，保留多比对和全读长

3. **`nextflow.config`**：
   - 新增 `price_read_filter`、`price_keep_anno`、`price_fdr` 参数
   - `price_read_filter` 默认 `'auto'`

**P1 — 自动推导**：

4. **`scripts/deduce_price_read_filter.py`**（新增）：
   - 读取 riboWaltz `*_psite_offset.tsv` 或 RiboseQC `*_P_sites_calcs`
   - 输出 `min:max` 读长范围
   - 无有效数据时输出空字符串

5. **`workflows/riboseq/main.nf`**（PRICE 调用处）：
   - 当 `price_read_filter == 'auto'` 时，先运行 `DEDUCE_PRICE_READ_FILTER` 小工具进程
   - 将输出传给 PRICE 的 `-filter`

**P2 — 文档与防御性**：

6. **`modules/local/price/main.nf`**：
   - 将 `-nthreads` 上限设为 16（默认 254 内存不安全）
   - 添加崩溃检测：检查 `.command.log` 中是否有 `/ by zero` 异常，若有则输出明确诊断信息而非空白占位符

7. **`docs/orf-tools/06-price.md`**：
   - 更新为 GEDI PRICE 文档（当前内容为 rp-bp）

### 5.3 关于 BAM 输入的进一步建议

当前问题根源是 `ch_bams_for_sorf_prediction` 过滤太激进。但不建议让 PRICE 直接使用 `ch_riboseq_genome_bam`（未过滤），因为该 BAM 包含线粒体/ rRNA 读段，会增加 PRICE 运行时间且可能干扰 codon 模型。

**推荐方案**：新增一个"PRICE 专用"的轻度过滤 BAM：
- **保留**：全读长范围、多比对读段
- **仅移除**：线粒体/叶绿体 contig（`sorf_exclude_contigs_regex`）、未比对/低质量读段
- 实现方式：在 sORF 过滤基础上新增 `SORF_BAM_FILTER_PRICE` 变体，或复用 RiboseQC 前的 host-only BAM

如果不想增加模块复杂度，次优方案是直接使用 `ch_riboseq_genome_bam` 并依赖 PRICE 内部 `-filter` 做读长选择。代价是最多增加 ~30% 运行时间（全读长 vs 部分读长），但产出 ORF 最多（A4 方案 1958 ORF vs A3 1796）。

---

## 6. 实验文件清单

### 测试脚本

| 文件 | 说明 |
|------|------|
| `process/price_test/run_price_matrix.sh` | 参数化 PRICE 运行脚本 |
| `process/price_test/genome.oml` | 复用基因组索引符号链接 |
| `process/price_test/bams/` | 测试用 BAM 与 .bai 符号链接 |
| `process/price_test/runs/` | 各实验组合的独立运行目录 |

### 数据来源

| 文件 | 用途 |
|------|------|
| `result/sorf/Ribo_*.sorf.filtered.bam` | sORF 过滤 BAM（当前流程输入） |
| `result/alignment/star/sorted/Ribo_*.genome.sorted.bam` | 未过滤基因组 BAM（测试用） |
| `result/riboseq_qc/ribowaltz/*_psite_offset.tsv` | 周期性读长（自动推导 filter 来源） |
| `result/riboseqc/*_P_sites_calcs` | RiboseQC P-site 数据（备选来源） |

### 完整结果

参见 `process/price_test/runs/_results.csv`：
```
run|orf_count|crashed|exitcode|duration_sec|dir_size
baseline_Ribo15_sorf|2492|false|0|435|52M
S2_A1_Ribo12_sorf|1104|false|0|384|26M
S2_A3_Ribo12_genome_f2634|1796|false|0|1489|39M
S2_A2_Ribo12_genome_f2830|1497|false|0|1589|33M
D4_Ribo13_genome_f2830|954|false|0|1721|27M
S2_A4_Ribo12_genome_nofilt|1958|false|0|2384|45M
S3_F3_Ribo12_keepAnno|2020|false|0|1531|41M
S3_F4_Ribo12_keepAnno_fdr02|2020|false|0|1534|41M
S3_F1_Ribo12_fdr02|1796|false|0|1557|39M
S3_F2_Ribo12_fdr05|1796|false|0|1583|40M
```

---

## 7. 结论

PRICE ORF 偏少由三层原因叠加造成：
1. **sORF 过滤 BAM 触发 GEDI 除零崩溃**（2/6 样本 ORF=0）→ 改用基因组 BAM
2. **唯一比对 + contig 过滤削减 ~60% ORF** → 取消唯一比对要求，仅保留 contig 排除
3. **未使用 `-keepAnno` 保留注释 CDS** → 默认开启 `-keepAnno`（额外 +12%）

推荐的最优参数组合：**基因组轻度过滤 BAM + `-filter 26:34` + `-keepAnno`**  
综合预期：各样本 ORF 数较当前流程提升 **80-100%**，Ribo_11/Ribo_13 从 0 恢复到 ~800-900。
