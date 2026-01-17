# README.md 更新摘要

## 更新日期
2026-01-17

## 更新目的
根据项目最新实现，将 README.md 中的"设计规范"和"计划功能"更新为"已实现"状态，反映实际代码中已完成的功能。

## 主要更新内容

### 1. QC 和 sORF 预测过滤策略（已实现）

**更新位置**: "Ribo-seq Quality Control" 章节

**更改内容**:
- 标题从 "Design specification (development note): QC and filtering policy" 改为 **"QC and sORF prediction filtering policy (IMPLEMENTED)"**
- 将所有"将会"(will)的描述改为现在时态，表示已实现
- 明确指出所有 sORF 预测工具使用过滤后的 BAM 文件
- 说明未过滤的 BAM 保留用于基线 QC 和未来的定量模块

**关键实现点**:
```
- RiboseQC 仅应用于 type=riboseq 的样本
- 管道对 type=riboseq 样本运行两次 QC：
  - 预过滤 QC：在未过滤的基因组 BAM 上（用于基线 QC）
  - 后过滤 QC：在过滤后的基因组 BAM 上（评估过滤效果）
- 所有 sORF 预测工具使用过滤后的 BAM
```

### 2. 单样本 ORF 预测策略（已实现）

**更新位置**: "ORF Prediction Tools" 章节

**更改内容**:
- 标题从 "Design specification (development note): per-sample prediction only" 改为 **"Per-sample ORF prediction policy (IMPLEMENTED)"**
- 明确所有工具默认使用单样本模式
- 说明合并预测模式默认禁用，但可通过参数启用
- 列出具体的 ORF 预测工具名称

**关键实现点**:
```
- 所有工具（Ribo-TISH, Ribotricer, RiboCode, rp-bp, ORFquant）默认单样本模式
- 合并预测默认禁用，可通过 --sorf_predict_pooled true 启用（仅 Ribo-TISH）
- 管道保留足够的中间输出以便后续合并
```

### 3. BAM 输入模式过滤行为（已实现）

**更新位置**: "Starting from BAM Files" 章节

**更改内容**:
- 标题从 "Planned behaviour (sORF filtering + QC before/after)" 改为 **"sORF filtering + QC before/after behavior"**
- 将所有"将会"改为现在时态
- 明确提到 `SORF_BAM_FILTER` 模块
- 添加过滤控制参数说明

**关键实现点**:
```
- BAM 会被坐标排序和索引
- 对 riboseq 样本在未过滤 BAM 上运行 RiboseQC 和 Ribo-TISH Quality
- 通过 SORF_BAM_FILTER 模块创建过滤后的 BAM
- 在过滤后的 BAM 上再次运行 QC
- 可通过 --sorf_filter 参数控制（默认: true）
```

### 4. ORF 预测工具选择（更新参数说明）

**更新位置**: "Selecting ORF Prediction Tools" 章节

**更改内容**:
- 将 `--run_ribocode` 和 `--run_rpbp` 改为使用 `--skip_*` 参数
- 明确 RiboCode 和 rp-bp 默认被跳过
- 添加关于 BAM 输入模式下 RiboCode 不可用的说明

**新参数语法**:
```bash
# 启用可选工具
--skip_ribocode false  # RiboCode（需要转录组 BAM，BAM 输入模式不可用）
--skip_rpbp false --contaminant_fasta /path/to/contaminants.fa  # rp-bp

# 跳过默认工具
--skip_ribotish
--skip_ribotricer
--skip_riboseqc
--skip_orfquant
```

### 5. sORF BAM 过滤（已实现）

**更新位置**: "sORF BAM filtering" 章节

**更改内容**:
- 标题从 "sORF BAM filtering (design specification)" 改为 **"sORF BAM Filtering (IMPLEMENTED)"**
- 将"计划参数"改为"关键参数"
- 添加具体的默认值和实现细节
- 重新组织参数说明为更清晰的格式

**关键实现点**:
```bash
# 过滤规则
1. 唯一比对读段（可配置策略：auto|nh|mapq）
2. 排除不需要的染色体/序列（线粒体、叶绿体、模糊序列）
3. 读段长度过滤（默认 28-30 nt）

# 核心参数
--sorf_filter true                    # 启用/禁用过滤（默认: true）
--sorf_unique_mode auto               # 唯一性模式（默认: auto）
--sorf_unique_mapq 60                 # MAPQ 阈值（默认: 60）
--sorf_read_len_min 28                # 最小读长（默认: 28）
--sorf_read_len_max 30                # 最大读长（默认: 30）
--sorf_exclude_contigs_regex '...'    # 排除序列的正则表达式
--sorf_predict_pooled false           # 启用合并预测（默认: false）
```

### 6. ORFquant 自定义容器说明（新增）

**更新位置**: "Selecting ORF Prediction Tools" 章节之后

**新增内容**:
添加了一个重要提示框，说明 ORFquant 的补丁版本：

```
ORFquant Custom Container:
- 问题：BiocGenerics 导出的 Position 和 combine 与 ggplot2/gridExtra 冲突
- 解决方案：修改 ORFquant NAMESPACE 使用选择性导入
- 容器：提供预构建的补丁容器或从 Singularity.orfquant.patched.def 构建
- 使用：--orfquant_container /path/to/orfquant_patched.sif
- 或提供补丁包：--orfquant_pkg /path/to/ORFquant-1.1.tar.gz
```

## 代码实现验证

所有更新都基于实际代码实现的验证：

1. **SORF_BAM_FILTER 模块**: `modules/local/sorf_bam_filter/main.nf`
   - 实现了唯一性过滤、序列排除和读长过滤
   - 支持 auto/nh/mapq 三种唯一性模式

2. **RiboseQC 双重运行**: `workflows/riboseq/main.nf` 第 448-466 行
   - RIBOSEQC_PREFILTER: 在未过滤 BAM 上运行
   - RIBOSEQC_POSTFILTER: 在过滤后 BAM 上运行（meta.id 添加后缀避免冲突）

3. **ORFquant 依赖**: `workflows/riboseq/main.nf` 第 473-486 行
   - 检查 RiboseQC 是否跳过，自动跳过 ORFquant
   - 使用 RiboseQC 的 orfquant 输出

4. **参数定义**: `nextflow.config` 第 101-121 行
   - 所有 sorf_* 参数已定义并设置默认值
   - skip_riboseqc, skip_orfquant, skip_ribocode, skip_rpbp 参数已定义

5. **ORFquant 补丁**: `patched_packages/ORFquant-1.1/NAMESPACE`
   - 修改了 import 顺序，BiocGenerics 最后导入以确保优先级

## 文档改进

除了反映实现状态，还进行了以下改进：

1. **时态一致性**: 所有已实现功能使用现在时态
2. **参数清晰化**: 提供具体的参数示例和默认值
3. **实现细节**: 添加模块名称和具体行为说明
4. **用户指导**: 更清晰的使用说明和配置示例
5. **标注状态**: 使用 "(IMPLEMENTED)" 标记已实现功能

## 验证建议

建议用户验证以下内容：

1. 运行测试管道确认过滤功能工作正常
2. 检查输出目录中的 pre-filter 和 post-filter QC 结果
3. 验证 sORF 预测工具使用的是过滤后的 BAM
4. 测试自定义过滤参数（如调整读长范围）
5. 如使用 ORFquant，确认使用了补丁版本容器

## 相关文件

- `README.md`: 主要文档文件（已更新）
- `CLAUDE.md`: 开发者指南（已创建，包含架构说明）
- `nextflow.config`: 参数配置
- `workflows/riboseq/main.nf`: 主工作流实现
- `modules/local/sorf_bam_filter/main.nf`: BAM 过滤模块
- `containers/Singularity.orfquant.patched.def`: ORFquant 补丁容器定义

## 统计信息

- 修改行数: 127 行（71 行新增，56 行删除）
- 主要章节更新: 6 个
- 新增说明框: 1 个（ORFquant 容器）
- 参数更新: ~15 个参数的说明更新或新增
