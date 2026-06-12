# 项目分析摘要 - Nextflow 流程同步计划

**生成日期**：2026-01-30  
**项目**：nf-core/riboseq GENCODE 注释工作流集成  
**状态**：✅ 分析完成

---

## 📌 执行摘要

已完成对当前 Nextflow 流程和 `scripts/singularity_single_tool_tests` 中的新增功能的全面分析，制定了详细的同步计划。

### 核心发现

| 方面 | 当前状态 | 新增功能 | 集成工作量 |
|------|--------|--------|----------|
| **ORF 预测工具** | 6 个工具已集成 | 无新增工具 | - |
| **格式转换** | 无 | Ribo-TISH, Ribotricer → GENCODE | 🔴 P0 |
| **多工具统一** | 无 | GENCODE ORF Mapper | 🔴 P0 |
| **定量分析** | 无 | P-site 计数矩阵 | 🟡 P1 |
| **工具比较** | 无 | 高置信度 ORF 识别 | 🟡 P1 |

### 关键数字

- **新增模块**：5 个
- **新增子工作流**：3 个
- **预计工作量**：57 小时（~7-8 个工作日）
- **预计完成**：4-5 周
- **新增参数**：8+ 个

---

## 📁 生成的文档

本分析已生成以下三份详细文档，存储在项目根目录：

### 1. [SYNC_PLAN_20260130.md](./SYNC_PLAN_20260130.md) - 完整同步计划

**内容**（10 个章节）：
- ✅ 执行摘要
- ✅ 详细差异分析（Nextflow vs 单工具脚本）
- ✅ 3 个关键脚本（14-16）的功能详解
- ✅ Python 转换脚本位置和功能
- ✅ 同步计划（从脚本到 Nextflow）
  - 核心模块创建（5 个新模块的规范）
  - 子工作流整合（3 个新子工作流）
  - 主工作流集成（在 main.nf 中添加）
- ✅ 5 周路线图
- ✅ 文件变更清单
- ✅ 关键考虑事项
- ✅ 验收标准
- ✅ 快速参考（关键输入文件格式）

**用途**：项目管理、架构规划、高层决策

---

### 2. [IMPLEMENTATION_CHECKLIST.md](./IMPLEMENTATION_CHECKLIST.md) - 逐周实现检查清单

**内容**（每周详细任务）：
- ✅ 快速概览（优先级划分）
- ✅ 第 1 周：基础准备（9 小时）
  - 模块框架创建
  - Python 脚本准备
  - 镜像验证
  - 测试数据准备
- ✅ 第 2 周：核心模块实现（12 小时）
  - 3 个关键模块完整实现
  - 单元测试
- ✅ 第 3 周：子工作流集成（14 小时）
  - 3 个子工作流创建
  - 集成测试
- ✅ 第 4 周：主工作流集成（10 小时）
  - nextflow_schema.json 更新
  - 主工作流导入和调用
  - 完整流程测试
- ✅ 第 5 周：文档和优化（12 小时）
  - 文档更新
  - 性能优化
  - 发布准备
- ✅ 时间和资源投入总结
- ✅ 风险评估和缓解措施

**用途**：开发人员参考、任务分配、进度跟踪

---

### 3. [TECHNICAL_REFERENCE.md](./TECHNICAL_REFERENCE.md) - 技术参考指南

**内容**（10 个技术主题）：
- ✅ 工作流架构（数据流图、处理阶段）
- ✅ 数据格式规范
  - Ribo-TISH predict 输出格式
  - Ribotricer translating_ORFs 格式
  - GENCODE FASTA/BED 格式
  - 统一注释 GTF 和 OUT 格式
- ✅ 模块接口定义（Groovy 代码示例）
- ✅ 子工作流接口定义
- ✅ Nextflow 通道设计
- ✅ 参数配置（nextflow_schema.json, modules.config）
- ✅ 容器镜像和依赖列表
- ✅ 错误处理策略
- ✅ 性能调优建议
- ✅ 故障排查指南

**用途**：开发人员技术参考、集成测试、调试

---

## 🔍 关键分析结果

### 1. 新增功能概览

```
singularity_single_tool_tests 中的新增功能（与 Nextflow 对比）

脚本 14-16（ORF 注释流程）：已完全缺失，需要集成
├── 14_ribotish_to_gencode.sh      → CONVERT_RIBOTISH_TO_GENCODE 模块
├── 15_ribotricer_to_gencode.sh    → CONVERT_RIBOTRICER_TO_GENCODE 模块
└── 16_gencode_orf_mapper.sh       → GENCODE_ORF_MAPPER 模块

脚本 17-19（定量和比较分析）：待确认，部分需要创建
├── 17_unify_predictions.sh        → 子工作流中的合并逻辑
├── 18_classify_orfs.sh            → ANALYZE_TOOL_COMPARISON 模块
└── (R 分析脚本)                   → QUANTIFY_ORFS 模块
```

### 2. 核心数据流

```
ORF 预测输出（多工具）
    ↓
格式转换（14-15）
    ↓
文件合并（collectFile）
    ↓
GENCODE 统一注释（16）
    ↓
├─→ 定量分析（17）→ orf_counts_matrix.csv ⭐
├─→ 工具比较（18）→ high_confidence_orfs.tsv ⭐
└─→ 其他分析
```

### 3. 集成架构

```
workflows/riboseq/main.nf
    ├── 包含：GENCODE_ANNOTATION（子工作流）
    ├── 包含：QUANTIFY_ORFS（子工作流）
    ├── 包含：TOOL_COMPARISON（子工作流）
    └── 调用参数：
        - skip_gencode_annotation
        - ensembl_dir ✅ 关键参数
        - orf_collapse_threshold
        - skip_orf_quantification
        - skip_tool_comparison
```

### 4. 新增输出文件

| 文件 | 用途 | 优先级 |
|------|------|--------|
| `*.orfs.gtf` | 用于 featureCounts 定量 | 🔴 P0 |
| `*.orfs.out` | 工具检测矩阵，用于比较 | 🔴 P0 |
| `orf_counts_matrix.csv` | DESeq2 输入 | 🟡 P1 |
| `high_confidence_orfs.tsv` | 高置信度 ORF（≥2工具） | 🟡 P1 |
| `upset_plot.pdf` | 工具交集可视化 | 🟢 P2 |

---

## 📊 资源和依赖

### 外部依赖
- ✅ gencode-riboseqORFs 库（包含转换脚本和 ORF Mapper）
- ✅ Ensembl 注释数据（用户提前准备或脚本下载）
- ✅ R 包：GenomicRanges, rtracklayer, tidyverse, UpSetR

### 容器镜像
- `biopython:1.81`（格式转换，~200 MB）
- `gencode-orf-mapper:v1.1`（ORF 统一注释，~500 MB）
- `rocker/tidyverse:4.3`（R 分析，~2 GB）

### Python 脚本位置
```
scripts/gencode_converters/
├── ribotish_to_gencode.py
├── ribotricer_to_gencode.py
└── orfquant_to_gencode.py（可选）
```

---

## 🎯 建议的实施步骤

### 第 1 阶段：批准和准备（2-3 天）
- [ ] 审批本计划
- [ ] 确保 gencode-riboseqORFs 库可用
- [ ] 准备/验证 Ensembl 注释数据

### 第 2 阶段：实施（4-5 周）
按照 [IMPLEMENTATION_CHECKLIST.md](./IMPLEMENTATION_CHECKLIST.md) 逐周执行

### 第 3 阶段：验收（1 周）
- [ ] 完整流程测试
- [ ] 文档审查
- [ ] 发布准备

---

## 📝 关键建议

### 1. 优先级管理
- 🔴 **P0（第 1-2 周）**：脚本 14-16（格式转换和统一注释）
  - 这是下游分析的基础
  - 强烈推荐首先集成
- 🟡 **P1（第 3 周）**：脚本 17-19（定量和比较分析）
  - 高价值功能，但不阻塞其他工作

### 2. 参数策略
- 所有新功能应默认 **禁用**（`skip_gencode_annotation=true`）
- 用户需显式启用和提供 `ensembl_dir`
- 这样保证向后兼容性

### 3. 测试策略
- 在 test profile 中包含 GENCODE 注释功能
- 提供最小化的 Ensembl 测试数据集
- 验收标准：test 运行 < 5 分钟

### 4. 文档优先级
- 必须：README, usage.md, output.md 更新
- 推荐：添加 GENCODE 工作流使用指南
- 可选：参考 gencode-riboseqORFs 文档

---

## 🚀 快速入门（面向开发人员）

1. **阅读本摘要**（5 分钟）

2. **阅读详细计划**（20 分钟）
   - [SYNC_PLAN_20260130.md](./SYNC_PLAN_20260130.md) 第 1-3 章

3. **审查实现清单**（10 分钟）
   - [IMPLEMENTATION_CHECKLIST.md](./IMPLEMENTATION_CHECKLIST.md)
   - 确认第 1 周任务

4. **参考技术指南**（按需）
   - [TECHNICAL_REFERENCE.md](./TECHNICAL_REFERENCE.md)
   - 模块接口定义、数据格式规范

5. **开始第 1 周工作**
   - 创建模块框架
   - 准备 Python 脚本

---

## 📞 联系和支持

**文档作者**：GitHub Copilot  
**生成时间**：2026-01-30  
**版本**：1.0  

如有问题，请参考：
- 技术问题 → [TECHNICAL_REFERENCE.md](./TECHNICAL_REFERENCE.md)
- 实施问题 → [IMPLEMENTATION_CHECKLIST.md](./IMPLEMENTATION_CHECKLIST.md)
- 架构问题 → [SYNC_PLAN_20260130.md](./SYNC_PLAN_20260130.md)

---

## ✅ 验收标准清单

项目完成时应满足：

- [ ] 所有 5 个模块在 `modules/local/` 中实现
- [ ] 所有 3 个子工作流在 `subworkflows/local/` 中实现
- [ ] `workflows/riboseq/main.nf` 中成功集成
- [ ] 所有新参数在 `nextflow_schema.json` 中定义
- [ ] test profile 包含 GENCODE 注释功能测试
- [ ] README.md, docs/usage.md, docs/output.md 更新
- [ ] 关键输出文件（.orfs.gtf, .orfs.out, orf_counts_matrix.csv）生成验证
- [ ] 向后兼容性测试通过（不启用新功能时流程正常运行）
- [ ] Nextflow linting 通过（nf-test）
- [ ] 文档完整且无错误

---

**说明**：这三份文档（SYNC_PLAN_20260130.md, IMPLEMENTATION_CHECKLIST.md, TECHNICAL_REFERENCE.md）共同构成了完整的项目计划和技术参考。建议按以下方式使用：

- **项目管理** → SYNC_PLAN_20260130.md（架构、风险、验收标准）
- **开发执行** → IMPLEMENTATION_CHECKLIST.md（周度任务、工时估计）
- **技术实现** → TECHNICAL_REFERENCE.md（API、格式、调试）

祝项目顺利！ 🚀
