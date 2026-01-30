# Nextflow 流程同步实现检查清单

**项目**：nf-core/riboseq  
**日期**：2026-01-30  
**目标**：逐步实现 singularity_single_tool_tests 中的新增功能

---

## 📋 快速概览

### 新增功能数量
- **新模块**：5 个
- **新子工作流**：3 个
- **新参数**：8+ 个
- **预计工作量**：3-4 周

### 优先级
- 🔴 **P0（第 1 周）**：14-16 脚本（格式转换和统一注释）
- 🟡 **P1（第 2 周）**：17-19 脚本（定量和比较分析）
- 🟢 **P2（第 3 周）**：性能优化和文档

---

## ✅ 逐周实现检查清单

### 第 1 周：基础准备

**目标**：创建模块框架，准备 Python 脚本和镜像

- [ ] **任务 1.1** - 创建模块目录结构
  ```bash
  mkdir -p modules/local/{convert_ribotish_to_gencode,convert_ribotricer_to_gencode,gencode_orf_mapper,quantify_orfs_psites,analyze_tool_comparison}
  ```
  - [ ] 创建 `main.nf` 文件（暂时为空）
  - [ ] 创建 `meta.yml` 文件
  - **预计时间**：30 分钟

- [ ] **任务 1.2** - 准备 Python 转换脚本
  - [ ] `scripts/gencode_converters/ribotish_to_gencode.py`
    - 输入：Ribo-TISH `*_pred.txt`
    - 输出：`*.gencode.fa`, `*.gencode.bed`
    - 测试：使用脚本 `14_ribotish_to_gencode.sh` 的逻辑
    - **预计时间**：2 小时
  - [ ] `scripts/gencode_converters/ribotricer_to_gencode.py`
    - 输入：Ribotricer `*_translating_ORFs.tsv`
    - 输出：`*.gencode.fa`, `*.gencode.bed`
    - 包含 phase score 过滤
    - **预计时间**：1.5 小时
  - [ ] `scripts/gencode_converters/README.md`
    - 说明每个脚本的用途和参数
    - **预计时间**：30 分钟

- [ ] **任务 1.3** - 验证 Singularity 镜像
  - [ ] 测试 `biopython:1.81` 镜像
    ```bash
    singularity pull biopython_1.81.sif docker://continuumio/miniconda3:latest
    singularity exec biopython_1.81.sif python3 -c "import Bio; print(Bio.__version__)"
    ```
    - **预计时间**：1 小时
  - [ ] 确认 gencode-orf-mapper 镜像可用
    - 检查 URL 或构建自定义镜像
    - **预计时间**：1.5-2 小时（取决于是否需要构建）

- [ ] **任务 1.4** - 单元测试数据准备
  - [ ] 从 `test_data/06_ribotish_predict/` 获取样本 Ribo-TISH 输出
  - [ ] 从 `test_data/08_ribotricer_detectorfs/` 获取样本 Ribotricer 输出
  - [ ] 创建小型测试数据集（1-2 ORF）
    - **预计时间**：1 小时

**周总结**
- 工作量：~9 小时
- 关键交付物：
  - 5 个模块框架
  - 2 个 Python 脚本（初稿）
  - 镜像验证完成
  - 测试数据准备

---

### 第 2 周：核心模块实现

**目标**：实现三个关键模块（14-16），单独测试

- [ ] **任务 2.1** - 实现 `modules/local/convert_ribotish_to_gencode/main.nf`
  - [ ] 创建 Groovy 模块框架
  ```groovy
  process CONVERT_RIBOTISH_TO_GENCODE {
      container 'depot.galaxyproject.org/singularity/biopython:1.81'
      input:
          tuple val(meta), path(predict_txt), path(fasta)
      output:
          tuple val(meta), path("*.gencode.fa"), path("*.gencode.bed")
      ...
  }
  ```
  - [ ] 调用 Python 脚本
  - [ ] 处理版本文件
  - [ ] 测试（使用 test data）
  - **预计时间**：2 小时

- [ ] **任务 2.2** - 实现 `modules/local/convert_ribotricer_to_gencode/main.nf`
  - [ ] 创建模块（类似 2.1）
  - [ ] 添加 phase score 过滤参数
  - [ ] 测试
  - **预计时间**：1.5 小时

- [ ] **任务 2.3** - 实现 `modules/local/gencode_orf_mapper/main.nf`
  - [ ] 创建模块
  ```groovy
  process GENCODE_ORF_MAPPER {
      container 'gencode-orf-mapper:custom'
      input:
          path(merged_fa)
          path(merged_bed)
          path(ensembl_dir)
      output:
          path("*.orfs.gtf"), path("*.orfs.out"), path("*.orfs.fa")
      ...
  }
  ```
  - [ ] 验证 Ensembl 目录结构
  - [ ] 处理 collapse parameters
  - [ ] 测试
  - **预计时间**：2.5 小时

- [ ] **任务 2.4** - 单元测试
  - [ ] 测试 14：Ribo-TISH 转换
    - [ ] 对比脚本 `14_ribotish_to_gencode.sh` 的输出
    - 验证 FASTA header 格式（`>ORF_NAME--STUDY_ID`）
    - 验证 BED 1-based 坐标
    - **预计时间**：1 小时
  - [ ] 测试 15：Ribotricer 转换
    - [ ] 验证 phase score 过滤生效
    - **预计时间**：1 小时
  - [ ] 测试 16：ORF Mapper
    - [ ] 验证 GTF 和 OUT 文件结构
    - [ ] 检查 collapse 结果
    - **预计时间**：1.5 小时

**周总结**
- 工作量：~12 小时
- 关键交付物：
  - 3 个完整模块（14-16）
  - 单元测试通过
  - 输出格式验证完成

---

### 第 3 周：子工作流集成

**目标**：创建 3 个子工作流，整合模块

- [ ] **任务 3.1** - 创建 `subworkflows/local/gencode_annotation/main.nf`
  ```groovy
  workflow GENCODE_ANNOTATION {
      take:
          ch_ribotish_pred
          ch_ribotricer_tsv
          ch_fasta
          ch_ensembl_dir
      
      main:
          CONVERT_RIBOTISH_TO_GENCODE(ch_ribotish_pred)
          CONVERT_RIBOTRICER_TO_GENCODE(ch_ribotricer_tsv)
          
          // 合并所有工具的 FASTA/BED
          merged_fa = ...
          merged_bed = ...
          
          GENCODE_ORF_MAPPER(merged_fa, merged_bed, ch_ensembl_dir)
      
      emit:
          orfs_gtf = GENCODE_ORF_MAPPER.out.gtf
          orfs_out = GENCODE_ORF_MAPPER.out.out
          versions = ...
  }
  ```
  - [ ] 实现合并逻辑（collectFile/concat）
  - [ ] 处理多样本场景
  - [ ] 添加元数据传播
  - **预计时间**：3 小时

- [ ] **任务 3.2** - 创建 `subworkflows/local/quantify_orfs/main.nf`
  - [ ] 包装 R 脚本 `scripts/R/quantify_orfs_from_psites.R`
    - 输入：GTF, bedgraph 目录
    - 输出：counts_matrix.csv, counts_tpm.tsv
  - [ ] 处理 Rscript 调用（类似 MultiQC）
  - **预计时间**：2 小时

- [ ] **任务 3.3** - 创建 `subworkflows/local/tool_comparison/main.nf`
  - [ ] 包装 R 脚本 `scripts/R/analyze_tool_comparison.R`
    - 输入：orfs.out, 工具列表
    - 输出：高置信度 ORF, upset plot, 特异性 ORF
  - [ ] 处理可视化输出
  - **预计时间**：2 小时

- [ ] **任务 3.4** - 子工作流集成测试
  - [ ] 测试 GENCODE_ANNOTATION 端到端
    - [ ] 提供 test 数据
    - [ ] 验证输出文件
    - **预计时间**：2 小时
  - [ ] 测试 QUANTIFY_ORFS
    - [ ] 验证 counts_matrix 格式
    - **预计时间**：1 小时
  - [ ] 测试 TOOL_COMPARISON
    - [ ] 验证高置信度 ORF 识别
    - **预计时间**：1 小时

**周总结**
- 工作量：~14 小时
- 关键交付物：
  - 3 个完整子工作流
  - 子工作流集成测试通过
  - R 脚本包装完成

---

### 第 4 周：主工作流集成

**目标**：在 `workflows/riboseq/main.nf` 中集成，完整流程测试

- [ ] **任务 4.1** - 更新 `nextflow_schema.json`
  ```json
  {
    "gencode_annotation": {
      "type": "object",
      "description": "GENCODE ORF annotation and quantification",
      "properties": {
        "skip_gencode_annotation": {
          "type": "boolean",
          "default": false,
          "description": "Skip GENCODE annotation"
        },
        "ensembl_dir": {
          "type": "string",
          "description": "Path to Ensembl annotation directory"
        },
        "orf_collapse_threshold": {
          "type": "number",
          "default": 0.9,
          "description": "ORF collapse threshold (0.0-1.0)"
        },
        "skip_orf_quantification": {
          "type": "boolean",
          "default": false
        },
        "skip_tool_comparison": {
          "type": "boolean",
          "default": false
        }
      }
    }
  }
  ```
  - **预计时间**：1 小时

- [ ] **任务 4.2** - 在 `workflows/riboseq/main.nf` 中添加导入
  ```groovy
  include { GENCODE_ANNOTATION } from '../../subworkflows/local/gencode_annotation'
  include { QUANTIFY_ORFS } from '../../subworkflows/local/quantify_orfs'
  include { TOOL_COMPARISON } from '../../subworkflows/local/tool_comparison'
  ```
  - **预计时间**：0.5 小时

- [ ] **任务 4.3** - 在工作流末尾添加 GENCODE 流程调用
  - [ ] 添加条件检查（`!params.skip_gencode_annotation`）
  - [ ] 传入正确的输入通道
  - [ ] 混合版本信息
  - [ ] 混合 MultiQC 文件
  - **预计时间**：2 小时

- [ ] **任务 4.4** - 添加定量和比较分析调用
  - [ ] 条件检查
  - [ ] 通道连接
  - **预计时间**：1.5 小时

- [ ] **任务 4.5** - 更新 `conf/modules.config`
  ```groovy
  withName: 'CONVERT_RIBOTISH_TO_GENCODE' {
      container = 'depot.galaxyproject.org/singularity/biopython:1.81'
  }
  withName: 'GENCODE_ORF_MAPPER' {
      container = 'path/to/gencode-orf-mapper.sif'
      cpus = 8
      memory = '32 GB'
  }
  ...
  ```
  - **预计时间**：1 小时

- [ ] **任务 4.6** - 完整流程测试
  - [ ] 使用现有 test profile 运行
    ```bash
    nextflow run main.nf -profile test,singularity
    ```
    - [ ] 验证所有新增步骤执行
    - [ ] 验证输出文件生成
    - **预计时间**：2 小时
  - [ ] 创建 test_gencode profile（可选）
    - [ ] 提供最小 Ensembl 注释目录
    - **预计时间**：1 小时

- [ ] **任务 4.7** - 回归测试
  - [ ] 运行不使用 GENCODE 注释的流程
    - [ ] 验证现有 ORF 预测结果不变
    - [ ] 验证流程仍能完成
    - **预计时间**：1 小时

**周总结**
- 工作量：~10 小时
- 关键交付物：
  - 主工作流集成完成
  - 完整流程测试通过
  - 配置文件更新
  - 回归测试通过

---

### 第 5 周：文档和优化

**目标**：文档更新，性能优化，发布准备

- [ ] **任务 5.1** - 更新文档
  - [ ] `README.md`
    - [ ] 工作流图添加 GENCODE 注释部分
    - [ ] 添加"ORF 预测工具"章节链接
    - **预计时间**：1 小时
  - [ ] `docs/usage.md`
    - [ ] 添加 GENCODE 注释工作流使用示例
    ```bash
    nextflow run nf-core/riboseq \
      -profile singularity \
      --input samplesheet.csv \
      --ensembl_dir /path/to/Ens110 \
      --outdir results
    ```
    - [ ] 说明如何准备 Ensembl 目录
    - **预计时间**：1.5 小时
  - [ ] `docs/output.md`
    - [ ] 添加新输出文件说明
      - `orfs.gtf`：用于定量
      - `orfs.out`：工具检测矩阵
      - `orf_counts_matrix.csv`：DESeq2 输入
      - `high_confidence_orfs.tsv`：高置信度 ORF
    - **预计时间**：1 小时

- [ ] **任务 5.2** - 更新 `CHANGELOG.md`
  ```markdown
  ## [Unreleased]
  
  ### Added
  - GENCODE ORF annotation workflow (14-16 scripts integration)
    - Format conversion for Ribo-TISH and Ribotricer
    - Multi-tool ORF unification and deduplication
    - P-site quantification module
    - Tool comparison analysis
  - New parameters: skip_gencode_annotation, ensembl_dir, orf_collapse_threshold
  - New outputs: orfs.gtf, orfs.out, orf_counts_matrix.csv, high_confidence_orfs.tsv
  ```
  - **预计时间**：0.5 小时

- [ ] **任务 5.3** - 性能优化
  - [ ] 并行化 FASTA/BED 转换
    - 当前：逐个转换
    - 优化：同时处理多个工具输出
    - **预计时间**：1.5 小时
  - [ ] 优化 ORF Mapper 内存使用
    - 检查是否可以分批处理大型数据集
    - **预计时间**：1 小时
  - [ ] 缓存优化
    - 确保 Ensembl 文件不重复下载
    - **预计时间**：1 小时

- [ ] **任务 5.4** - 最后检查
  - [ ] Linting（nf-test, ShellCheck）
    - **预计时间**：1 小时
  - [ ] 文档完整性检查
    - **预计时间**：0.5 小时
  - [ ] 功能测试清单
    ```
    - [ ] 单样本 GENCODE 注释
    - [ ] 多样本 GENCODE 注释
    - [ ] 跳过 GENCODE 注释（params.skip_gencode_annotation=true）
    - [ ] 无 Ensembl 目录时的错误处理
    - [ ] ORF 定量
    - [ ] 工具比较分析
    ```
    - **预计时间**：2 小时

- [ ] **任务 5.5** - 版本发布
  - [ ] 更新版本号（e.g., v2.1.0）
  - [ ] 创建 Git tag
  - [ ] 准备发布说明
  - **预计时间**：1 小时

**周总结**
- 工作量：~12 小时
- 关键交付物：
  - 完整文档更新
  - 性能优化完成
  - 发布准备就绪

---

## 📊 总结

### 时间投入
- **第 1 周**：~9 小时（基础准备）
- **第 2 周**：~12 小时（核心模块）
- **第 3 周**：~14 小时（子工作流）
- **第 4 周**：~10 小时（主工作流集成）
- **第 5 周**：~12 小时（文档和优化）
- **总计**：~57 小时（约 7-8 个工作日）

### 交付物
| 周次 | 交付物数量 | 质量指标 |
|-----|----------|--------|
| 1 | 5 个模块框架 + 2 个 Python 脚本 | 脚本通过初步测试 |
| 2 | 3 个完整模块 | 单元测试 100% 通过 |
| 3 | 3 个子工作流 | 集成测试通过 |
| 4 | 主工作流集成 | 完整流程通过 |
| 5 | 文档 + 优化 | 发布就绪 |

### 风险和缓解措施

| 风险 | 影响 | 缓解措施 |
|-----|------|---------|
| Ensembl 镜像不可用 | 无法运行 ORF Mapper | 提前准备备用镜像或本地构建 |
| R 包依赖冲突 | 定量/比较分析失败 | 创建专用 R 环境/容器 |
| 大型数据集性能问题 | 流程超时 | 实施分批处理和缓存 |
| 输出格式兼容性 | 下游工具无法读取 | 与 gencode-riboseqORFs 库作者沟通 |

---

## 🎯 下一步

1. **批准此计划** - 获得项目管理部门同意
2. **分配资源** - 指派开发人员到各个任务
3. **创建 GitHub Issue** - 为每个任务创建跟踪问题
4. **开始第 1 周工作** - 启动模块框架和脚本准备

---

**文档版本**：1.0  
**最后更新**：2026-01-30  
**状态**：✅ 检查清单就绪，等待执行
