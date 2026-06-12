# gencode-riboseqORFs 集成可行性评估报告

## 项目概述

**项目**: [gencode-riboseqORFs](https://github.com/jorruior/gencode-riboseqORFs)
**作者**: Jorge Ruiz-Orera
**版本**: v1.1.0 (2022-09-19)
**发表文章**: https://doi.org/10.1038/s41587-022-01369-0
**用途**: 统一多个 Ribo-seq 研究的 ORF 预测结果，并将其映射到 GENCODE/Ensembl 注释

## 核心功能

### 主要目标
1. **统一 ORF 注释**: 合并来自不同 Ribo-seq 研究的 ORF 预测结果
2. **去冗余**: 识别并合并相似/重叠的 ORF 变体
3. **标准化映射**: 将所有 ORF 映射到特定的 Ensembl 版本（如 v101）
4. **分类注释**: 根据 ORF 类型进行分类（uORF, dORF, lncRNA-ORF 等）

### ORF 分类系统（6 大类）

根据 README，工具支持以下 ORF 生物型分类：

1. **uORFs** (Upstream ORFs): 上游 ORF
2. **uoORFs** (Upstream overlapping ORFs): 上游重叠 ORF
3. **dORFs** (Downstream ORFs): 下游 ORF
4. **doORFs** (Downstream overlapping ORFs): 下游重叠 ORF
5. **intORFs** (Internal out-of-frame ORFs): 内部框外 ORF
6. **lncRNA-ORFs**: 长非编码 RNA 的 ORF

## 技术架构

### 依赖项
```
- Python 3
- Biopython
- gffread v0.10.1+
- bedtools v2.27.1+
```

### 输入要求

#### 1. Ensembl 注释文件（通过脚本自动下载）
```bash
bash scripts/retrieve_ensembl_data.sh <ENSEMBL_RELEASE> <GENOME_ASSEMBLY>
```
需要的文件：
- `PROTEOME_FASTA`: 蛋白质序列
- `TRANSCRIPTOME_FASTA`: 转录本序列（mRNA + ncRNA）
- `SORTED_TRANSCRIPTOME_GTF`: 排序的 GTF 文件
- `TRANSCRIPT_SUPPORT`: TSL 和 APPRIS 支持信息（从 Biomart 下载）
- `PSITES_BED`: P-site 坐标（通过脚本生成）

#### 2. ORF 预测文件（来自各个 Ribo-seq 研究）

**FASTA 格式**（蛋白质序列）:
```
>A1BG_58858387_46aa--6_Chen2020
MQPRAQGAVGVLRSAGDSGLAPSPPVAAQGRGLWGAGEASLIPPRN*
>A1BG_58858945_25aa--6_Chen2020
MPSCAARDPSPTSPSSCCARARRRP*
```
- 命名格式: `{ORF_NAME}--{STUDY_ID}`
- 必须包含终止密码子 `*`

**BED 格式**（基因组坐标，1-based）:
```
19  58346882  58347022  A1BG_58858387_46aa  6_Chen2020  -
19  58347503  58347580  A1BG_58858945_25aa  6_Chen2020  -
```
- 第 4 列: ORF 名称
- 第 5 列: 研究标识
- 必须是 1-based 坐标

### 输出文件

1. **`<OUT_NAME>.orfs.fa`**: 统一后的 ORF 蛋白质序列
2. **`<OUT_NAME>.orfs.bed`**: 统一后的 ORF 基因组坐标（1-based）
3. **`<OUT_NAME>.orfs.gtf`**: 统一后的 ORF GTF 注释（不含假基因和 CDS）
4. **`<OUT_NAME>.orfs.out`**: 详细的 ORF 特征表格

#### 输出表格字段（关键信息）
- `orf_id`: 唯一 ORF 标识符
- `orf_biotype`: ORF 类型（6 大类）
- `gene_biotype`: 宿主基因类型
- `trans`: 主要宿主转录本 ID（Ensembl）
- `gene`: 宿主基因 ID
- `gene_name`: 宿主基因名称
- `n_datasets`: 检测到该 ORF 的数据集数量
- `X_study`: 每个研究的检测标记（1/0）
- `n_variants`: 变体数量
- `all_orf_names`: 原始 ORF 名称（来自各研究）

### 核心算法

#### ORF 去冗余方法（`-m` 参数）

**方法 1: `longest_string`**（默认，用于 Phase I）
- 比较重叠 ORF 对的最长共同氨基酸序列
- 如果共同序列长度/短 ORF 长度 ≥ 阈值（默认 0.9），则合并
- 要求起始或终止密码子相同
- 选择最长的 ORF 作为代表

**方法 2: `psite_overlap`**
- 比较 P-site 位置的重叠比例
- 如果共享 P-site 比例 ≥ 阈值（默认 0.9），则合并
- 速度较慢但更精确

## 与 nf-core/riboseq 的集成分析

### ✅ 集成优势

#### 1. **完美的功能互补**
- **nf-core/riboseq**: 专注于从原始数据到 ORF 预测的完整流程
- **gencode-riboseqORFs**: 专注于跨研究 ORF 统一注释和标准化
- 集成后可实现：数据分析 → ORF 预测 → 统一注释 → 标准化输出

#### 2. **输入格式兼容性高**

**nf-core/riboseq 现有工具输出**:

| 工具 | 输出格式 | 与 gencode-riboseqORFs 兼容性 |
|------|---------|-------------------------------|
| **Ribo-TISH** | BED + 序列信息 | ✅ 高度兼容（需格式转换） |
| **Ribotricer** | TSV（含坐标和序列） | ✅ 可转换为 BED + FASTA |
| **RiboCode** | GTF + 序列 | ✅ 可转换为 BED + FASTA |
| **rp-bp** | BED + 序列 | ✅ 高度兼容 |
| **ORFquant** | GTF + 定量信息 | ✅ 可转换（已有坐标） |

#### 3. **标准化的 ORF 分类**
- 提供 6 种标准 ORF 类型分类
- 与 GENCODE 官方注释兼容
- 便于后续功能注释和比较分析

#### 4. **跨样本整合能力**
- 自动识别不同样本间的相同 ORF
- 提供每个 ORF 在多个数据集中的检测情况
- 支持 meta-analysis 和可重复性评估

### ⚠️ 集成挑战

#### 1. **物种和注释版本限制**

**当前限制**:
- 主要针对人类（Homo sapiens）
- 脚本硬编码了 Ensembl 人类数据库 URL
- 默认使用 Ensembl v101

**解决方案**:
```bash
# 需要修改脚本支持其他物种
# 或为每个物种准备相应的注释文件
- 小鼠: ftp.ensembl.org/pub/release-XXX/gtf/mus_musculus/
- 其他物种: 类似路径结构
```

**影响**: ⚠️ 中等 - 需要适配脚本或准备额外注释

#### 2. **坐标系统差异**

**gencode-riboseqORFs 要求**:
- **1-based** BED 坐标（不同于标准 BED 的 0-based）

**nf-core/riboseq 工具输出**:
- 大多数工具使用 **0-based** 坐标（标准）
- 需要转换步骤

**解决方案**:
```python
# 转换函数（0-based → 1-based）
def convert_to_1based(bed_0based):
    # start: 0-based → 1-based (加 1)
    # end: 0-based 半开区间已经等于 1-based 闭区间
    bed_1based_start = bed_0based_start + 1
    bed_1based_end = bed_0based_end  # 保持不变
```

**影响**: ✅ 低 - 可通过模块自动处理

#### 3. **命名格式要求**

**要求格式**: `{ORF_NAME}--{STUDY_ID}`

**nf-core/riboseq 当前输出**:
- Ribo-TISH: `chr1:12345-67890:+`
- Ribotricer: `ENST00000123456_ORF_1`
- RiboCode: 自定义格式

**解决方案**:
```bash
# 创建标准化模块
- 从 meta 对象提取样本 ID 作为 STUDY_ID
- 标准化 ORF 命名：{GENE}_{POSITION}_{LENGTH}aa--{SAMPLE_ID}
```

**影响**: ✅ 低 - 通过预处理模块解决

#### 4. **多种 ORF 预测工具的整合复杂性**

**当前情况**:
- nf-core/riboseq 支持 5 种不同的 ORF 预测工具
- 每种工具输出格式不同
- 需要为每种工具创建格式转换器

**解决方案**:
```nextflow
// 创建统一的格式转换子工作流
include { CONVERT_RIBOTISH_TO_GENCODE_FORMAT } from './modules/local/format_converters'
include { CONVERT_RIBOTRICER_TO_GENCODE_FORMAT } from './modules/local/format_converters'
include { CONVERT_RIBOCODE_TO_GENCODE_FORMAT } from './modules/local/format_converters'
// ...
```

**影响**: ⚠️ 中等 - 需要开发多个转换模块

#### 5. **计算资源需求**

**gencode-riboseqORFs 特点**:
- Python 脚本，单线程运行
- 对大量 ORF 的去冗余计算密集
- P-site overlap 方法特别慢

**优化策略**:
```groovy
// Nextflow 配置
process {
    withName: 'GENCODE_ORF_MAPPER' {
        cpus = 1  // 单线程
        memory = { 8.GB * task.attempt }
        time = { 4.h * task.attempt }
    }
}
```

**影响**: ⚠️ 中等 - 大规模数据集可能耗时较长

#### 6. **依赖项管理**

**额外依赖**:
- gffread (0.10.1+)
- bedtools (v2.27.1+)
- Biopython

**解决方案**:
```dockerfile
# 创建专用容器
FROM biocontainers/bedtools:2.30.0
RUN apt-get update && apt-get install -y \
    python3-pip \
    gffread
RUN pip3 install biopython

# 或使用 Conda
conda:
  - bedtools=2.30.0
  - gffread=0.12.7
  - biopython=1.79
```

**影响**: ✅ 低 - 已有容器化方案

### 🎯 推荐的集成策略

#### 方案 1: 完整集成（推荐）✨

**实现步骤**:

1. **创建格式转换模块**（`modules/local/orf_format_converters/`）
   ```
   - ribotish_to_gencode.nf
   - ribotricer_to_gencode.nf
   - ribocode_to_gencode.nf
   - rpbp_to_gencode.nf
   - orfquant_to_gencode.nf
   ```

2. **创建 gencode-riboseqORFs 包装模块**（`modules/local/gencode_orf_mapper/`）
   ```nextflow
   process GENCODE_ORF_MAPPER {
       conda "bioconda::bedtools=2.30.0 bioconda::gffread=0.12.7 conda-forge::biopython=1.79"
       container "biocontainers/gencode-orf-mapper:1.1.0"

       input:
       path orfs_fasta
       path orfs_bed
       path ensembl_dir
       val study_id

       output:
       path "*.orfs.fa", emit: fasta
       path "*.orfs.bed", emit: bed
       path "*.orfs.gtf", emit: gtf
       path "*.orfs.out", emit: table

       script:
       """
       python3 /opt/ORF_mapper_to_GENCODE_v1.1.py \\
           -d ${ensembl_dir} \\
           -f ${orfs_fasta} \\
           -b ${orfs_bed} \\
           -o ${study_id} \\
           -l 16 \\
           -c 0.9 \\
           -m longest_string
       """
   }
   ```

3. **创建 Ensembl 注释准备子工作流**（`subworkflows/local/prepare_ensembl_annotation.nf`）
   ```nextflow
   workflow PREPARE_ENSEMBL_ANNOTATION {
       take:
       ensembl_release
       genome_assembly

       main:
       // 下载并准备 Ensembl 文件
       DOWNLOAD_ENSEMBL_FILES(ensembl_release, genome_assembly)
       PREPARE_PSITE_BED(DOWNLOAD_ENSEMBL_FILES.out.gtf)

       emit:
       ensembl_dir = DOWNLOAD_ENSEMBL_FILES.out.dir
   }
   ```

4. **集成到主工作流**（`workflows/riboseq/main.nf`）
   ```nextflow
   // 在所有 ORF 预测完成后
   if (!params.skip_gencode_annotation) {

       // 准备 Ensembl 注释（只运行一次）
       PREPARE_ENSEMBL_ANNOTATION(
           params.ensembl_release,
           params.genome_assembly
       )

       // 转换各工具输出为统一格式
       if (!params.skip_ribotish) {
           CONVERT_RIBOTISH_TO_GENCODE_FORMAT(
               RIBOTISH_PREDICT.out.orfs
           )
       }

       if (!params.skip_ribotricer) {
           CONVERT_RIBOTRICER_TO_GENCODE_FORMAT(
               RIBOTRICER_DETECTORFS.out.orfs
           )
       }

       // 收集所有转换后的 ORF
       ch_all_orfs_fasta = Channel.empty()
           .mix(CONVERT_RIBOTISH_TO_GENCODE_FORMAT.out.fasta)
           .mix(CONVERT_RIBOTRICER_TO_GENCODE_FORMAT.out.fasta)
           .collect()

       ch_all_orfs_bed = Channel.empty()
           .mix(CONVERT_RIBOTISH_TO_GENCODE_FORMAT.out.bed)
           .mix(CONVERT_RIBOTRICER_TO_GENCODE_FORMAT.out.bed)
           .collect()

       // 合并所有 ORF 文件
       MERGE_ORF_FILES(ch_all_orfs_fasta, ch_all_orfs_bed)

       // 运行 GENCODE ORF 映射器
       GENCODE_ORF_MAPPER(
           MERGE_ORF_FILES.out.fasta,
           MERGE_ORF_FILES.out.bed,
           PREPARE_ENSEMBL_ANNOTATION.out.ensembl_dir,
           params.project_id
       )
   }
   ```

5. **添加新参数**（`nextflow.config`）
   ```groovy
   params {
       // GENCODE ORF annotation
       skip_gencode_annotation    = false
       ensembl_release            = 110  // 或根据参考基因组自动推断
       genome_assembly            = 'GRCh38'  // 或 GRCm39 (小鼠)
       gencode_orf_min_length     = 16
       gencode_collapse_threshold = 0.9
       gencode_collapse_method    = 'longest_string'  // 或 'psite_overlap'
       project_id                 = null  // 研究标识符
   }
   ```

**优点**:
- ✅ 完全自动化，从原始数据到标准化 ORF 注释
- ✅ 跨工具 ORF 去冗余和统一
- ✅ 符合 GENCODE 标准，便于发表和共享
- ✅ 支持 meta-analysis（多数据集比较）

**缺点**:
- ⚠️ 开发工作量较大（约 5-7 个新模块）
- ⚠️ 需要为每个支持的物种准备 Ensembl 注释
- ⚠️ 增加运行时间（约 10-30 分钟，取决于 ORF 数量）

#### 方案 2: 可选模块集成（灵活）

**实现方式**:
- 将 gencode-riboseqORFs 作为**可选的后处理步骤**
- 默认跳过（`--skip_gencode_annotation true`）
- 用户需要时手动启用

**适用场景**:
- 用户只想快速获得 ORF 预测结果
- 非人类/小鼠物种（暂不支持）
- 单样本分析（无需跨样本整合）

**优点**:
- ✅ 不影响现有流程
- ✅ 用户可选择是否使用
- ✅ 降低默认运行时间

**缺点**:
- ⚠️ 用户需要了解何时启用
- ⚠️ 不默认提供标准化输出

#### 方案 3: 外部工具推荐（最简单）

**实现方式**:
- 在文档中推荐 gencode-riboseqORFs 作为后处理工具
- 提供使用示例和格式转换脚本
- 不集成到流程中

**优点**:
- ✅ 无开发工作
- ✅ 保持流程简洁

**缺点**:
- ❌ 用户需要手动操作
- ❌ 失去自动化优势

### 📊 集成价值评估

| 评估维度 | 分数 (1-10) | 说明 |
|---------|------------|------|
| **科学价值** | 9/10 | 提供 GENCODE 标准化注释，便于发表和共享 |
| **用户需求** | 8/10 | 多样本研究和 meta-analysis 有强需求 |
| **技术可行性** | 7/10 | 需要适度开发工作，但无重大技术障碍 |
| **维护成本** | 6/10 | 需要跟踪 Ensembl 版本更新 |
| **性能影响** | 7/10 | 增加少量运行时间，可接受 |
| **兼容性** | 8/10 | 与现有工具高度兼容 |

**总体评分**: **7.5/10** - **强烈推荐集成**

### 🔧 具体实现建议

#### 阶段 1: 基础集成（2-3 周）

1. **Week 1**: 创建容器和基础模块
   - 构建 gencode-riboseqORFs 容器（Docker + Singularity）
   - 实现 GENCODE_ORF_MAPPER 模块
   - 创建 PREPARE_ENSEMBL_ANNOTATION 子工作流

2. **Week 2**: 格式转换器开发
   - 实现 Ribo-TISH 格式转换器
   - 实现 Ribotricer 格式转换器
   - 添加单元测试

3. **Week 3**: 集成测试
   - 端到端测试（test profile）
   - 文档更新
   - 参数调优

#### 阶段 2: 扩展功能（1-2 周）

4. **Week 4**: 其他工具支持
   - RiboCode 格式转换器
   - rp-bp 格式转换器
   - ORFquant 格式转换器

5. **Week 5**:多物种支持
   - 小鼠（Mus musculus）支持
   - 自动推断 Ensembl 版本
   - 物种特异性配置

#### 阶段 3: 优化和文档（1 周）

6. **Week 6**:
   - 性能优化
   - 用户文档和示例
   - 集成到 nf-core 模板

### ⚡ 快速启动方案（POC）

如果你想快速验证可行性，可以先实现一个**概念验证版本**:

```bash
# 1. 手动测试 gencode-riboseqORFs
cd /tmp/gencode-riboseqORFs
bash scripts/retrieve_ensembl_data.sh 110 GRCh38

# 2. 从 nf-core/riboseq 导出一个 ORF 预测结果
# (假设已有 Ribo-TISH 输出)

# 3. 转换格式并运行
python3 scripts/convert_ribotish_to_gencode.py \
    ribotish_output.txt \
    > orfs.fa

python3 ORF_mapper_to_GENCODE_v1.1.py \
    -d Ens110 \
    -f orfs.fa \
    -b orfs.bed \
    -a ATG

# 4. 检查输出
head orfs.orfs.out
```

如果 POC 成功，再考虑完整集成。

## 总结和建议

### ✅ 可行性结论

**高度可行**，推荐集成到 nf-core/riboseq 流程中。

### 🎯 推荐方案

**方案 1（完整集成）**，因为：

1. **科学价值高**: 提供 GENCODE 标准化的 ORF 注释
2. **技术可行**: 无重大技术障碍，主要是工程工作
3. **用户体验好**: 一键式从原始数据到标准化注释
4. **未来扩展性强**: 支持 meta-analysis 和跨研究比较

### 📋 行动清单

- [ ] 与 nf-core/riboseq 维护者讨论集成计划
- [ ] 创建 GitHub Issue 跟踪开发进度
- [ ] 构建 gencode-riboseqORFs 容器（Docker + Singularity）
- [ ] 实现格式转换模块（从最常用的 Ribo-TISH 开始）
- [ ] 创建测试数据集
- [ ] 编写集成文档和用户指南
- [ ] 提交 Pull Request 到 nf-core/riboseq

### ⚠️ 注意事项

1. **版本兼容性**: 需要明确支持的 Ensembl 版本范围
2. **物种支持**: 初期聚焦人类和小鼠，其他物种逐步添加
3. **许可证**: 确认 gencode-riboseqORFs 许可证兼容 nf-core
4. **性能测试**: 大规模数据集（>1000 个 ORF）的性能验证
5. **用户反馈**: Beta 测试阶段收集用户意见

### 🤝 需要的帮助

如果你决定推进这个集成，我建议：

1. **联系原作者**: 与 Jorge Ruiz-Orera 沟通，获取技术支持和建议
2. **社区讨论**: 在 nf-core Slack #riboseq 频道讨论集成计划
3. **测试用户**: 招募几个用户进行 beta 测试
4. **文献调研**: 查看其他研究如何使用 gencode-riboseqORFs

---

**评估完成日期**: 2026-01-17
**评估人**: Claude (Anthropic AI)
**建议审阅**: nf-core/riboseq 维护团队
