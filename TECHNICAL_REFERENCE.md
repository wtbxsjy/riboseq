# GENCODE ORF 注释工作流 - 技术参考指南

**项目**：nf-core/riboseq  
**版本**：1.0  
**日期**：2026-01-30

---

## 目录

1. [工作流架构](#工作流架构)
2. [数据格式规范](#数据格式规范)
3. [模块接口定义](#模块接口定义)
4. [子工作流接口定义](#子工作流接口定义)
5. [通道设计](#通道设计)
6. [参数配置](#参数配置)
7. [镜像和依赖](#镜像和依赖)
8. [错误处理](#错误处理)
9. [性能调优](#性能调优)
10. [故障排查](#故障排查)

---

## 工作流架构

```
输入（BAM 过滤后）
    ↓
[Ribo-TISH 预测] → 格式转换 (14) ⟶ ┐
[Ribotricer 预测] → 格式转换 (15) ⟶ ├→ 合并 (collectFile) → GENCODE ORF Mapper (16) → 统一注释 (GTF+OUT)
[ORFquant 预测] → 格式转换 (可选)   ⟶ ┘                                             ↓
                                                                         ┌→ 定量分析 (17) → 计数矩阵
                                                                         │
                                                                    [P-site Bedgraph]
                                                                         │
                                                                         └→ 工具比较 (18) → 高置信度 ORF
                                                                         
    → MultiQC 报告
```

### 处理阶段

| 阶段 | 组件 | 输入 | 输出 | 时间复杂度 |
|------|------|------|------|----------|
| 1. 格式转换 | 14, 15 模块 | 工具输出 | GENCODE FA/BED | O(n) per sample |
| 2. 文件合并 | collectFile | ORF 文件 | 单个合并文件 | O(1) |
| 3. 统一注释 | 16 模块 | 合并 FA/BED | GTF/OUT | O(n²) worst-case* |
| 4. 定量分析 | 17 模块 | GTF + Bedgraph | 计数矩阵 | O(m×n) |
| 5. 比较分析 | 18 模块 | OUT + 工具列表 | 比较表 | O(n) |

*O(n²)：ORF deduplication 涉及成对比较

---

## 数据格式规范

### 1. Ribo-TISH predict 输出 (`*_pred.txt`)

**源文件示例**：
```
transcript_id	strand	orf_id	start_codon_pos	orf_start	orf_end	orf_length
NM_001234567	+	ORF.1	100	100	450	350
NM_001234567	+	ORF.2	500	500	800	300
NM_001234568	-	ORF.1	200	200	650	450
```

**Nextflow 期望**：
- 制表符分隔
- 必需列：`transcript_id`, `strand`, `orf_id`, `orf_start`, `orf_end`
- 源自 Ribo-TISH predict 工具

**处理步骤**（ribotish_to_gencode.py）：
1. 解析输入文件
2. 按 `transcript_id` + `orf_id` 生成 ORF 名称
3. 从基因组 FASTA 中提取序列（坐标映射）
4. 翻译为蛋白质序列
5. 生成输出文件

### 2. Ribotricer translating_ORFs 输出 (`*_translating_ORFs.tsv`)

**源文件示例**：
```
transcript_id	strand	orf_id	orf_start	orf_end	orf_length	fscore	pscore	phase_score
NM_001234567	+	ORF_1	100	450	350	0.95	0.87	0.85
NM_001234568	-	ORF_1	200	650	450	0.92	0.80	0.88
```

**Nextflow 期望**：
- 制表符分隔
- 必需列：`transcript_id`, `strand`, `orf_id`, `orf_start`, `orf_end`, `pscore`（或 `phase_score`）
- 源自 Ribotricer detect-orfs 工具

**处理步骤**（ribotricer_to_gencode.py）：
1. 解析输入文件
2. **过滤**：`pscore >= --min-phase-score`（默认 0.5）
3. 坐标转换和序列提取（同 14）
4. 生成输出文件

### 3. GENCODE 格式 FASTA (`*.gencode.fa`)

**格式规范**：
```
>ORF_NAME--STUDY_ID
MPPQVALITQLCLGDSSCSPS...
>ORF_NAME2--STUDY_ID
MVPAPLSLLS...
```

**序列头规则**：
- `ORF_NAME`：唯一的 ORF 标识符（e.g., `NM_001234567_ORF.1`）
- `--`：分隔符（二个减号）
- `STUDY_ID`：样本 ID（来自 `--sample` 参数）

**序列内容**：
- 蛋白质序列（20 个标准氨基酸字母）
- 单行或多行（FASTA 格式）
- 最小长度检查（默认 ≥ 16 aa，脚本参数：`--min-length`）

### 4. GENCODE 格式 BED (`*.gencode.bed`)

**格式规范**（1-based，chr start-1 end）：
```
chr1	1000	1300	ORF_NAME--STUDY_ID	.	+
chr1	1500	1950	ORF_NAME2--STUDY_ID	.	+
chr2	5000	5450	ORF_NAME3--STUDY_ID	.	-
```

**字段说明**：
- chrom：染色体
- chromStart：1-based 起始位置
- chromEnd：1-based 结束位置（不含）
- name：ORF_NAME--STUDY_ID（与 FASTA header 一致）
- score：`.`（reserved）
- strand：`+` 或 `-`

**注意**：与标准 BED 不同，这里是 1-based（正常 BED 是 0-based）

### 5. 统一注释 GTF (`*.orfs.gtf`)

**格式规范**（GENCODE GTF 扩展）：
```
chr1	gencode	exon	1001	1300	.	+	0	gene_id "ORF_001--S1"; transcript_id "ORF_001--S1_exon1"; tool "Ribo-TISH"; tool_count "2";
chr1	gencode	CDS	1001	1300	.	+	0	gene_id "ORF_001--S1"; transcript_id "ORF_001--S1_exon1"; tool "Ribo-TISH,Ribotricer"; phase "0";
```

**关键属性**：
- `gene_id`：ORF 标识符
- `transcript_id`：exon 标识符（gene_id_exonN）
- `tool`：检测该 ORF 的工具列表（逗号分隔）
- `tool_count`：工具数量
- `phase`：CDS 框架（0, 1, 2）

**用途**：用于 featureCounts 或其他定量工具

### 6. 统一注释 OUT (`*.orfs.out`)

**格式规范**（制表符分隔）：
```
ORF_NAME	seq_length	seq_hash	ORFquant	RiboTISH	Ribotricer	RiboCode	rp-bp	mapped_transcript	mapped_gene
ORF_001--S1	300	abc123	1	1	0	1	0	NM_001234567	ENSG00000000001
ORF_002--S1	150	def456	0	1	1	0	1	NM_001234568	ENSG00000000002
ORF_003--S1	450	ghi789	1	1	1	1	1	NM_001234569	ENSG00000000003
```

**字段说明**：
- `ORF_NAME`：ORF 标识符
- `seq_length`：蛋白质长度（氨基酸）
- `seq_hash`：序列哈希（用于重复检测）
- `TOOL_1, TOOL_2, ...`：0/1 指示是否由该工具检测
- `mapped_transcript`：映射到的转录本 ID
- `mapped_gene`：映射到的基因 ID

**用途**：工具比较、高置信度 ORF 筛选

---

## 模块接口定义

### 模块 1：CONVERT_RIBOTISH_TO_GENCODE

**位置**：`modules/local/convert_ribotish_to_gencode/main.nf`

```groovy
process CONVERT_RIBOTISH_TO_GENCODE {
    tag "$meta.id"
    label 'process_low'
    container 'depot.galaxyproject.org/singularity/biopython:1.81'
    
    input:
    tuple val(meta), path(predict_txt), path(fasta)
    
    output:
    tuple val(meta), path("*.gencode.fa"), path("*.gencode.bed"), emit: gencode
    path "versions.yml", emit: versions
    
    script:
    """
    python3 ${projectDir}/scripts/gencode_converters/ribotish_to_gencode.py \\
        --pred ${predict_txt} \\
        --fasta ${fasta} \\
        --sample ${meta.id} \\
        --min-length ${params.min_orf_length} \\
        --output ${meta.id}
    
    cat > versions.yml <<-'EOF'
    "${task.process}":
        python: \$(python3 --version 2>&1 | awk '{print \$2}')
        biopython: \$(python3 -c 'import Bio; print(Bio.__version__)')
    EOF
    """
}
```

**参数**：
- `params.min_orf_length`：最小 ORF 长度（氨基酸，默认 16）

**预期行为**：
- 输入：Ribo-TISH `*_pred.txt`
- 输出：`${meta.id}.gencode.fa`, `${meta.id}.gencode.bed`
- 错误处理：文件不存在、格式错误时退出

---

### 模块 2：CONVERT_RIBOTRICER_TO_GENCODE

**位置**：`modules/local/convert_ribotricer_to_gencode/main.nf`

```groovy
process CONVERT_RIBOTRICER_TO_GENCODE {
    tag "$meta.id"
    label 'process_low'
    container 'depot.galaxyproject.org/singularity/biopython:1.81'
    
    input:
    tuple val(meta), path(tsv), path(fasta)
    
    output:
    tuple val(meta), path("*.gencode.fa"), path("*.gencode.bed"), emit: gencode
    path "versions.yml", emit: versions
    
    script:
    """
    python3 ${projectDir}/scripts/gencode_converters/ribotricer_to_gencode.py \\
        --tsv ${tsv} \\
        --fasta ${fasta} \\
        --sample ${meta.id} \\
        --min-length ${params.min_orf_length} \\
        --min-phase-score ${params.min_phase_score} \\
        --output ${meta.id}
    
    cat > versions.yml <<-'EOF'
    "${task.process}":
        python: \$(python3 --version 2>&1 | awk '{print \$2}')
        biopython: \$(python3 -c 'import Bio; print(Bio.__version__)')
    EOF
    """
}
```

**参数**：
- `params.min_orf_length`：最小 ORF 长度（默认 16）
- `params.min_phase_score`：最小 phase score（默认 0.5）

---

### 模块 3：GENCODE_ORF_MAPPER

**位置**：`modules/local/gencode_orf_mapper/main.nf`

```groovy
process GENCODE_ORF_MAPPER {
    label 'process_high'
    container 'URL_TO_GENCODE_ORF_MAPPER:latest'
    publishDir "${params.outdir}/gencode_orf_mapper", mode: 'copy'
    
    input:
    path merged_fa
    path merged_bed
    path ensembl_dir
    val project_id
    
    output:
    path "${project_id}.orfs.fa", emit: fa
    path "${project_id}.orfs.bed", emit: bed
    path "${project_id}.orfs.gtf", emit: gtf
    path "${project_id}.orfs.out", emit: out
    path "${project_id}.altmapped", emit: alt, optional: true
    path "${project_id}.unmapped", emit: unmapped, optional: true
    path "versions.yml", emit: versions
    
    script:
    """
    python3 ${ensembl_dir}/ORF_mapper_to_GENCODE_v1.1.py \\
        --fasta ${merged_fa} \\
        --bed ${merged_bed} \\
        --project ${project_id} \\
        --ensembl-dir ${ensembl_dir} \\
        --min-length ${params.min_orf_length} \\
        --collapse-threshold ${params.orf_collapse_threshold} \\
        --collapse-method ${params.orf_collapse_method}
    
    cat > versions.yml <<-'EOF'
    "${task.process}":
        gencode-orf-mapper: "v1.1"
    EOF
    """
}
```

**参数**：
- `params.min_orf_length`：最小长度（默认 16）
- `params.orf_collapse_threshold`：collapse 阈值（默认 0.9）
- `params.orf_collapse_method`：collapse 方法（`longest_string` 或 `psite_overlap`）

**验证**：
- Ensembl 目录必须包含 PROTEOME_FASTA, TRANSCRIPTOME_FASTA, SORTED_TRANSCRIPTOME_GTF

---

## 子工作流接口定义

### 子工作流 1：GENCODE_ANNOTATION

**位置**：`subworkflows/local/gencode_annotation/main.nf`

```groovy
workflow GENCODE_ANNOTATION {
    
    take:
    ch_ribotish_pred        // [ meta, predict_txt ]
    ch_ribotricer_tsv       // [ meta, tsv ]
    ch_fasta                // path(genome.fa)
    ch_ensembl_dir          // path(ensembl_dir)
    project_id              // string
    
    main:
    ch_versions = Channel.empty()
    
    // 格式转换
    CONVERT_RIBOTISH_TO_GENCODE(
        ch_ribotish_pred.combine([ch_fasta])
    )
    ch_versions = ch_versions.mix(CONVERT_RIBOTISH_TO_GENCODE.out.versions)
    
    CONVERT_RIBOTRICER_TO_GENCODE(
        ch_ribotricer_tsv.combine([ch_fasta])
    )
    ch_versions = ch_versions.mix(CONVERT_RIBOTRICER_TO_GENCODE.out.versions)
    
    // 收集和合并
    merged_fa = CONVERT_RIBOTISH_TO_GENCODE.out.gencode
        .map { meta, fa, bed -> fa }
        .mix(CONVERT_RIBOTRICER_TO_GENCODE.out.gencode.map { meta, fa, bed -> fa })
        .collectFile(name: "merged_all_tools.fa")
    
    merged_bed = CONVERT_RIBOTISH_TO_GENCODE.out.gencode
        .map { meta, fa, bed -> bed }
        .mix(CONVERT_RIBOTRICER_TO_GENCODE.out.gencode.map { meta, fa, bed -> bed })
        .collectFile(name: "merged_all_tools.bed")
    
    // ORF Mapper
    GENCODE_ORF_MAPPER(
        merged_fa,
        merged_bed,
        ch_ensembl_dir,
        project_id
    )
    ch_versions = ch_versions.mix(GENCODE_ORF_MAPPER.out.versions)
    
    emit:
    orfs_fa = GENCODE_ORF_MAPPER.out.fa
    orfs_bed = GENCODE_ORF_MAPPER.out.bed
    orfs_gtf = GENCODE_ORF_MAPPER.out.gtf
    orfs_out = GENCODE_ORF_MAPPER.out.out
    versions = ch_versions
}
```

**输入通道**：
- `ch_ribotish_pred`：Ribo-TISH 预测结果（可选）
- `ch_ribotricer_tsv`：Ribotricer 预测结果（可选）
- `ch_fasta`：基因组 FASTA
- `ch_ensembl_dir`：Ensembl 注释目录

**输出通道**：
- `orfs_gtf`：⭐ 用于定量的 GTF
- `orfs_out`：⭐ 用于工具比较的 OUT 表

---

## 通道设计

### 通道命名约定

```
ch_[data_type]_[stage]_[variant]

示例：
ch_ribotish_pred              # Ribo-TISH 预测输出
ch_ribotricer_tsv             # Ribotricer TSV 输出
ch_bams_for_sorf_prediction   # 用于 sORF 预测的 BAM
ch_orfs_gtf                   # 统一 ORF GTF
```

### 元数据(meta)结构

```groovy
meta = [
    id: "sample_id",           // 样本唯一标识
    sample_type: "riboseq",    // 样本类型
    strandedness: "forward",   // 链特异性
    // ... 其他字段
]
```

### 通道操作

**合并 FA/BED 文件**：
```groovy
// 从多个模块收集 FASTA 文件
merged_fa = CONVERT_RIBOTISH_TO_GENCODE.out.gencode
    .map { meta, fa, bed -> fa }
    .mix(
        CONVERT_RIBOTRICER_TO_GENCODE.out.gencode.map { meta, fa, bed -> fa }
    )
    .collectFile(name: "merged_all_tools.fa")
```

**条件流**：
```groovy
if (!params.skip_gencode_annotation && params.ensembl_dir) {
    GENCODE_ANNOTATION(...)
    ch_orfs_gtf = GENCODE_ANNOTATION.out.orfs_gtf
} else {
    ch_orfs_gtf = Channel.empty()
}
```

---

## 参数配置

### nextflow_schema.json 配置

```json
{
  "gencode_annotation": {
    "type": "object",
    "description": "GENCODE ORF annotation and quantification",
    "properties": {
      "skip_gencode_annotation": {
        "type": "boolean",
        "description": "Skip GENCODE annotation step",
        "default": false
      },
      "ensembl_dir": {
        "type": "string",
        "description": "Path to Ensembl annotation directory (required if skip_gencode_annotation=false)",
        "pattern": "^(/|~).*$"
      },
      "min_orf_length": {
        "type": "integer",
        "description": "Minimum ORF length (amino acids)",
        "default": 16,
        "minimum": 1,
        "maximum": 500
      },
      "min_phase_score": {
        "type": "number",
        "description": "Minimum phase score for Ribotricer (0.0-1.0)",
        "default": 0.5,
        "minimum": 0.0,
        "maximum": 1.0
      },
      "orf_collapse_threshold": {
        "type": "number",
        "description": "ORF similarity threshold for merging (0.0-1.0)",
        "default": 0.9,
        "minimum": 0.0,
        "maximum": 1.0
      },
      "orf_collapse_method": {
        "type": "string",
        "description": "Method for ORF deduplication",
        "default": "longest_string",
        "enum": ["longest_string", "psite_overlap"]
      },
      "skip_orf_quantification": {
        "type": "boolean",
        "description": "Skip P-site based quantification",
        "default": false
      },
      "skip_tool_comparison": {
        "type": "boolean",
        "description": "Skip tool comparison analysis",
        "default": false
      }
    }
  }
}
```

### conf/modules.config 配置

```groovy
// GENCODE 注释模块
withName: 'CONVERT_RIBOTISH_TO_GENCODE' {
    container = 'depot.galaxyproject.org/singularity/biopython:1.81'
    cpus = 2
    memory = '8 GB'
    time = '2 h'
}

withName: 'CONVERT_RIBOTRICER_TO_GENCODE' {
    container = 'depot.galaxyproject.org/singularity/biopython:1.81'
    cpus = 2
    memory = '8 GB'
    time = '2 h'
}

withName: 'GENCODE_ORF_MAPPER' {
    container = 'gencode-orf-mapper:v1.1'
    cpus = 8
    memory = '32 GB'
    time = '4 h'
    errorStrategy = { task.exitStatus == 1 ? 'retry' : 'terminate' }
    maxRetries = 2
}
```

---

## 镜像和依赖

### 容器镜像

| 模块 | 镜像 | 大小 | 来源 | 备注 |
|------|------|------|------|------|
| 14, 15 | biopython:1.81 | ~200 MB | Docker Hub | 支持 Singularity 转换 |
| 16 | gencode-orf-mapper:v1.1 | ~500 MB | 自定义 | 需要 gencode-riboseqORFs 库 |
| 17, 18 | rocker/tidyverse:4.3 | ~2 GB | Rocker Project | 包含 R + 常用包 |

### Python 依赖

**ribotish_to_gencode.py**：
- biopython >= 1.78
- pandas >= 1.0

**ribotricer_to_gencode.py**：
- biopython >= 1.78
- pandas >= 1.0

### R 依赖

**quantify_orfs_from_psites.R**：
```R
require(GenomicRanges)
require(rtracklayer)
require(data.table)
require(tidyverse)
```

**analyze_tool_comparison.R**：
```R
require(tidyverse)
require(UpSetR)
```

---

## 错误处理

### 预期错误和处理

| 错误 | 症状 | 检测方式 | 解决方案 |
|------|------|--------|--------|
| 文件不存在 | `FileNotFoundException` | 脚本检查 `if [[ ! -f $file ]]` | 验证输入路径 |
| 格式错误 | `ValueError` 或解析失败 | Python 异常捕获 | 检查源文件格式 |
| Ensembl 缺失 | `OSError` 在 ORF Mapper | 目录验证 | 运行 `prepare_ensembl_annotation.sh` |
| 内存不足 | Java/Python `OutOfMemoryError` | 进程监控 | 增加 `memory` 配置 |
| 镜像拉取失败 | `Singularity error` | 容器运行时异常 | 提前 pull 镜像或使用本地缓存 |

### 错误恢复策略

```groovy
// 在 modules.config 中
withName: 'GENCODE_ORF_MAPPER' {
    errorStrategy = { task.exitStatus == 1 ? 'retry' : 'terminate' }
    maxRetries = 2
}

// 主工作流中的条件检查
if (!params.ensembl_dir && !params.skip_gencode_annotation) {
    error "GENCODE annotation enabled but ensembl_dir not specified"
}
```

---

## 性能调优

### 计算资源配置

```
模块 14-15（格式转换）：
  - 2 CPUs
  - 8 GB RAM
  - 2 h 时限
  - 可并行运行（多样本）

模块 16（ORF Mapper）：
  - 8 CPUs（collapse 算法可并行）
  - 32 GB RAM（大型数据集需要更多）
  - 4 h 时限
  - 单工作流程（不可跳过）

模块 17-18（R 分析）：
  - 4 CPUs
  - 16 GB RAM
  - 2 h 时限
```

### 并行化机会

1. **格式转换并行**：多个工具同时转换
   ```groovy
   CONVERT_RIBOTISH_TO_GENCODE(...) // 可与下一行并行
   CONVERT_RIBOTRICER_TO_GENCODE(...)
   ```

2. **样本级并行**：多个样本独立处理
   ```groovy
   ch_ribotish_pred.map { ... }  // 每个样本分别处理
   ```

3. **collapse 算法优化**：使用哈希表加速序列比较
   - 当前：O(n²) 成对比较
   - 优化：使用 k-mer 前过滤 → O(n log n)

### 内存优化

**大型数据集分批处理**：
```bash
# ORF Mapper 支持分批参数（如果源码支持）
--batch-size 10000  # 每次处理 10K ORF
```

**流式处理而非一次性加载**：
```python
# 避免
orfs = list(read_fasta(file))  # 全部加载到内存

# 改进
for record in read_fasta(file):  # 逐条处理
    process(record)
```

---

## 故障排查

### 常见问题

**问题 1：Ensembl 目录验证失败**
```
[ERROR] Missing required Ensembl file: /path/to/Ens110/PROTEOME_FASTA
```
**解决方案**：
```bash
# 检查目录内容
ls -la /path/to/Ens110/
# 应包含：PROTEOME_FASTA, TRANSCRIPTOME_FASTA, SORTED_TRANSCRIPTOME_GTF 等
```

**问题 2：格式转换输出为空**
```
[ERROR] Output FASTA not generated
```
**解决方案**：
- 检查输入文件格式
- 验证最小长度过滤是否过严格
  ```bash
  grep ">" output.fa | wc -l  # 应 > 0
  ```

**问题 3：ORF Mapper 内存溢出**
```
java.lang.OutOfMemoryError: Java heap space
```
**解决方案**：
- 增加 `memory` 配置
- 启用分批处理
- 简化 collapse 参数（提高阈值）

**问题 4：R 包依赖缺失**
```
Error in library(GenomicRanges) : there is no package called 'GenomicRanges'
```
**解决方案**：
- 使用 Rocker 镜像（已包含依赖）
- 或在自定义脚本中调用 `BiocManager::install()`

### 调试模式

**启用详细日志**：
```bash
nextflow run main.nf -profile singularity --trace > trace.txt

# 检查每个步骤的资源使用
cat trace.txt | grep -E "task_id|cpu|memory|status"
```

**检查临时文件**：
```bash
# 查看 Nextflow work 目录
ls -la work/*/

# 检查特定步骤的输入/输出
ls -la work/52/*CONVERT_RIBOTISH*/
```

---

**文档版本**：1.0  
**最后更新**：2026-01-30  
**状态**：✅ 技术参考完成
