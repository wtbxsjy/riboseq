# gencode-riboseqORFs 集成实施计划

## 📅 项目时间线

**总预计时间**: 4-6 周
**优先级**: 中高
**复杂度**: 中等

---

## 🎯 第一阶段：容器和基础模块（Week 1-2）

### 任务 1.1: 构建 gencode-riboseqORFs 容器

#### Dockerfile 创建

**位置**: `containers/gencode-orf-mapper/Dockerfile`

```dockerfile
FROM continuumio/miniconda3:latest

LABEL \
    authors="nf-core/riboseq" \
    description="Container for gencode-riboseqORFs integration"

# 安装系统依赖
RUN apt-get update && apt-get install -y \
    build-essential \
    wget \
    git \
    && rm -rf /var/lib/apt/lists/*

# 安装 Conda 包
RUN conda install -c bioconda -c conda-forge \
    bedtools=2.30.0 \
    gffread=0.12.7 \
    && conda clean -a

# 安装 Python 依赖
RUN pip install biopython==1.79

# 克隆 gencode-riboseqORFs 项目
WORKDIR /opt
RUN git clone https://github.com/jorruior/gencode-riboseqORFs.git
WORKDIR /opt/gencode-riboseqORFs

# 设置环境变量
ENV PATH="/opt/gencode-riboseqORFs:${PATH}"

# 验证安装
RUN python3 --version && \
    bedtools --version && \
    gffread --version && \
    python3 -c "import Bio; print('Biopython:', Bio.__version__)"

WORKDIR /data
```

#### Singularity 定义文件

**位置**: `containers/gencode-orf-mapper/Singularity.def`

```singularity
Bootstrap: docker
From: continuumio/miniconda3:latest

%labels
    Author nf-core/riboseq
    Description Container for gencode-riboseqORFs integration

%post
    # 安装系统依赖
    apt-get update && apt-get install -y \
        build-essential \
        wget \
        git

    # 安装 Conda 包
    conda install -c bioconda -c conda-forge \
        bedtools=2.30.0 \
        gffread=0.12.7 \
        -y
    conda clean -a

    # 安装 Python 依赖
    pip install biopython==1.79

    # 克隆项目
    cd /opt
    git clone https://github.com/jorruior/gencode-riboseqORFs.git

%environment
    export PATH="/opt/gencode-riboseqORFs:${PATH}"

%runscript
    exec "$@"

%help
    Container for gencode-riboseqORFs integration with nf-core/riboseq

    Tools included:
    - Python 3 with Biopython
    - bedtools 2.30.0
    - gffread 0.12.7
    - gencode-riboseqORFs scripts
```

#### 构建脚本

**位置**: `containers/gencode-orf-mapper/build.sh`

```bash
#!/bin/bash
set -euo pipefail

# 构建 Docker 容器
docker build -t nfcore/gencode-orf-mapper:1.1.0 .

# 推送到 Docker Hub (可选)
# docker push nfcore/gencode-orf-mapper:1.1.0

# 构建 Singularity 容器
singularity build gencode-orf-mapper_1.1.0.sif Singularity.def

echo "✅ Container build completed!"
echo "Docker: nfcore/gencode-orf-mapper:1.1.0"
echo "Singularity: gencode-orf-mapper_1.1.0.sif"
```

---

### 任务 1.2: 创建 Ensembl 注释准备模块

**位置**: `modules/local/prepare_ensembl_annotation/main.nf`

```groovy
process DOWNLOAD_ENSEMBL_FILES {
    tag "Ensembl_${ensembl_release}_${genome_assembly}"
    label 'process_low'

    conda "conda-forge::wget=1.20.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/wget:1.20.3' :
        'biocontainers/wget:1.20.3' }"

    input:
    val ensembl_release
    val genome_assembly
    val species  // 'homo_sapiens' or 'mus_musculus'

    output:
    path "Ens${ensembl_release}/", emit: ensembl_dir
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def species_lower = species.toLowerCase()
    def species_caps = species.toLowerCase().replaceAll('_', ' ').split(' ').collect { it.capitalize() }.join('_')

    """
    mkdir -p Ens${ensembl_release}
    cd Ens${ensembl_release}

    # 下载 GTF
    wget -q http://ftp.ensembl.org/pub/release-${ensembl_release}/gtf/${species_lower}/${species_caps}.${genome_assembly}.${ensembl_release}.gtf.gz

    # 下载蛋白质序列
    wget -q http://ftp.ensembl.org/pub/release-${ensembl_release}/fasta/${species_lower}/pep/${species_caps}.${genome_assembly}.pep.all.fa.gz

    # 下载 cDNA 序列
    wget -q http://ftp.ensembl.org/pub/release-${ensembl_release}/fasta/${species_lower}/cdna/${species_caps}.${genome_assembly}.cdna.all.fa.gz

    # 下载 ncRNA 序列
    wget -q http://ftp.ensembl.org/pub/release-${ensembl_release}/fasta/${species_lower}/ncrna/${species_caps}.${genome_assembly}.ncrna.fa.gz

    # 解压
    gunzip -f *.gz

    # 排序 GTF
    sort -k1,1 -k4,4n -k5,5n ${species_caps}.${genome_assembly}.${ensembl_release}.gtf > ${species_caps}.${genome_assembly}.sorted.gtf

    # 合并转录本序列
    cat ${species_caps}.${genome_assembly}.cdna.all.fa ${species_caps}.${genome_assembly}.ncrna.fa | cut -d"." -f1,1 > ${species_caps}.${genome_assembly}.trans.fa

    # 处理蛋白质序列
    cut -d"." -f1,1 ${species_caps}.${genome_assembly}.pep.all.fa > tmpfile
    mv tmpfile ${species_caps}.${genome_assembly}.pep.all.fa

    # 下载 Transcript Support 信息
    wget -q -O ENST_support.txt 'http://ensembl.org/biomart/martservice?query=<Query virtualSchemaName="default" formatter="TSV" header="1" uniqueRows="0" count="" datasetConfigVersion="0.6"><Dataset name="${species_lower}_gene_ensembl" interface="default"><Attribute name="ensembl_transcript_id"/><Attribute name="transcript_tsl"/><Attribute name="transcript_appris"/></Dataset></Query>'

    cat <<-END_VERSIONS > ../versions.yml
    "${task.process}":
        ensembl_release: ${ensembl_release}
        genome_assembly: ${genome_assembly}
        wget: \$(wget --version | head -n 1 | sed 's/GNU Wget //')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p Ens${ensembl_release}
    touch Ens${ensembl_release}/stub.gtf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ensembl_release: ${ensembl_release}
        genome_assembly: ${genome_assembly}
    END_VERSIONS
    """
}

process CALCULATE_PSITE_BED {
    tag "Ensembl_${ensembl_release}"
    label 'process_low'

    conda "conda-forge::python=3.9 conda-forge::biopython=1.79"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.79' :
        'biocontainers/biopython:1.79' }"

    input:
    path ensembl_dir
    val ensembl_release

    output:
    path "Ens${ensembl_release}/", emit: ensembl_dir_complete
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # 复制脚本
    cp ${moduleDir}/assets/calculate_frame_bed.py .

    # 运行脚本
    cd ${ensembl_dir}
    python3 ../calculate_frame_bed.py *.sorted.gtf
    cd ..

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
        biopython: \$(python3 -c "import Bio; print(Bio.__version__)")
    END_VERSIONS
    """

    stub:
    """
    mkdir -p Ens${ensembl_release}
    touch Ens${ensembl_release}/psites.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "3.9.0"
    END_VERSIONS
    """
}
```

**辅助脚本**: `modules/local/prepare_ensembl_annotation/assets/calculate_frame_bed.py`

这个脚本从 gencode-riboseqORFs 项目复制：

```bash
cp /tmp/gencode-riboseqORFs/scripts/calculate_frame_bed.py \
   modules/local/prepare_ensembl_annotation/assets/
```

---

### 任务 1.3: 创建格式转换模块

#### Ribo-TISH 转换器

**位置**: `modules/local/convert_ribotish_to_gencode/main.nf`

```groovy
process CONVERT_RIBOTISH_TO_GENCODE {
    tag "${meta.id}"
    label 'process_low'

    conda "conda-forge::python=3.9 conda-forge::biopython=1.79"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.79' :
        'biocontainers/biopython:1.79' }"

    input:
    tuple val(meta), path(ribotish_predict), path(ribotish_quality)

    output:
    tuple val(meta), path("*.gencode.fa"), emit: fasta
    tuple val(meta), path("*.gencode.bed"), emit: bed
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def study_id = meta.study_id ?: meta.id

    """
    python3 ${moduleDir}/bin/ribotish_to_gencode.py \\
        --predict ${ribotish_predict} \\
        --quality ${ribotish_quality} \\
        --study_id ${study_id} \\
        --output_prefix ${prefix} \\
        --min_length ${params.gencode_orf_min_length ?: 16}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
        biopython: \$(python3 -c "import Bio; print(Bio.__version__)")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.gencode.fa
    touch ${prefix}.gencode.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "3.9.0"
    END_VERSIONS
    """
}
```

**转换脚本**: `modules/local/convert_ribotish_to_gencode/bin/ribotish_to_gencode.py`

```python
#!/usr/bin/env python3
"""
Convert Ribo-TISH output to gencode-riboseqORFs compatible format
"""

import argparse
from Bio import SeqIO
from Bio.Seq import Seq
import re

def parse_ribotish_predict(predict_file):
    """Parse Ribo-TISH predict output"""
    orfs = {}
    with open(predict_file, 'r') as f:
        for line in f:
            if line.startswith('#') or line.startswith('GenomePos'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 10:
                continue

            # 提取关键信息
            genome_pos = parts[0]  # chr:start-end:strand
            tid = parts[1]
            tis_type = parts[2]
            tis_group = parts[3]
            length = parts[4]
            # ... 其他字段

            orfs[genome_pos] = {
                'tid': tid,
                'type': tis_type,
                'length': length,
                'genome_pos': genome_pos
            }

    return orfs

def extract_sequence_from_bed(bed_file, fasta_file):
    """从 BED 坐标提取序列（需要基因组 FASTA）"""
    # 这里需要使用 bedtools getfasta 或类似工具
    # 简化版本：假设 Ribo-TISH 已经提供序列
    pass

def convert_to_gencode_format(orfs, study_id, output_prefix, min_length=16):
    """
    转换为 gencode-riboseqORFs 格式

    FASTA format: >{ORF_NAME}--{STUDY_ID}
    BED format (1-based): chr start end {ORF_NAME} {STUDY_ID} strand
    """

    fasta_output = f"{output_prefix}.gencode.fa"
    bed_output = f"{output_prefix}.gencode.bed"

    with open(fasta_output, 'w') as fa, open(bed_output, 'w') as bed:
        for orf_id, orf_data in orfs.items():
            # 解析基因组位置: chr:start-end:strand
            match = re.match(r'(.+):(\d+)-(\d+):([+-])', orf_id)
            if not match:
                continue

            chrom, start, end, strand = match.groups()
            start = int(start)
            end = int(end)

            # 过滤长度
            orf_length = (end - start) // 3
            if orf_length < min_length:
                continue

            # 构建 ORF 名称
            gene_name = orf_data.get('tid', 'unknown').split('.')[0]
            orf_name = f"{gene_name}_{start}_{orf_length}aa"

            # 写入 FASTA (需要序列信息)
            # 这里简化，实际需要从基因组提取
            sequence = orf_data.get('sequence', 'M' * orf_length + '*')
            fa.write(f">{orf_name}--{study_id}\n")
            fa.write(f"{sequence}\n")

            # 写入 BED (1-based)
            bed_start = start + 1  # 转换为 1-based
            bed_end = end  # end 在 0-based 半开区间中已经是正确的
            bed.write(f"{chrom}\t{bed_start}\t{bed_end}\t{orf_name}\t{study_id}\t{strand}\n")

def main():
    parser = argparse.ArgumentParser(description='Convert Ribo-TISH to gencode format')
    parser.add_argument('--predict', required=True, help='Ribo-TISH predict output')
    parser.add_argument('--quality', required=False, help='Ribo-TISH quality output')
    parser.add_argument('--study_id', required=True, help='Study identifier')
    parser.add_argument('--output_prefix', required=True, help='Output prefix')
    parser.add_argument('--min_length', type=int, default=16, help='Minimum ORF length')

    args = parser.parse_args()

    # 解析 Ribo-TISH 输出
    orfs = parse_ribotish_predict(args.predict)

    # 转换格式
    convert_to_gencode_format(orfs, args.study_id, args.output_prefix, args.min_length)

    print(f"✅ Converted {len(orfs)} ORFs to gencode format")

if __name__ == '__main__':
    main()
```

---

### 任务 1.4: 创建 GENCODE ORF 映射器模块

**位置**: `modules/local/gencode_orf_mapper/main.nf`

```groovy
process GENCODE_ORF_MAPPER {
    tag "${project_id}"
    label 'process_medium'

    conda "bioconda::bedtools=2.30.0 bioconda::gffread=0.12.7 conda-forge::biopython=1.79"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'file://gencode-orf-mapper_1.1.0.sif' :
        'nfcore/gencode-orf-mapper:1.1.0' }"

    input:
    path orfs_fasta
    path orfs_bed
    path ensembl_dir
    val project_id

    output:
    path "*.orfs.fa", emit: unified_fasta
    path "*.orfs.bed", emit: unified_bed
    path "*.orfs.gtf", emit: unified_gtf
    path "*.orfs.out", emit: unified_table
    path "*.altmapped", optional: true, emit: altmapped
    path "*.unmapped", optional: true, emit: unmapped
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def min_len = params.gencode_orf_min_length ?: 16
    def max_len = params.gencode_orf_max_length ?: 999999
    def collapse_thr = params.gencode_collapse_threshold ?: 0.9
    def collapse_method = params.gencode_collapse_method ?: 'longest_string'
    def add_cds = params.gencode_add_cds ? 'yes' : 'no'

    """
    python3 /opt/gencode-riboseqORFs/ORF_mapper_to_GENCODE_v1.1.py \\
        -d ${ensembl_dir} \\
        -f ${orfs_fasta} \\
        -b ${orfs_bed} \\
        -o ${project_id} \\
        -l ${min_len} \\
        -L ${max_len} \\
        -c ${collapse_thr} \\
        -m ${collapse_method} \\
        -C ${add_cds} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gencode_orf_mapper: "1.1.0"
        python: \$(python3 --version | sed 's/Python //')
        bedtools: \$(bedtools --version | sed 's/bedtools v//')
        gffread: \$(gffread --version 2>&1 | head -n1)
    END_VERSIONS
    """

    stub:
    """
    touch ${project_id}.orfs.fa
    touch ${project_id}.orfs.bed
    touch ${project_id}.orfs.gtf
    touch ${project_id}.orfs.out

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gencode_orf_mapper: "1.1.0"
    END_VERSIONS
    """
}
```

---

## 🎯 第二阶段：集成到主工作流（Week 3）

### 任务 2.1: 创建子工作流

**位置**: `subworkflows/local/gencode_orf_annotation.nf`

```groovy
//
// Subworkflow: Unified ORF annotation using gencode-riboseqORFs
//

include { DOWNLOAD_ENSEMBL_FILES } from '../../modules/local/prepare_ensembl_annotation/main'
include { CALCULATE_PSITE_BED } from '../../modules/local/prepare_ensembl_annotation/main'
include { CONVERT_RIBOTISH_TO_GENCODE } from '../../modules/local/convert_ribotish_to_gencode/main'
include { CONVERT_RIBOTRICER_TO_GENCODE } from '../../modules/local/convert_ribotricer_to_gencode/main'
include { GENCODE_ORF_MAPPER } from '../../modules/local/gencode_orf_mapper/main'

workflow GENCODE_ORF_ANNOTATION {
    take:
    ch_ribotish_orfs      // channel: [ meta, predict, quality ]
    ch_ribotricer_orfs    // channel: [ meta, orfs ]
    ensembl_release       // val: Ensembl release number
    genome_assembly       // val: Genome assembly (e.g., GRCh38)
    species               // val: Species name (e.g., homo_sapiens)
    project_id            // val: Project identifier

    main:
    ch_versions = Channel.empty()

    //
    // Prepare Ensembl annotation (run once)
    //
    DOWNLOAD_ENSEMBL_FILES(
        ensembl_release,
        genome_assembly,
        species
    )
    ch_versions = ch_versions.mix(DOWNLOAD_ENSEMBL_FILES.out.versions)

    CALCULATE_PSITE_BED(
        DOWNLOAD_ENSEMBL_FILES.out.ensembl_dir,
        ensembl_release
    )
    ch_versions = ch_versions.mix(CALCULATE_PSITE_BED.out.versions)

    //
    // Convert ORF formats to gencode-compatible format
    //
    ch_converted_orfs_fasta = Channel.empty()
    ch_converted_orfs_bed = Channel.empty()

    if (!params.skip_ribotish && ch_ribotish_orfs) {
        CONVERT_RIBOTISH_TO_GENCODE(
            ch_ribotish_orfs
        )
        ch_converted_orfs_fasta = ch_converted_orfs_fasta.mix(
            CONVERT_RIBOTISH_TO_GENCODE.out.fasta.map { meta, fa -> fa }
        )
        ch_converted_orfs_bed = ch_converted_orfs_bed.mix(
            CONVERT_RIBOTISH_TO_GENCODE.out.bed.map { meta, bed -> bed }
        )
        ch_versions = ch_versions.mix(CONVERT_RIBOTISH_TO_GENCODE.out.versions.first())
    }

    if (!params.skip_ribotricer && ch_ribotricer_orfs) {
        CONVERT_RIBOTRICER_TO_GENCODE(
            ch_ribotricer_orfs
        )
        ch_converted_orfs_fasta = ch_converted_orfs_fasta.mix(
            CONVERT_RIBOTRICER_TO_GENCODE.out.fasta.map { meta, fa -> fa }
        )
        ch_converted_orfs_bed = ch_converted_orfs_bed.mix(
            CONVERT_RIBOTRICER_TO_GENCODE.out.bed.map { meta, bed -> bed }
        )
        ch_versions = ch_versions.mix(CONVERT_RIBOTRICER_TO_GENCODE.out.versions.first())
    }

    //
    // Merge all ORF files
    //
    ch_merged_fasta = ch_converted_orfs_fasta.collectFile(
        name: "${project_id}_all_orfs.fa",
        newLine: false
    )

    ch_merged_bed = ch_converted_orfs_bed.collectFile(
        name: "${project_id}_all_orfs.bed",
        newLine: true
    )

    //
    // Run GENCODE ORF mapper
    //
    GENCODE_ORF_MAPPER(
        ch_merged_fasta,
        ch_merged_bed,
        CALCULATE_PSITE_BED.out.ensembl_dir_complete,
        project_id
    )
    ch_versions = ch_versions.mix(GENCODE_ORF_MAPPER.out.versions)

    emit:
    unified_fasta = GENCODE_ORF_MAPPER.out.unified_fasta
    unified_bed = GENCODE_ORF_MAPPER.out.unified_bed
    unified_gtf = GENCODE_ORF_MAPPER.out.unified_gtf
    unified_table = GENCODE_ORF_MAPPER.out.unified_table
    altmapped = GENCODE_ORF_MAPPER.out.altmapped
    unmapped = GENCODE_ORF_MAPPER.out.unmapped
    versions = ch_versions
}
```

---

### 任务 2.2: 集成到主工作流

**修改**: `workflows/riboseq/main.nf`

在文件末尾的 ORF 预测部分后添加：

```groovy
//
// SUBWORKFLOW: GENCODE ORF annotation (unified cross-tool annotation)
//
if (!params.skip_gencode_annotation) {

    // 收集所有 ORF 预测结果
    ch_ribotish_for_gencode = Channel.empty()
    ch_ribotricer_for_gencode = Channel.empty()

    if (!params.skip_ribotish) {
        ch_ribotish_for_gencode = ch_bams_for_sorf_prediction
            .join(RIBOTISH_QUALITY_RIBOSEQ.out.offset)
            .join(RIBOTISH_PREDICT_INDIVIDUAL.out.orfs)
            .map { meta, bam, bai, quality, predict ->
                [ meta + [ study_id: params.project_id ?: meta.id ], predict, quality ]
            }
    }

    if (!params.skip_ribotricer) {
        ch_ribotricer_for_gencode = RIBOTRICER_DETECTORFS.out.orfs
            .map { meta, orfs ->
                [ meta + [ study_id: params.project_id ?: meta.id ], orfs ]
            }
    }

    // 运行 GENCODE 注释
    GENCODE_ORF_ANNOTATION(
        ch_ribotish_for_gencode,
        ch_ribotricer_for_gencode,
        params.ensembl_release ?: 110,
        params.genome_assembly ?: 'GRCh38',
        params.species ?: 'homo_sapiens',
        params.project_id ?: 'riboseq_project'
    )
    ch_versions = ch_versions.mix(GENCODE_ORF_ANNOTATION.out.versions)

    // 添加到 MultiQC 报告
    ch_multiqc_files = ch_multiqc_files.mix(
        GENCODE_ORF_ANNOTATION.out.unified_table
    )
}
```

---

### 任务 2.3: 添加配置参数

**修改**: `nextflow.config`

```groovy
params {
    // ... 现有参数 ...

    // GENCODE ORF annotation
    skip_gencode_annotation     = false
    project_id                  = null  // 项目标识符，默认使用 outdir basename
    ensembl_release             = null  // 自动从参考基因组推断
    genome_assembly             = null  // 自动从参考基因组推断
    species                     = 'homo_sapiens'  // homo_sapiens, mus_musculus, etc.
    gencode_orf_min_length      = 16
    gencode_orf_max_length      = null
    gencode_collapse_threshold  = 0.9
    gencode_collapse_method     = 'longest_string'  // 或 'psite_overlap'
    gencode_add_cds             = false
}
```

---

## 🎯 第三阶段：测试和文档（Week 4）

### 任务 3.1: 创建测试配置

**位置**: `conf/test_gencode.config`

```groovy
/*
 * Test config for GENCODE ORF annotation
 */

params {
    config_profile_name = 'Test profile with GENCODE annotation'
    config_profile_description = 'Test nf-core/riboseq with GENCODE ORF annotation'

    // 输入数据
    input = 'https://raw.githubusercontent.com/nf-core/test-datasets/riboseq/samplesheet/samplesheet.csv'
    fasta = 'https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/riboseq_expression/Homo_sapiens.GRCh38.dna.chromosome.20.fa.gz'
    gtf = 'https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/riboseq_expression/Homo_sapiens.GRCh38.111_chr20.gtf'

    // GENCODE 注释设置
    skip_gencode_annotation = false
    project_id = 'test_riboseq'
    ensembl_release = 111
    genome_assembly = 'GRCh38'
    species = 'homo_sapiens'
    gencode_orf_min_length = 10  // 降低阈值以捕获更多测试 ORF

    // 其他测试设置
    skip_contaminant_filter = true
    min_trimmed_reads = 100
}

process {
    withName: 'GENCODE_ORF_MAPPER' {
        memory = 4.GB
        time = 1.h
    }
}
```

### 任务 3.2: 创建 nf-test 测试

**位置**: `subworkflows/local/gencode_orf_annotation/tests/main.nf.test`

```groovy
nextflow_workflow {

    name "Test Workflow GENCODE_ORF_ANNOTATION"
    script "../main.nf"
    workflow "GENCODE_ORF_ANNOTATION"

    test("Should run without failures") {

        when {
            workflow {
                """
                input[0] = Channel.empty() // ribotish_orfs
                input[1] = Channel.empty() // ribotricer_orfs
                input[2] = 111             // ensembl_release
                input[3] = 'GRCh38'        // genome_assembly
                input[4] = 'homo_sapiens'  // species
                input[5] = 'test_project'  // project_id
                """
            }
        }

        then {
            assert workflow.success
            assert workflow.trace.tasks().size() > 0
        }

    }
}
```

### 任务 3.3: 更新文档

**更新**: `README.md`

在 "Selecting ORF Prediction Tools" 章节后添加：

```markdown
### GENCODE ORF Annotation (Optional)

The pipeline can optionally integrate with [gencode-riboseqORFs](https://github.com/jorruior/gencode-riboseqORFs)
to provide unified, standardized ORF annotation across all prediction tools.

**Features**:
- ✅ Unified ORF annotation across multiple prediction tools
- ✅ ORF de-duplication and variant clustering
- ✅ Standardized classification (uORF, dORF, lncRNA-ORF, etc.)
- ✅ GENCODE/Ensembl compatible output
- ✅ Meta-analysis support (track ORFs across studies)

**Usage**:

```bash
# Enable GENCODE annotation
nextflow run nf-core/riboseq \
    -profile docker \
    --input samplesheet.csv \
    --outdir results \
    --skip_gencode_annotation false \
    --project_id my_study \
    --ensembl_release 110 \
    --genome_assembly GRCh38
```

**Parameters**:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--skip_gencode_annotation` | `false` | Skip GENCODE ORF annotation |
| `--project_id` | `null` | Project identifier (used in output) |
| `--ensembl_release` | `null` | Ensembl release number (auto-detected) |
| `--genome_assembly` | `null` | Genome assembly (auto-detected) |
| `--species` | `homo_sapiens` | Species name |
| `--gencode_orf_min_length` | `16` | Minimum ORF length (aa) |
| `--gencode_collapse_threshold` | `0.9` | Threshold for ORF merging |
| `--gencode_collapse_method` | `longest_string` | Method: `longest_string` or `psite_overlap` |

**Output Files**:

- `<project_id>.orfs.fa`: Unified ORF protein sequences
- `<project_id>.orfs.bed`: Unified ORF genomic coordinates
- `<project_id>.orfs.gtf`: Unified ORF GTF annotation
- `<project_id>.orfs.out`: Detailed ORF feature table

**ORF Classification**:

The tool classifies ORFs into 6 categories:
- **uORFs**: Upstream ORFs
- **uoORFs**: Upstream overlapping ORFs
- **dORFs**: Downstream ORFs
- **doORFs**: Downstream overlapping ORFs
- **intORFs**: Internal out-of-frame ORFs
- **lncRNA-ORFs**: Long non-coding RNA ORFs

**Supported Species**:
- Human (Homo sapiens) ✅
- Mouse (Mus musculus) ✅
- Other species: Contact developers for support

> [!NOTE]
> GENCODE annotation requires ~10-30 minutes additional runtime and downloads
> ~500MB of Ensembl annotation files (cached after first run).
```

---

## 📊 验证清单

### 功能测试

- [ ] Docker 容器构建成功
- [ ] Singularity 容器构建成功
- [ ] Ensembl 注释下载正常
- [ ] P-site BED 文件生成正确
- [ ] Ribo-TISH 格式转换正确
- [ ] Ribotricer 格式转换正确
- [ ] ORF 合并无重复
- [ ] GENCODE 映射器运行成功
- [ ] 输出文件格式正确

### 集成测试

- [ ] 完整流程运行（test_gencode profile）
- [ ] 与现有流程无冲突
- [ ] MultiQC 报告包含 GENCODE 结果
- [ ] 版本信息正确记录
- [ ] 错误处理正常

### 性能测试

- [ ] 小数据集（<100 ORFs）: <5分钟
- [ ] 中等数据集（100-1000 ORFs）: <15分钟
- [ ] 大数据集（>1000 ORFs）: <30分钟
- [ ] 内存使用 <8GB

### 文档完整性

- [ ] README 更新
- [ ] 参数说明完整
- [ ] 使用示例清晰
- [ ] CLAUDE.md 更新
- [ ] CHANGELOG.md 记录

---

## 🚀 快速启动指南

### 1. 构建容器

```bash
cd containers/gencode-orf-mapper
bash build.sh
```

### 2. 测试基础功能

```bash
# 测试 Ensembl 下载
nextflow run . -profile test_gencode --skip_ribotish --skip_ribotricer

# 测试格式转换（需要先有 ORF 预测结果）
nextflow run . -profile test_gencode --skip_gencode_annotation false
```

### 3. 完整测试

```bash
# 运行完整流程
nextflow run . -profile test_gencode,docker --outdir test_gencode_results

# 检查输出
ls test_gencode_results/gencode_annotation/
```

---

## 📝 下一步行动

1. **立即开始**: 创建容器 Dockerfile
2. **本周完成**: 基础模块实现
3. **下周目标**: 格式转换器开发
4. **月底交付**: 完整集成测试

需要我开始实现第一个模块吗？
