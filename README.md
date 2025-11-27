<h1>
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/images/nf-core-riboseq_logo_dark.png">
    <img alt="nf-core/riboseq" src="docs/images/nf-core-riboseq_logo_light.png">
  </picture>
</h1>

[![GitHub Actions CI Status](https://github.com/nf-core/riboseq/actions/workflows/ci.yml/badge.svg)](https://github.com/nf-core/riboseq/actions/workflows/ci.yml)
[![GitHub Actions Linting Status](https://github.com/nf-core/riboseq/actions/workflows/linting.yml/badge.svg)](https://github.com/nf-core/riboseq/actions/workflows/linting.yml)[![AWS CI](https://img.shields.io/badge/CI%20tests-full%20size-FF9900?labelColor=000000&logo=Amazon%20AWS)](https://nf-co.re/riboseq/results)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.10966364-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.10966364)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A524.04.2-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/nf-core/riboseq)

[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23riboseq-4A154B?labelColor=000000&logo=slack)](https://nfcore.slack.com/channels/riboseq)[![Follow on Twitter](http://img.shields.io/badge/twitter-%40nf__core-1DA1F2?labelColor=000000&logo=twitter)](https://twitter.com/nf_core)[![Follow on Mastodon](https://img.shields.io/badge/mastodon-nf__core-6364ff?labelColor=FFFFFF&logo=mastodon)](https://mstdn.science/@nf_core)[![Watch on YouTube](http://img.shields.io/badge/youtube-nf--core-FF0000?labelColor=000000&logo=youtube)](https://www.youtube.com/c/nf-core)

## Introduction

**nf-core/riboseq** is a bioinformatics pipeline for analysis of Ribo-seq data. It borrows heavily from nf-core/rnaseq in the preprocessing stages.

### Pipeline Overview

```mermaid
flowchart TB
    subgraph Input
        A[📁 FASTQ Files]
        B[📋 Sample Sheet]
    end

    subgraph Preprocessing["🔧 Preprocessing"]
        C[Merge FastQ<br/>cat]
        D[Read QC<br/>FastQC]
        E[UMI Extraction<br/>UMI-tools]
        F[Trimming<br/>Trim Galore! / fastp]
        G[Contaminant Removal<br/>Bowtie / Bowtie2]
        H[Strandedness Inference<br/>fq + Salmon]
    end

    subgraph Alignment["🧬 Alignment"]
        I{Aligner}
        J[STAR]
        K[HISAT2]
        L[Genome BAM +<br/>Transcriptome BAM]
        M[Sort & Index<br/>SAMtools]
        N[UMI Deduplication<br/>UMI-tools]
    end

    subgraph QC["📊 Quality Control"]
        O[RiboseQC<br/>Comprehensive QC +<br/>P-site Analysis]
        P[Ribo-TISH Quality<br/>P-site Offset]
    end

    subgraph ORF["🔬 ORF Prediction"]
        Q[Ribo-TISH Predict]
        R[Ribotricer]
        S[RiboCode]
        T[rp-bp]
        OQ[ORFquant<br/>Splice-aware<br/>Quantification]
    end

    subgraph Analysis["📈 Downstream Analysis"]
        U[anota2seq<br/>Translational Efficiency]
    end

    subgraph Output["📤 Output"]
        V[MultiQC Report]
        W[ORF Predictions]
        X[QC Reports]
    end

    A --> C
    B --> C
    C --> D --> E --> F --> G --> H
    H --> I
    I -->|--aligner star| J
    I -->|--aligner hisat2| K
    J --> L
    K --> L
    L --> M --> N

    N --> O
    N --> P
    N --> Q
    N --> R
    N --> S
    N --> T
    N --> U

    O --> OQ
    O --> X
    P --> X
    Q --> W
    R --> W
    S --> W
    T --> W
    OQ --> W
    U --> V
    X --> V
```

### Preprocessing Steps

1. Merge re-sequenced FastQ files ([`cat`](http://www.linfo.org/cat.html))
2. Sub-sample FastQ files and auto-infer strandedness ([`fq`](https://github.com/stjude-rust-labs/fq), [`Salmon`](https://combine-lab.github.io/salmon/))
3. Read QC ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
4. UMI extraction ([`UMI-tools`](https://github.com/CGATOxford/UMI-tools))
5. Adapter and quality trimming ([`Trim Galore!`](https://github.com/FelixKrueger/TrimGalore) or [`fastp`](https://github.com/OpenGene/fastp))
6. Removal of ribosomal/contaminant reads ([`Bowtie`](http://bowtie-bio.sourceforge.net/index.shtml) or [`Bowtie2`](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml))
7. Genome alignment with [`STAR`](https://github.com/alexdobin/STAR) or [`HISAT2`](http://daehwankimlab.github.io/hisat2/), outputting both genome and transcriptome alignments
8. Sort and index alignments ([`SAMtools`](https://sourceforge.net/projects/samtools/files/samtools/))
9. UMI-based deduplication ([`UMI-tools`](https://github.com/CGATOxford/UMI-tools))

### Ribo-seq Quality Control

1. **RiboseQC**: Comprehensive quality control including read length distribution, P-site analysis, metagene profiles, codon periodicity, and frame bias ([`RiboseQC`](https://github.com/ohlerlab/RiboseQC))
2. **Ribo-TISH Quality**: Check reads distribution around annotated protein coding regions, show frame bias and estimate P-site offset ([`Ribo-TISH`](https://github.com/zhpn1024/ribotish))

### ORF Prediction Tools

1. **Ribo-TISH** (default): Predict translated ORFs and translation initiation sites _de novo_ from alignment data ([`Ribo-TISH`](https://github.com/zhpn1024/ribotish))
2. **Ribotricer** (default): Derive candidate ORFs from reference data and detect translated ORFs ([`Ribotricer`](https://github.com/smithlabcode/ribotricer))
3. **RiboCode** (optional): Detect actively translating ORFs using transcriptome-aligned reads ([`RiboCode`](https://github.com/xztcwang/RiboCode))
4. **rp-bp** (optional): Ribosome profiling with Bayesian predictions for translated ORFs ([`rp-bp`](https://github.com/dieterich-lab/rp-bp))
5. **ORFquant** (default, requires RiboseQC): Splice-aware ORF detection and quantification at the single-ORF level ([`ORFquant`](https://github.com/lcalviell/ORFquant))

### Translational Efficiency Analysis

- **anota2seq** (optional): Study the dynamics of transcription and translation using matched RNA-seq and Ribo-seq data ([`anota2seq`](https://bioconductor.org/packages/release/bioc/html/anota2seq.html))

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

First, prepare a samplesheet with your input data that looks as follows:

`samplesheet.csv`:

```csv
sample,fastq_1,fastq_2,strandedness,type
CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz,forward,riboseq
```

Each row represents a fastq file (single-end) or a pair of fastq files (paired end). Each row should have a 'type' value of `riboseq`, `tiseq` or `rnaseq`. Future iterations of the workflow will conduct paired analysis of matched riboseq and rnaseq samples to accomplish analysis types such as 'translational efficiency, but in the current version you should set this to `riboseq` or `tiseq` for reglar Ribo-seq or TI-seq data respectively.

Now, you can run the pipeline using:

```bash
nextflow run nf-core/riboseq \
   -profile <docker/singularity/.../institute> \
   --input samplesheet.csv \
   --outdir <OUTDIR>
```

### Choosing an Aligner

The pipeline supports two aligners:

- **STAR** (default): `--aligner star`
- **HISAT2**: `--aligner hisat2`

Both aligners produce genome and transcriptome alignments. Pre-built indexes can be provided using:

```bash
# For STAR
--star_index /path/to/star/index

# For HISAT2
--hisat2_index /path/to/hisat2/genome/index
```

> [!NOTE]
> The HISAT2 transcriptome index is always auto-built from your GTF/FASTA to ensure compatibility.

### Selecting ORF Prediction Tools

By default, the pipeline runs Ribo-TISH, Ribotricer, RiboseQC, and ORFquant. Additional tools can be enabled:

```bash
# Enable RiboCode (requires STAR or HISAT2 transcriptome alignments)
--run_ribocode

# Enable rp-bp (requires contaminant FASTA)
--run_rpbp --contaminant_fasta /path/to/contaminants.fa
```

To skip specific tools:

```bash
--skip_ribotish
--skip_ribotricer
--skip_riboseqc
--skip_orfquant    # Note: ORFquant requires RiboseQC, skipping RiboseQC will also skip ORFquant
```

> [!NOTE]
> **ORFquant** uses the P-site analysis output from **RiboseQC** (`*_for_ORFquant` files). If you skip RiboseQC, ORFquant will also be automatically skipped.

### Including a translational efficiency analysis

![anota2seq - fold change plot](docs/images/fc.png)

In the translational efficiency analysis provided by [anota2seq](https://bioconductor.org/packages/release/bioc/html/anota2seq.html), we use matched pairs of Ribo-seq and RNA-seq data to study the relationship between transcription and translation as they differ between two treatment groups. For example the test data for this workflow has a contrasts file like:

```csv
id,variable,reference,target,batch,pair
treated_vs_control,treatment,control,treated,,pair
```

This describes how to compare groups of samples between treament groups, and between RNA-seq and Ribo-seq. In order the columns are:

- `id`: a unique identifier to use for the contrast
- 'variable`: which vaiable (column) of the sample sheet should be used to separate the treatment groups?
- `reference`: which value of the variable column should be used to select samples to be used as the reference/ base group?
- `target`: which value of the variable column should be used to select samples to be used as the target/treated group?
- `batch`: (optional) specify a variable in the sample sheet that defines sample batches
- `pair`: (optional) specify a variable in the sample shet that defines sample pairing between RNA-seq and Ribo-seq samples. If not specified, it is assumed that the two types of sample are ordered the same.

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/usage/getting_started/configuration#custom-configuration-files).

For more details and further functionality, please refer to the [usage documentation](https://nf-co.re/riboseq/usage) and the [parameter documentation](https://nf-co.re/riboseq/parameters).

## Pipeline output

To see the results of an example test run with a full size dataset refer to the [results](https://nf-co.re/riboseq/results) tab on the nf-core website pipeline page.
For more details about the output files and reports, please refer to the
[output documentation](https://nf-co.re/riboseq/output).

## Credits

nf-core/riboseq was originally written by [Jonathan Manning](https://github.com/pinin4fjords) (Bioinformatics Engineer at Seqera) with support from [Altos Labs](https://www.altoslabs.com/) and in discussion with [Felix Krueger](https://github.com/FelixKrueger) and [Christel Krueger](https://github.com/ChristelKrueger). We thank the following people for their input:

- Anne Bresciani (ZS)
- [Felipe Almeida](https://github.com/fmalmeida) (ZS)
- [Mikhail Osipovitch](https://github.com/mosi223) (ZS)
- [Edward Wallace](https://github.com/ewallace) (University of Edinburgh)
- [Jack Tierney](https://github.com/JackCurragh) (University College Cork)
- [Maxime U Garcia](https://github.com/maxulysse) (Seqera)

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#riboseq` channel](https://nfcore.slack.com/channels/riboseq) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citations

If you use nf-core/riboseq for your analysis, please cite it using the following doi: [10.5281/zenodo.10966364](https://doi.org/10.5281/zenodo.10966364)

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
