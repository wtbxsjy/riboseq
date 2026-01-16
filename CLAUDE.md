# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Repository Overview

**nf-core/riboseq** is a Nextflow DSL2 bioinformatics pipeline for analyzing ribosome profiling (Ribo-seq) data. The pipeline performs preprocessing, alignment, quality control, and ORF (open reading frame) prediction from ribosome footprinting experiments.

## Essential Commands

### Running the Pipeline

```bash
# Basic run with test profile
nextflow run . -profile test,docker --outdir results

# Run with custom samplesheet
nextflow run . -profile docker --input samplesheet.csv --outdir results

# Run with Singularity (common for HPC environments)
nextflow run . -profile test,singularity --outdir results

# Run with specific aligner (STAR is default, HISAT2 available)
nextflow run . -profile test,docker --aligner hisat2 --outdir results

# Resume a failed run
nextflow run . -profile test,docker --outdir results -resume
```

### Testing

```bash
# Run all nf-test tests
nf-test test --profile debug,test,docker --verbose

# Run tests for a specific module
nf-test test modules/local/orfquant/main.nf.test --profile debug,test,docker

# Run pipeline-level tests
nf-test test tests/pipeline/ --profile debug,test,docker
```

### Development Commands

```bash
# Lint the pipeline (nf-core standards compliance)
nf-core pipelines lint .

# Update the pipeline schema (after adding/modifying parameters)
nf-core pipelines schema build

# Format code with Prettier
prettier --write .

# Run pre-commit hooks
pre-commit run --all-files

# Clean test artifacts
rm -rf work/ .nf-test/ results*/
```

### Container Management

```bash
# Build custom ORFquant container with patches
apptainer build --fakeroot -F orfquant_patched.sif containers/Singularity.orfquant.patched.def

# Use custom container in pipeline
nextflow run . -profile test,singularity \
  --orfquant_container /path/to/orfquant_patched.sif
```

## Architecture Overview

### Pipeline Structure

The pipeline follows nf-core conventions with a modular architecture:

1. **Entry point**: `main.nf` - Orchestrates the main workflow
2. **Main workflow**: `workflows/riboseq/main.nf` - Contains the core RIBOSEQ workflow logic
3. **Subworkflows**: `subworkflows/local/` and `subworkflows/nf-core/` - Reusable workflow components
4. **Modules**: `modules/local/` and `modules/nf-core/` - Individual process definitions
5. **Configuration**: `conf/` directory contains resource configs, test profiles, and module-specific settings

### Data Flow

The pipeline supports two input modes:

1. **FASTQ mode** (default): Raw sequencing reads → preprocessing → alignment → QC/ORF prediction
2. **BAM mode**: Pre-aligned BAM files → QC/ORF prediction (skips preprocessing/alignment)

**Key processing stages:**

```
FASTQ Input → Merge → QC (FastQC) → UMI Extract → Trim → Contaminant Filter →
Strandedness Inference → Alignment (STAR/HISAT2) → Sort/Index → UMI Dedup →
BAM Filtering (sORF) → QC (RiboseQC, Ribo-TISH) → ORF Prediction → MultiQC Report
```

### Critical Components

**sORF BAM Filtering** (`modules/local/sorf_bam_filter/`):
- Applied before ORF prediction tools to ensure consistent input
- Filters for unique mapping reads (NH:i:1 tag or MAPQ threshold)
- Removes reads from mitochondrial/chloroplast/ambiguous contigs
- Filters by read length (default 28-30 nt for ribosome footprints)
- Controlled by `--sorf_filter*` parameters

**RiboseQC** (`subworkflows/local/riboseqc.nf`):
- Comprehensive Ribo-seq QC: P-site analysis, metagene profiles, periodicity
- Runs twice on `type=riboseq` samples: pre-filter (baseline) and post-filter QC
- Generates `*_for_ORFquant` files required by ORFquant

**ORFquant** (`subworkflows/local/orfquant.nf`):
- Splice-aware ORF detection and quantification
- **Important**: Requires RiboseQC output; skipping RiboseQC auto-skips ORFquant
- Uses custom patched container to avoid BiocGenerics::Position/combine namespace conflicts
- Patch details: Modified NAMESPACE to use selective ggplot2/gridExtra imports

### Sample Type Handling

The pipeline distinguishes between sample types via the `type` column in samplesheets:
- `riboseq`: Regular ribosome profiling data (undergoes full QC including RiboseQC)
- `tiseq`: Translation initiation sequencing (TI-seq)
- `rnaseq`: RNA-seq data (for future translational efficiency analysis)

**Important**: Only `type=riboseq` samples receive RiboseQC analysis and filtered/unfiltered QC comparisons.

### Aligner Support

Two aligners are available:
- **STAR** (default): Generates genome + transcriptome BAMs simultaneously
- **HISAT2**: Lower memory footprint (~2GB vs STAR's 30GB+); transcriptome index auto-built from GTF

Both produce genome and transcriptome alignments. RiboCode requires transcriptome alignments.

### ORF Prediction Tools

The pipeline integrates multiple ORF prediction tools:
- **Ribo-TISH** (default): De novo ORF prediction from alignment data
- **Ribotricer** (default): Reference-guided ORF detection
- **RiboseQC** (default): QC + P-site analysis
- **ORFquant** (default): Splice-aware quantification (requires RiboseQC)
- **RiboCode** (optional): Transcriptome-based ORF detection (`--run_ribocode`)
- **rp-bp** (optional): Bayesian ORF predictions (`--run_rpbp`)

**Design principle**: All ORF predictors run in **per-sample mode only** (no pooled/all-samples mode) to maintain runtime/memory control at scale.

## Key Configuration Files

- `nextflow.config`: Main pipeline configuration with all parameters
- `conf/base.config`: Default resource requirements (CPU/memory/time) with process labels
- `conf/modules.config`: Module-specific argument overrides
- `conf/test*.config`: Test profile configurations for different environments
- `nf-test.config`: nf-test framework configuration

## Important Implementation Notes

### Custom ORFquant Container

The pipeline uses a patched ORFquant package to resolve namespace conflicts:
- **Issue**: BiocGenerics exports `Position` and `combine`, which conflict with ggplot2/gridExtra
- **Solution**: Modified `NAMESPACE` to use selective imports (`importFrom()`) instead of full imports
- **Location**: `patched_packages/ORFquant-1.02/` or `ORFquant-1.1/`
- **Container**: Built via `containers/Singularity.orfquant.patched.def`

### BAM Input Mode

When providing pre-aligned BAMs:
- Set samplesheet with `bam` and `bam_index` columns instead of `fastq_1`/`fastq_2`
- Strandedness must be explicitly specified (no `auto` mode)
- UMI deduplication and RiboCode are automatically skipped
- sORF filtering still applies to ensure consistent ORF prediction inputs

### Species-Specific Contig Filtering

The default `--sorf_exclude_contigs_regex` targets common mitochondrial/chloroplast contigs:
- **Animals** (Gencode): Excludes `chrM`, `MT`, `chrUn_*`, `*_random`, `*_alt`, `*_fix`
- **Plants** (Ensembl): Also excludes `Mt` (mitochondrion), `Pt` (plastid/chloroplast)

**Override this parameter** if your reference uses different naming conventions.

### Test Profiles for Different Environments

- `test`: Minimal test with Docker/Singularity (2 CPU, 6GB RAM)
- `test_codespace`: GitHub Codespaces (HISAT2, low memory)
- `test_colab`: Google Colab environment
- `test_local_singularity`: Local Singularity testing with custom containers
- `test_full`: Full-size dataset test

## Development Workflow

### Adding New Parameters

1. Add default value to `params {}` block in `nextflow.config`
2. Run `nf-core pipelines schema build` to update `nextflow_schema.json`
3. Add validation logic if needed
4. Document in help text and update docs

### Adding New Modules

1. Check if module exists in nf-core/modules first
2. For custom modules: Create in `modules/local/<tool>/`
3. Follow nf-core module structure (main.nf, meta.yml, tests/)
4. Add process resource labels in `conf/base.config`
5. Add module-specific arguments in `conf/modules.config`
6. Write nf-test tests

### Submitting Changes

- Target branch: `dev` (not `main`/`master`)
- Ensure `nf-core pipelines lint` passes
- Ensure `nf-test test` passes with `--profile debug,test,docker`
- Update `CHANGELOG.md` following existing format
- PRs trigger GitHub Actions CI/linting checks

## Common Gotchas

1. **ORFquant requires RiboseQC**: If you skip RiboseQC (`--skip_riboseqc`), ORFquant is auto-skipped
2. **RiboCode needs transcriptome BAMs**: Only works with STAR or HISAT2 aligners
3. **STAR memory requirements**: Needs 30GB+ for human genome; use HISAT2 for low-memory environments
4. **sORF filtering is pre-QC**: Unfiltered BAMs run through first QC round, then filtered BAMs go to ORF predictors
5. **Gencode references**: Set `--gencode true` to handle Gencode-specific GTF attributes
6. **Nextflow version**: Requires Nextflow >= 24.04.2 (specified in manifest)
7. **BAM input strandedness**: Cannot use `auto` - must specify `forward`, `reverse`, or `unstranded`

## File Locations Reference

- Custom scripts: `bin/` (e.g., `filter_gtf.py`, `gtf2bed`)
- Helper test scripts: `scripts/singularity_single_tool_tests/`
- Pipeline tests: `tests/` and `*.nf.test` files throughout
- Example data: `example/` directory
- Documentation: `docs/` (usage.md, output.md)
