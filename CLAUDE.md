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

**riboWaltz** (`modules/local/ribowaltz/`):
- P-site offset calculation, metagene/codon/CDS analysis (complementary to RiboseQC)
- Runs **before** RiboseQC — its per-length P-site offsets serve as fallback when RiboseQC `P_sites_calcs` is empty
- **Uses transcriptome BAM** in alignment mode (STAR/HISAT2), genome BAM in BAM-input mode
- **Bioc 3.20**: Replaces `GenomicFeatures::makeTxDbFromGFF()` with `txdbmaker::makeTxDbFromGFF()`; custom annotation builder strips transcript ID version suffixes (e.g. `.11`) so BAM reads match GTF transcripts
- Container: build from `containers/Singularity.ribowaltz.def`; the runtime install script adds `txdbmaker` if missing from the base image
- Parameters: `--skip_ribowaltz`, `--ribowaltz_read_lengths [28,29,30]`, `--extra_ribowaltz_args`, `--ribowaltz_container`

**P-site offset fallback chain** (`EXTRACT_RL_CUTOFF` → `PREPARE_FOR_ORFQUANT_CORRECTED` → ORFquant):
1. RiboseQC `P_sites_calcs` valid → use RiboseQC offset (unchanged)
2. RiboseQC data empty + riboWaltz available → use riboWaltz `corrected_offset_from_5` per read length (improved accuracy)
3. Neither available → hardcoded defaults `28-32 → 12` (behaviour unchanged from original)

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

### ORF Unification and Classification (Post-Processing)

After per-sample ORF prediction, the pipeline runs two additional stages:

**ORF Unification** (`scripts/unify_orf_predictions.py`):
- Merges Ribo-TISH / Ribotricer / ORFquant results across all samples into a single non-redundant set.
- Deduplication order: exact-match → frame-aware → overlap grouping (selects representative per group).
- Outputs: `unified_orfs.bed` (BED12), `unified_orfs.gtf`, `unified_orfs.metadata.tsv`, `unified_orfs.stats.txt`.
- Tool names stored in `ORFCandidate.sources` use **capitalised/hyphenated** forms: `'Ribo-TISH'`, `'Ribotricer'`, `'ORFquant'` — use these exact strings when parsing the `sources` set or `tools` metadata column.
- Skip with `--skip_unify_orf_predictions true`.

**ORF Classification** (`modules/local/classify_orfs/`, `scripts/classify_orfs_wrapper.py`):
Three classifiers run in parallel (all enabled by default):

1. **GENCODE/Ensembl mode** (`CLASSIFY_ORFS_GENCODE`, `scripts/gencode-riboseqORFs/ORF_mapper_to_GENCODE_v1.1.py`):
   - Maps ORFs against transcriptome; assigns `orf_biotype`: `CDS`, `dORF`, `uORF`, `doORF`, `uoORF`, `intORF`, `lncRNA`.
   - Requires `--orf_classify_ensembl_dir` pointing to a directory with standardised symlinks (`TRANSCRIPTOME_FASTA`, `SORTED_TRANSCRIPTOME_GTF`, `PROTEOME_FASTA`, `TRANSCRIPT_SUPPORT`, `PSITES_BED`).
   - Container: `--gencode_orf_mapper_container` (needs bedtools + BioPython).
   - **Input format**: `classify_orfs_wrapper.py` auto-converts the BED12 output to BED6 (sample ID in col[4]) and translates nucleotide sequences to protein FASTA with key `{orf_id}--{sample_id}` as required by the mapper.
   - Output: `gencode_results.orfs.out`, `gencode_results.orfs.gtf`.

2. **ORFquant mode** (`CLASSIFY_ORFS_ORFQUANT`, `scripts/class_orf/run_orfquant_classify.R` → `orfquant_orf_classify.R`):
   - Classifies ORFs against the reference GTF at genomic (`ORF_category_Gen`), transcript (`ORF_category_Tx`), and best-isoform (`ORF_category_Tx_compatible`) levels.
   - **Transcript-space projection**: `project_to_tx_coords()` (top-level function in `orfquant_orf_classify.R`) maps ORF genomic blocks through the transcript exon chain to 1-based transcript coordinates, correctly handling multi-exon ORFs for both strands. This replaces a previous genomic-approximation that was wrong for ~18% of ORFs.
   - `normalize_annotation()` loads both `exon` and `CDS` features; returns `exon_txs` (per-transcript exon GRanges) and `cds_txs_tx_coords` (CDS bounds in transcript space).
   - `ORF_category_Tx_compatible` = best classification across all transcripts of the gene.
   - Output: `orfquant_classification.tsv`.

3. **ORF-type mode** (`CLASSIFY_ORFS_ORF_TYPE`, `scripts/class_orf/class_ORFtype.py`):
   - Gene-level classification: `canonical_CDS`, `uORF`, `dORF`, `overlap_uORF`, etc.
   - Output: `orftype_classification.tsv`.

Skip all classification with `--skip_orf_classification true`.

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
8. **GENCODE classifier requires its own container** (`--gencode_orf_mapper_container`): needs bedtools and BioPython; the `unify_orf` container does NOT include bedtools.
9. **GENCODE Ensembl directory** (`--orf_classify_ensembl_dir`): must contain five standardised symlinks — `TRANSCRIPTOME_FASTA`, `SORTED_TRANSCRIPTOME_GTF`, `PROTEOME_FASTA`, `TRANSCRIPT_SUPPORT`, `PSITES_BED` — created by the reference-preparation scripts.
10. **ORF tool name capitalisation**: `ORFCandidate.sources` and the unified metadata `tools` column use `'Ribo-TISH'`, `'Ribotricer'`, `'ORFquant'` (not lowercase). Use these exact strings when filtering or counting by tool.
11. **riboWaltz needs transcriptome BAM** for accurate per-transcript P-site analysis. In alignment mode it receives transcriptome BAMs; in BAM-input mode it falls back to genome BAMs.
12. **riboWaltz Bioc 3.20 compatibility**: `create_annotation()` calls `GenomicFeatures::makeTxDbFromGFF()` which is defunct in Bioc 3.20. The patched R script uses `txdbmaker::makeTxDbFromGFF()` instead. If the container lacks `txdbmaker`, it is installed at runtime (~2-3 min overhead). Rebuild the container (`containers/Singularity.ribowaltz.def`) to eliminate this.
13. **Transcript ID version mismatch**: GTF IDs (e.g. `ENST00001008.11`) have different version suffixes than BAM IDs (e.g. `ENST00001008.6`). The patched R script strips version numbers from both sides for 100% match rate.
14. **FastQC/TrimGalore deadlock with Singularity**: When running many samples (>4) with Singularity, FastQC JVMs started by TrimGalore hang in `futex_do_wait` because `/tmp/hsperfdata_*` files collide across containers (Singularity shares host `/tmp` by default, unlike Docker). **Fix**: (a) `JAVA_TOOL_OPTIONS=-XX:-UsePerfData` is set in `nextflow.config` `env` block to disable JVM perfdata files; (b) `maxForks` is capped at 4 in `conf/modules.config` for TrimGalore/FastQC processes to limit concurrency; (c) FastQC threads are pinned to `-t 2` to prevent resource exhaustion. If running on a system with very slow shared `/tmp` (e.g. NFS), consider setting `TMPDIR` to a local disk as well.
15. **Nextflow 26.x strict schema validation**: CLI-passed params are parsed as strings. Boolean params must use `"type": ["boolean", "string"]` with `"enum": [true, false, "true", "false"]` in `nextflow_schema.json`. The old `nf-validation` `validation {}` config block is incompatible with `nf-schema@2.3.0+` and must be removed.
16. **Dual-genome `--additional_fasta` and RiboseQC**: The auto-generated GTF for `--additional_fasta` only contains `exon` features which breaks RiboseQC's TxDb construction (`Error: subscript contains invalid names`). **Fix**: pre-concatenate host + pathogen genomes and GTFs manually, then pass the combined files via `--fasta` and `--gtf` directly. Use `--pathogen_gtf` for pathogen-specific ORF prediction.
17. **Session lock on `-resume` after force-kill**: If Nextflow is killed abruptly (`kill -9`), the session cache at `.nextflow/cache/<uuid>/db/LOCK` persists and blocks the next `-resume`. **Fix**: `rm -rf .nextflow/cache/<session-uuid>/`.

## File Locations Reference

- Custom scripts: `bin/` (e.g., `filter_gtf.py`, `gtf2bed`)
- ORF unification: `scripts/unify_orf_predictions.py`
- ORF classification wrapper: `scripts/classify_orfs_wrapper.py`
- ORFquant classification library: `scripts/class_orf/orfquant_orf_classify.R`
- GENCODE mapper: `scripts/gencode-riboseqORFs/ORF_mapper_to_GENCODE_v1.1.py` + `functions.py`
- ORF-type classifier: `scripts/class_orf/class_ORFtype.py`
- Helper test scripts: `scripts/singularity_single_tool_tests/`
- Pipeline tests: `tests/` and `*.nf.test` files throughout
- Example data: `example/` directory
- Documentation: `docs/` (usage.md, output.md)
