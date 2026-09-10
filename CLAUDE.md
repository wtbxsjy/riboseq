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
13. **Transcript ID version mismatch**: GTF IDs (e.g. `ENST00001008.11`) have different version suffixes than BAM IDs (e.g. `ENST00001008.6`). The patched R script strips version numbers from both sides for 100% match rate. **Critical fix (2026-06-07)**: `renameSeqlevels()` cannot be used for version stripping because multiple isoforms (e.g. `AT1G01020.1`, `AT1G01020.2`) collapse to the same gene-level ID, violating `seqlevels` uniqueness. Instead, modify `seqnames(ga)` directly at the read level via character vector, then set `seqlevels(ga) <- unique(new_seqnames)`. The annotation table also deduplicates by gene-level ID (`!duplicated(transcript)`), keeping the first isoform's UTR/CDS lengths per gene — acceptable for P-site offset determination.
14. **FastQC/TrimGalore deadlock with Singularity**: When running many samples (>4) with Singularity, FastQC JVMs started by TrimGalore hang in `futex_do_wait` because `/tmp/hsperfdata_*` files collide across containers (Singularity shares host `/tmp` by default, unlike Docker). **Fix**: (a) `JAVA_TOOL_OPTIONS=-XX:-UsePerfData` is set in `nextflow.config` `env` block to disable JVM perfdata files; (b) `maxForks` is capped at 4 in `conf/modules.config` for TrimGalore/FastQC processes to limit concurrency; (c) FastQC threads are pinned to `-t 2` to prevent resource exhaustion. If running on a system with very slow shared `/tmp` (e.g. NFS), consider setting `TMPDIR` to a local disk as well.
15. **Nextflow 26.x strict schema validation**: CLI-passed params are parsed as strings. Boolean params must use `"type": ["boolean", "string"]` with `"enum": [true, false, "true", "false"]` in `nextflow_schema.json`. The old `nf-validation` `validation {}` config block is incompatible with `nf-schema@2.3.0+` and must be removed.
16. **Dual-genome `--additional_fasta` and RiboseQC**: The auto-generated GTF for `--additional_fasta` only contains `exon` features which breaks RiboseQC's TxDb construction (`Error: subscript contains invalid names`). **Fix**: pre-concatenate host + pathogen FASTA manually and pass via `--fasta`; keep `--gtf` host-only so host ORF tools/RiboseQC stay clean.
17. **Pathogen side is TE/known-gene quantification only** (design decision 2026-08-26): pathogen BAMs go through `SPLIT_BAM_BY_CONTIG` → `SORF_BAM_FILTER_PATHOGEN` → `TE_ANALYSIS_PATHOGEN` (GTF2BED → featureCounts counts → optional DESeq2 deltaTE). No de novo ORF prediction/unify/classify exists on the pathogen side; do not "restore" the removed `ch_pathogen_fasta_gtf` channel.
18. **Pathogen GTF/GFF3 format requirements**: `bin/gtf2bed` (used by `GTF2BED_PATHOGEN`) groups by `transcript_id` (GTF) or `ID` (GFF3, unquoted) and uses `exon` blocks with a `CDS` fallback (added 2026-08-26) — Prokka GFF3 (CDS-only) now works, RefSeq GTF works, but a GTF/GFF3 with neither exon nor CDS yields an empty BED and featureCounts fails. Seqnames must exactly match contigs in the combined `--fasta`.
19. **TE counts-only mode**: `TE_ANALYSIS` (host, lncRNA, and pathogen) runs without `--contrasts` (empty contrast channel → `DESEQ2_DELTATE` never executes; `MERGE_COUNTS` still produces the counts matrix). deltaTE also auto-skips (warn) when the samplesheet has no `rnaseq` samples. `--skip_te_analysis_pathogen` must be defined in `nextflow_schema.json` (boolean+string enum per gotcha 15).
20. **SPLIT_BAM_BY_CONTIG zero-match is a hard error**: if `--pathogen_contig_pattern` (POSIX ERE against the FAI first column) matches no contigs, the process exits 1 (previously warned and passed everything through as host). Check that pathogen sequences are concatenated into `--fasta`.
21. **Force-kill leaves the LOCK file but does NOT block `-resume`**: After `kill -9`, `.nextflow/cache/<uuid>/db/LOCK` persists, but NF 26.04.3's RocksDB recovers on reopen — verified 2026-09-02 23:52: PRJEB26593 was SIGKILLed at 23:37 (OOM storm) and `-resume` at 23:52 hit thousands of `Cached process` records, continuing normally. **Do NOT delete `.nextflow/cache/<session-uuid>/` after a force-kill** — that dir holds the resume DB, and deleting it is what makes the next `-resume` find no task records and **re-execute the entire DAG from scratch** (this is exactly what happened in round 3: 24/24 FQ_LINT + all STAR re-submitted after a cache deletion). After any kill, just relaunch with `-resume` and check the first minute of `.nextflow.log` for `Cached process` lines; only touch the cache if `-resume` actually errors.
22. **`.collect()` on an empty channel emits NOTHING** (NF 26): `channel.collect().map { it ?: [] }` never fires on empty input because `.collect()` silently emits nothing — the `.map` guard is dead code, the downstream `combine` chain starves, and the pipeline still reports SUCCESS. **Fix**: `.collect().ifEmpty([])` (idiom already used by COLLECT_QC_STATS wiring). This was why UNIFY_ORF_PREDICTIONS + all three classifiers + host TE_ANALYSIS silently never ran.
    **Root cause of the empty psites channels (verified 2026-09-04 with task logs + container source)**: RiboseQC `RiboseQC_analysis()` with the default `readlength_choice_method="max_coverage"` selects the max-coverage read length; if its frame_preference < 50% (choose_readlengths gate), chosen_rl becomes empty and the success path prints "Not enough signal or low frame preference, skipped P-sites & codon occupancy calculation" and writes NO `*_P_sites_*.bedgraph` (placeholders only on failure). PRJEB26593: all 12 Ribo-seq samples hit this (RL29 chosen, 47.1% < 50%). Samples that pass the gate DO emit bedgraphs (GSE208041 SRR20106164: RL28 86.5% → 841K-line P_sites bedgraphs). **Fix**: new process `RIBOSEQC_PSITES_ALL` (`modules/local/riboseqc/psites_all/`) runs `readlength_choice_method="all"` — P-sites positions computed for every read length using each length's cutoff from the P_sites_calcs table (better than `rescue_all_rls=TRUE`, which fills rescued lengths with hardcoded cutoff 12). UNIFY (and per-tool unify) now consume `RIBOSEQC_POSTFILTER.out.psites_bedgraph_all`. RIBOSEQC_ANALYSIS is deliberately untouched so its tasks — and the whole ORFquant chain — stay cache-valid under `-resume` (a module edit there would have re-run 12× ORFquant ≈ 10-21h for byte-identical outputs; P_sites_calcs is byte-identical between the modes since `choice` only affects final_choice selection). The header-only expression fallback (commit 2cbdc61) remains as a safety net for genuinely empty cases. Related NF fact: a `path "*_P_sites_*.bedgraph"` output glob matching 2 files emits ONE tuple per task with a file LIST (meta.id not duplicated in sample lists).
23. **Queue×queue process-input pairing is racy when one channel is combine-derived** (NF 26.04, verified with minimal repros): `DESEQ2_DELTATE(ch_contrasts(3 items), ch_merged_data(1 item))` silently formed only the FIRST contrast's task (real run: only `infection_2h_vs_0h`). **Fix**: `.first()` on the single-emission channel converts it to a value channel that broadcasts to every contrast (3/3 tasks in repeated runs, even when a sibling TE consumer starves). Use `.first()` for any "single matrix/annotation × per-item channel" pairing.
24. **`-stub-run` crashes on PREPARE_GENOME:BOWTIE_BUILD**: its stub writes an empty `versions.yml`, and `CUSTOM_DUMPSOFTWAREVERSIONS` dies with `collectEntries() on null object` as soon as that channel emits — unrelated to your edits, but it prevents stub-run from validating late-stage wiring. Stub-run reaches it after ~70 fast tasks; treat that NPE as expected.
25. **UNIFY_ORF_PREDICTIONS hits the 8h `process_medium` time limit with real P_sites bedgraphs** (verified 2026-09-05, PRJEB26593): 1.2M ORF candidates × 96 bedgraph files (1.5 GB) took >8h at `--threads 2` (was hardcoded in `modules/local/unify_orf_predictions/main.nf`); the task died at exactly task-start+8h (exit 143). Worse, **NF 26.04.3's local executor aborted the WHOLE session** with `IllegalThreadStateException: process hasn't exited` while handling the timeout kill (race in `LocalTaskHandler` polling a just-SIGTERMed process) instead of just failing the task. Fix (commit 9472fe4): `--threads ${task.cpus}` + `withName: 'UNIFY_ORF_PREDICTIONS|UNIFY_ORF_PREDICTIONS_PER_TOOL' { cpus = 16; time = { 48.h * task.attempt } }`. Note the internal caps in `unify_orf_predictions.py`: stats streaming is GIL-bound ThreadPool capped at 4 workers; sequence/CDS annotation uses `args.threads` via multiprocessing Pool (the phase that was killed at 37.5% through).
26. **DESeq2 deltaTE: fit the interaction model per contrast, and keep the heatmap annotations in sync** (verified 2026-09-05, PRJEB26593, commits a5bcfd8 + e118a54). (a) With all timepoints in one combined model (`~ contrast_var + seq_type + contrast_var:seq_type`), there are 3 interaction coefficients and the template always tested `interaction_coef[1]` (the first non-reference group, 2h) for EVERY contrast — all 3 contrasts produced IDENTICAL DTEG lists. Fix: subset `count_table`/`sample_sheet` to the contrast's two groups (`contrast_levels <- c(reference_level, target_level)`) before `DESeqDataSetFromMatrix`, so the single interaction coefficient tests this contrast. (b) Follow-on bug: `plot_heatmap` built `HeatmapAnnotation` from the FULL sample sheet (24 samples) while `dds_combined`/`mat_scaled` now had only the contrast's 12 columns → `number of observations in top annotation should be as same as ncol of the matrix`. Fix: derive annotations from `sample_sheet[colnames(mat_scaled), , drop = FALSE]`. Any per-contrast subset of `dds_combined` must be mirrored in every consumer of `sample_sheet`. Also: in the `DESEQ2_DELTATE` publishDir `saveAs` closure, `"plots/\$filename"` renders a LITERAL `plots/$filename` path (all plot files overwrite one another); use `"plots/${filename}"` (commit a5bcfd8, conf/modules.config:971).
27. **Post-UNIFY tail: hash determinism, GENCODE mapper tmp pollution, and host-OOM storms** (verified 2026-09-06, PRJEB26593, commit 08c4c9b). Three interacting failure modes, all fixed at once: (a) **Multi-file channel lists make task hashes order-sensitive**: UNIFY/ORF_QC take `.collect()`ed lists (tool files, 96 P-site bedgraphs, sample names) whose order is nondeterministic per session → different hash every resume → the whole tail re-ran each time (outputs were md5-identical, so content was never the trigger). Fix: `.map { files -> (files instanceof List ? files : [files]).sort { a, b -> a.name <=> b.name } }` after every `.collect()` feeding these tasks. (b) **GENCODE mapper writes tmp INTO the ensembl dir**: `ORF_mapper_to_GENCODE_v1.1.py:148` does `os.mkdir(folder + "/tmp/")` and dumps ~250-500 GB of `.ov`/`.bed` files there. With symlink staging the SOURCE `gencode_v49/` mutated → its dir-input hash changed after every run → the classifier could never cache-hit, and 732 GB piled up (checked: du). Fix: `stageInMode = 'copy'` in the `CLASSIFY_ORFS_GENCODE` withName block (tmp then lands in the task workdir); delete accumulated `gencode_v49/tmp/*.ov *.bed` from the reference dir. (c) **Host OOM storms on the shared 188G box** (~95 GB chronically used by other tenants): after UNIFY finished, ORF_QC (peak_rss 90 GB! — compare_orf_tools.py on the 3.3 GB metadata.tsv) + GENCODE + ORFquant/ORF-type + EXPRESSION_QUANT ran concurrently → kernel SIGKILLed compare_orf_tools.py ×2 and the GENCODE mapper (exit 137, "died with <Signals.SIGKILL: 9>"). No dmesg access → diagnose via exit-137 cluster + `free -g`. Fix: per-run `executor.queueSize = 1` (executor_serial.config in the run's process/ dir, not the repo — the local executor otherwise submits everything runnable at once), OOM-retry `errorStrategy = { task.exitStatus in 137..143 ? 'retry' : 'terminate' }` + `maxRetries 2` for ORF_QC/classifiers, and honest memory budgets (`ORF_QC` 96 GB measured peak; note NF's local executor does NOT enforce `memory` — it's a hint). Split long tails into phases: phase 1 `--skip_te_analysis true` (serial heavy tail), phase 2 resume with `queueSize = 8` for the cheap TE featureCounts ×24. **Watch out**: killing a hung NF daemon does NOT kill its submitted task trees (they reparent to init and keep running — a 66 GB compare_orf_tools.py survived its daemon by hours); kill the task PIDs explicitly.
28. **`{sample}_pN` in the expression summary is structurally ≥ 1 and cannot be used as a rice-style threshold** (verified 2026-09-07, PRJEB26593): pN = `max_bg_value × n_intervals / sum_bg_values`, and since max ≥ sum/n_intervals, pN ≥ 1 whenever the ORF has any reads (10.1M reads>9 records: min 1.0, 0% below 1). The rice/maize post_analysis threshold `pN > 0.5` is vacuous under this definition. For a Stage-1 filter use `{sample}_reads > 9` (≥1 sample), and put the "≥50% ratio" requirement on the P-site purity metric instead. **Two caveats about `post_analysis/scripts/compute_psite_purity.py` when using it for that ratio** (verified 2026-09-08, PRJEB26593): (a) its per-sample `{sample}_p_site_pos` column is buggy — it divided the value-sum (`_p_site_pos_wt`, key fixed now) by psite so it was always 1.0; use the ORF-level `global_p_site_pos` (value×pos weighted, correct) as the position criterion. (b) Its `p_site_pct` = raw P-site counts / RiboseQC **coverage** bedgraph values is scale-invalid: the coverage bedgraphs are RPM-normalized (genome-wide integral ≈ 1e6, verified) while the P_sites bedgraphs are raw counts. Compute the ratio on raw counts instead: `pct = {sample}_p_site_GSE (purity) / {sample}_reads (expression summary, raw psite+psite_uniq)`. With this definition (PRJEB26593, stage1 set): ratio q50 = 0.49, and `psite > 9 AND pct ≥ 0.5` in ≥1 sample is the working Stage-2 filter.


29. **PRICE's contig names lack the `chr` prefix → 26.9% of the unified ORFs are silently zeroed out** (verified 2026-09-10, PRJEB26593; fix uncommitted at time of writing). PRICE (GEDI) emits Ensembl-style contig names (`1` … `X`, `Y`, `MT`) in `{sample}.orfs.tsv` / `{sample}_Detected_ORFs.gtf`, while the reference FASTA/BAM and the other four tools (Ribo-TISH, Ribotricer, ORFquant, RiboCode) use `chr`-prefixed names. Consequences for the 327,119 PRICE ORFs (= 26.9% of 1,215,479 unified ORFs in PRJEB26593): (a) `QUANTIFY_ORFS`/featureCounts can never match a BAM contig → counts identically 0 (all 10,641 ORFs with any count in the published `merged_counts.tsv` are chr-prefixed); (b) the UNIFY P-site quantification is zeroed the same way — `total_psites > 0` for **0.00%** of PRICE ORFs vs 99.98% for other tools; (c) `extract_sequence()` does `genome_fasta[cand.chrom]` → PRICE ORFs get an all-`N` `sequence` column (unusable for AMP analysis); (d) exact-match merging with an identical ORF called by another tool never fires (`id_key` contains `chrom`). Proof it is naming, not missing signal: after adding the `chr` prefix, 55.1% of PRICE ORFs have ≥1 P-site in a single sample (other tools: 63.5%). **Fix**: `scripts/unify_orf_predictions.py` gained `_make_chrom_normalizer(gtf_index)` (maps onto the reference's own naming via `GTFIndex.chrom_names` + the existing `_chrom_aliases()`; idempotent; GL/KI scaffolds are left unchanged) and main() now normalises every candidate after parsing and **before** merging / GTF lookups / sequence extraction, rebuilding `cand.id_key` (that key caches `chrom` in `__init__`, so a stale key would corrupt dedup). Rerunning UNIFY is mandatory after the fix and it **renumbers all ORF IDs** — post_analysis pass lists must be re-mapped by coordinates. Whenever a new predictor is added, check its contig convention against the reference.
30. **featureCounts unique assignment throws away ~88% of reads on the redundant unified-ORF annotation — use `-O` for ORF-level counting** (verified 2026-09-10, PRJEB26593). The unified ORF set is massively overlapping (1,215,479 ORFs fall into 20,612 overlap clusters; only 0.46% are isolated, 99.5% share a cluster with ≥2 ORFs), so `QUANTIFY_ORFS`'s default counting marks **6,074,672 of 6,880,139 reads (88.3%) as `Unassigned_Ambiguity`** on ERR2603016 (assigned 259,683, matches the published matrix exactly) and only 1,664 of the 13,552 AMP-scored ORFs (12.3%) receive any count — the "quantified" subset is decided by locus crowding, not by translation. **Fix for ORF-level DE (post_analysis, not yet in the pipeline)**: add `-O` (count a read for every feature it overlaps). This makes each count equal the number of reads overlapping that ORF's span — verified against an independent `samtools` region count (4,130 = 4,130) — and restores coverage to 99.7–100% of the QC-passed sets. Caveats: absolute values are inflated (≈27× vs `--fraction`) and overlap-sharing makes neighbouring ORFs' counts non-independent, so treat DE of ORFs inside dense clusters with care; a `-O --fraction` sensitivity run gives per-ORF log2FC r = 0.948 (Spearman 0.963, sign agreement 89.6%, only 4.2% of ORFs differ by >0.5). DESeq2 (this version) rejects non-integer matrices, so fractional counts cannot be fed to `DESeqDataSetFromMatrix` directly. Working implementation: `post_analysis/PRJEB26593/de_quant/run_orf_de.R` (same BAMs/annotation as `QUANTIFY_ORFS`, `sfType="poscounts"`, per-contrast sample subsetting, TE prefilter rule `max(te_prefilter_min_nonzero, ceil(frac × N))`).

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

## Root Directory File Organization

Keep the root directory lean — only standard project-level files belong there. All other files should be placed in appropriate subdirectories.

### Files that stay in root

- `README.md` — project homepage
- `CHANGELOG.md` — version changelog
- `CLAUDE.md` — Claude Code project instructions
- `AGENTS.md` — AI agent instructions
- `CITATIONS.md` — citation information
- `CODE_OF_CONDUCT.md` — community code of conduct
- `nextflow.config` — pipeline configuration
- `main.nf` — pipeline entry point

### Shell scripts → `scripts/`

Test-runner shell scripts (e.g. `run_test_codespace.sh`, `run_test_colab.sh`) go in `scripts/`, not root.

### Development logs & planning docs → `docs/devlog/`

One-off planning documents, feasibility reports, session summaries, implementation checklists, sync plans, and temporary changelogs are archived in `docs/devlog/`. This includes files like:

- `GENCODE_*_SUMMARY.md`, `GENCODE_*_FEASIBILITY.md` — integration planning
- `IMPLEMENTATION_CHECKLIST.md`, `SYNC_PLAN_*.md` — implementation notes
- `TECHNICAL_REFERENCE.md`, `SUMMARY.md` — technical references
- `CHANGELOG_*.md` (date-stamped) — temporary changelogs between releases

### Tutorial/demo data → tracked in `test_data/`

`test_data/` contains bundled demo data used by `docs/notebooks/` tutorials. It is **tracked in git** (not gitignored) so that notebooks are self-contained and reproducible for users who clone the repo.
