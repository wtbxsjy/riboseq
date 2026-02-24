# Copilot instructions for nf-core/riboseq

## Build, test, and lint
- Lint: `nf-core pipelines lint .`
- Full test suite: `nf-test test --profile debug,test,docker --verbose`
- Single test: `nf-test test modules/local/orfquant/main.nf.test --profile debug,test,docker`
- Pipeline tests: `nf-test test tests/pipeline/ --profile debug,test,docker`
- Smoke run: `nextflow run . -profile test,docker --outdir results`

## High-level architecture
- Nextflow DSL2 pipeline: entry point `main.nf`, core workflow in `workflows/riboseq/main.nf`.
- Subworkflows live in `subworkflows/local/` and `subworkflows/nf-core/`; processes in `modules/local/` and `modules/nf-core/`.
- Two input modes: FASTQ (full preprocessing/alignment/QC/ORF prediction) and BAM (skips preprocessing/alignment and requires explicit strandedness).
- Main flow (FASTQ): preprocessing → alignment (STAR/HISAT2) → sort/index → UMI dedup → sORF BAM filtering → QC (RiboseQC/Ribo-TISH) → ORF prediction (Ribo-TISH, Ribotricer, ORFquant; optional RiboCode, rp-bp) → MultiQC.
- RiboseQC produces inputs required by ORFquant; skipping RiboseQC auto-skips ORFquant.
- After per-sample ORF prediction, two post-processing stages run automatically:
  1. **ORF Unification** (`scripts/unify_orf_predictions.py`): merges Ribo-TISH / Ribotricer / ORFquant outputs into a non-redundant set; outputs `.bed` (BED12), `.gtf`, `.metadata.tsv`, `.stats.txt` under `orf_unification/`. Tool names in `sources` / metadata `tools` column are `'Ribo-TISH'`, `'Ribotricer'`, `'ORFquant'` (capitalised/hyphenated).
  2. **ORF Classification** (`scripts/classify_orfs_wrapper.py`): three parallel classifiers — `gencode` (transcriptome biotype, needs `--gencode_orf_mapper_container` with bedtools+BioPython and `--orf_classify_ensembl_dir`), `orfquant` (gene+transcript-level categories via transcript-space projection in `scripts/class_orf/orfquant_orf_classify.R`), and `orf_type` (gene-level). Outputs under `orf_classification/gencode/`, `orf_classification/orfquant/`, `orf_classification/orf_type/`.

## Key conventions
- Channel naming: `ch_output_from_<process>` for initial channels, `ch_<previousprocess>_for_<nextprocess>` for downstream channels.
- Parameters live in `nextflow.config` with defaults; run `nf-core pipelines schema build` to update `nextflow_schema.json` after changes.
- Resource defaults are defined by `withLabel:` selectors in `conf/base.config`; processes should pass `${task.cpus}`/`${task.memory}` to tools.
- sORF filtering produces filtered BAMs used by all ORF predictors; unfiltered BAMs are retained for baseline QC.
- Sample `type` drives behavior: only `type=riboseq` runs RiboseQC, and ORF prediction tools run per-sample only (no pooled mode by default).
- ORFquant runs with a patched container (`containers/Singularity.orfquant.patched.def`) to avoid R namespace conflicts.
