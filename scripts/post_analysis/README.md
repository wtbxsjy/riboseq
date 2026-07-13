# ORF Post-Analysis Module

Reusable analysis pipeline for nf-core/riboseq ORF outputs.

## Quick Start

```bash
# 1. Copy config template to your project
cp scripts/post_analysis/config_template.yaml run/YOUR_PROJECT/post_analysis/project_config.yaml

# 2. Edit the config with your paths

# 3. Run
bash scripts/post_analysis/run_all.sh run/YOUR_PROJECT/post_analysis/project_config.yaml
```

## Input Requirements

All files are standard pipeline outputs:

| File | Source Module | Description |
|------|--------------|-------------|
| `unified_orfs.metadata.tsv` | UNIFY_ORF_PREDICTIONS | ORF coordinates + metadata |
| `unified_orfs_orf_confidence.tsv` | ORF_QC | OCS scores + tier |
| `unified_orfs_expression_summary.tsv` | EXPRESSION_QUANT | Per-sample reads + pN |
| `gencode_results.orfs.out.gz` | CLASSIFY_ORFS_GENCODE | ORF biotype classification |
| `unified_orfs.bed.gz` | UNIFY_ORF_PREDICTIONS | BED12 ORF coordinates |
| `riboseqc/*_P_sites_{plus,minus}.bedgraph` | RiboseQC | P-site density |

## Two-Stage Filtering Pipeline

The pipeline uses a **two-stage** strategy to balance accuracy and performance:

### Stage 1: Preliminary (expression-based, fast)
`01_prelim_analysis.qmd` — Uses `{sample}_reads` and `{sample}_pN` from the expression
summary to eliminate ORFs with insufficient signal. Reads correlate strongly with
real P-site counts (r=0.99); pN measures frame periodicity.

**Output:** `prelim_orfs_for_psite.bed` — reduced ORF set for P-site computation.

### Stage 2: Real P-site (bedgraph-based, exact)
`compute_psite_purity.py` — Backtracks through actual RiboseQC bedgraph files
(`_P_sites_{plus,minus}.bedgraph` and `_coverage_{plus,minus}.bedgraph`),
computing exact `p_site_GSE`, `p_site_pct`, and `p_site_pos` per ORF per sample.
Runs on the Stage 1 reduced set, dramatically cutting runtime.

`02_psite_filtering.qmd` — Applies real P-site thresholds and selects final candidates.

### Stage 3: ggRibo Visualization
`03_generate_ggribo.qmd` — Generates per-ORF ribosome footprint coverage plots
with reading-frame shading. Top N samples sorted by a configurable metric
(default: `p_site_GSE`).

## Output

```
{output_dir}/
├── 01_prelim_analysis.html       # Stage 1: preliminary analysis
├── prelim_orfs_for_psite.bed      # ORFs for P-site computation
├── psite_purity.tsv               # Real P-site purity data
├── 02_psite_filtering.html        # Stage 2: P-site filtering report
├── final_orfs_for_ggribo.tsv      # Final selected ORFs
├── 03_generate_ggribo.html        # ggRibo plots + index
├── ggribo_plots/                  # Per-biotype ggRibo PNGs
│   ├── uORF/
│   ├── dORF/
│   ├── intergenic/
│   └── ...
├── tmp_gtf/                       # Temp GTF files (can delete)
└── logs/
```

## Intergenic ORFs

ORFs that don't overlap any annotated transcript in the GENCODE reference are
labelled **"intergenic"**. They are inherently non-CDS and are included in
downstream analysis alongside classified non-CDS ORFs (uORF, dORF, etc.).
This recovers the ~53% of unified ORFs that would otherwise be silently dropped.

## Adding a New Project

```bash
mkdir -p run/NEW_PROJECT/post_analysis/output/{logs,ggribo_plots}
cp config_template.yaml run/NEW_PROJECT/post_analysis/project_config.yaml
# Edit paths in project_config.yaml
bash run_all.sh run/NEW_PROJECT/post_analysis/project_config.yaml
```
