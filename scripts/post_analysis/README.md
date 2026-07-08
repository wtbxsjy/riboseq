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
| `riboseqc/*_P_sites_{plus,minus}.bedgraph` | RiboseQC | P-site density (optional) |

## Output

```
{output_dir}/
├── 01_filtering_analysis.html    # Cross-analysis report
├── 02_generate_ggribo.html       # ggRibo coverage plots
├── final_orfs_for_ggribo.tsv     # Selected ORF list
├── ggribo_plots/                 # Per-biotype ggRibo PNGs
│   ├── index.html
│   ├── uORF/
│   ├── dORF/
│   └── ...
└── logs/
```

## Steps

1. **Expression-based filtering** — Uses `_reads` and `_pN` columns from
   expression summary as proxies for P-site quality.
2. **Cross-analysis** — ORF biotype × OCS tier × tool agreement × expression filters.
3. **ggRibo coverage plots** — Per-ORF ribosome footprint density with
   reading-frame shading.

## Adding a New Project

```bash
mkdir -p run/NEW_PROJECT/post_analysis/output/{logs,ggribo_plots}
cp config_template.yaml run/NEW_PROJECT/post_analysis/project_config.yaml
# Edit paths in project_config.yaml
bash run_all.sh run/NEW_PROJECT/post_analysis/project_config.yaml
```
