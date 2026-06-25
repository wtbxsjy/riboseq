# Post-Analysis Integrated ORF Scoring — Development Record

> **Date**: 2026-06-25 | **Branch**: `dev` | **Status**: Implementation complete, awaiting rice re-run for full validation

## Overview

Designed and implemented a unified 4-dimensional ORF confidence scoring framework that integrates three previously independent filtering standards into a single comprehensive post-analysis pipeline. The key insight driving this work: the pipeline already produces most metrics needed for quality assessment — the post_analysis scripts were redundantly re-computing them from raw bedgraph files.

## Background: Three Independent Systems

Prior to this work, three separate ORF quality assessment approaches existed in parallel:

| System | Core Dimension | Status Before | 
|--------|---------------|---------------|
| **Yuanliang rules** (`yuanliang_rules.md`) | P-site signal purity | Documentation only, unused |
| **post_analysis** (`plot_sorf_strict_batch.R`) | Expression quality + reproducibility | Active, but re-computed from bedgraph |
| **ORF_QC module** (`compare_orf_tools.py`) | Multi-tool cross-validation confidence | Active in pipeline |

### Detailed Comparison

See `design_integrated_scoring.md` for the full 3-way comparison with threshold matrices.

Key metrics that differ across systems:

- **Yuanliang**: `p_site_GSE` (P-site reads only), `p_site_percentage` (P-site/total), `p_site_postion_GSE` (position offset). Requires distinguishing P-site vs non-P-site reads at the exon level.
- **post_analysis**: `total_reads` (all reads), `max_pn` (peak/mean ratio per sample), `avg_ratio` (mean coverage per position), `n_samples`. Does NOT distinguish P-site from non-P-site.
- **ORF_QC**: OCS = weighted combination of S_translation (tool significance), S_agreement (cross-tool overlap), S_coverage (P-sites/codon), S_periodicity (pN density), S_readlevel (global quality).

## Design: 4-Tier Scoring Framework

### Sub-Scores

```
S_purity (0.25)      — "Is the P-site signal real?"         Yuanliang metrics
S_expression (0.30)  — "Is translation strong + reproducible?" post_analysis metrics
S_confidence (0.25)  — "Do multiple tools agree?"            ORF_QC OCS
S_biology (0.20)     — "Is this ORF biologically interesting?" GENCODE class
```

### Composite Score

```
Composite = 0.25 × S_purity + 0.30 × S_expression + 0.25 × S_confidence + 0.20 × S_biology
```

### Tiers

| Tier | Score | Interpretation |
|------|-------|----------------|
| **Gold** | ≥ 0.70 | High-confidence, suitable for experimental validation |
| **Silver** | 0.50–0.69 | Good candidate, minor weaknesses acceptable |
| **Bronze** | 0.30–0.49 | Plausible but uncertain |
| **Weak** | < 0.30 | Likely noise or very low expression |

### Key Design Decisions

1. **S_expression gets highest weight (0.30)** — strong, reproducible translation signal is the most direct physical evidence
2. **S_purity = 0.25** — P-site quality is essential to rule out degradation artifacts
3. **S_confidence = 0.25** — multi-tool agreement provides orthogonal statistical validation
4. **S_biology = 0.20** — biological context informs but should not dominate physical evidence

## Implementation

### New Files

| File | Lines | Purpose |
|------|-------|---------|
| `config_thresholds.yaml` | 61 | All thresholds, configurable per organism |
| `scripts/compute_psite_purity.py` | 296 | Stream P-site + coverage bedgraphs to compute Yuanliang metrics per ORF per sample |
| `scripts/integrate_orf_scores.R` | 579 | Merge all pipeline TSVs, compute 4 sub-scores + composite + legacy flags |
| `design_integrated_scoring.md` | 380 | Full design document with data flow diagrams, formula derivations, migration plan |

### Data Flow

```
PIPELINE OUTPUTS (read directly, zero redundant computation)
═══════════════════════════════════════════════════════════
metadata.tsv            ─┐
expression_summary.tsv   ┤
expression_rpkm_tpm.tsv  ┤
orf_confidence.tsv       ┤──→ integrate_orf_scores.R  →  integrated_orf_scores.tsv
classification.tsv       ┤
                         ┤
P-site bedgraphs ────────┤──→ compute_psite_purity.py  →  psite_purity.tsv
coverage bedgraphs ──────┘    (only new computation needed)
```

### What Was Eliminated

The existing post_analysis pipeline re-computed these metrics from bedgraph files:

- `extract_amp_expression_omp.py` — per-sample `_reads`, `_pN`, `_mean_cov` → **redundant**: pipeline's `expression_summary.tsv` already has `{sample}_reads` and `{sample}_pN`
- `calc_orf_rpkm_omp.py` — per-sample RPKM/TPM → **redundant**: pipeline's `expression_rpkm_tpm.tsv` already has these

Estimated savings: ~80% of post_analysis runtime (~15-45 min → <1 min + purity computation).

### What Still Needs Bedgraph Scanning

Only Yuanliang's P-site purity metrics require new computation from bedgraph files:

- `p_site_GSE` per sample — sum of P-site bedgraph values overlapping ORF
- `p_site_percentage` per sample — `p_site_GSE / reads_GSE`
- `p_site_postion_GSE` — weighted average P-site position relative to ORF start

**Long-term optimization**: Add `{sample}_psites` to pipeline's `expression_summary.tsv` (~3 lines in `unify_orf_predictions.py`). Once done, `compute_psite_purity.py` becomes completely obsolete — all metrics come directly from pipeline TSVs.

### Legacy Compatibility

`strict_sorf_pass` flag preserved with identical logic to `plot_sorf_strict_batch.R`:
- `orf_biotype ∈ {uORF, dORF, doORF, uoORF, intORF, lncRNA}`
- `orf_len_aa ≤ 50`
- `total_reads ≥ 20`
- `max_pn ≥ 5`
- `avg_ratio ≥ 2`
- `n_samples ≥ 3`

### Test Results (Rice, partial data)

Tested with 16,668 ORFs from rice post_analysis data (using expression_summary as metadata — a worst-case fallback since `length_aa` must be approximated from `end-start+1` which over-counts multi-exon ORFs):

| Tier | Count | % |
|------|-------|---|
| Gold | 60 | 0.4% |
| Silver | 7,536 | 45.2% |
| Bronze | 8,765 | 52.6% |
| Weak | 307 | 1.8% |
| **strict_sorf_pass** | 5 | — |

> **Note**: These numbers are with incomplete data (no pipeline metadata.tsv with accurate `length_aa`, no orf_confidence.tsv). The Gold count is artificially low because S_confidence is NA (no OCS data). Real numbers expected after rice re-run.

## Pending Validation

- [ ] Rice pipeline re-run with all fixes → obtain complete `metadata.tsv`, `orf_confidence.tsv`, `expression_summary.tsv`, `expression_rpkm_tpm.tsv`
- [ ] Run `compute_psite_purity.py` on rice RiboseQC outputs
- [ ] Run `integrate_orf_scores.R` with full pipeline data
- [ ] Cross-validate: compare new Gold-tier ORFs with existing `rice_sorf_strict_list.tsv`
- [ ] Tune thresholds based on biological ground truth
- [ ] Repeat for maize (100 samples)

## Related Commits

Branch `dev-rice_run` → `dev` (22 commits, fast-forward merge, 2026-06-25):
- ORF_QC: channel single-consumer fix, real OCS output
- ORFquant: Bioc 3.20 compat, BSgenome forge, FaFile parallelism
- PRICE: absolute `-genomic` path fix
- RiboCode: `_open()` gz fallback
- EXPRESSION_QUANT: merged into UNIFY streaming (27h→<10s)
- MULTIQC: PoisonPill filter

## Usage

```bash
# Step 1: Compute P-site purity (Yuanliang metrics)
python3 scripts/compute_psite_purity.py \
  --bed results/unify/unified_orfs.bed \
  --riboseqc-dir results/riboseqc \
  --output psite_purity.tsv \
  --workers 8

# Step 2: Integrated scoring
Rscript scripts/integrate_orf_scores.R \
  --metadata results/unify/unified_orfs.metadata.tsv \
  --expression results/unify/expression_summary.tsv \
  --rpkm results/unify/expression_rpkm_tpm.tsv \
  --confidence results/orf_qc/orf_confidence.tsv \
  --purity psite_purity.tsv \
  --classification results/classification/result_with_classification.tsv \
  --config config_thresholds.yaml \
  --output integrated_orf_scores.tsv
```
