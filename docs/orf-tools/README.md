# ORF Prediction Tools — Analysis Documentation

This directory contains detailed analysis of all ORF prediction tools integrated in the riboseq pipeline, and a systematic design for a unified ORF QC module.

## Tool Analysis Documents

| # | Tool | Type | Key QC Contribution |
|---|------|------|---------------------|
| [01](01-ribocode.md) | **RiboCode** | De novo ORF detection | Per-ORF p-values, frame-0/1/2 coverage, f0% periodicity |
| [02](02-riboseqc.md) | **RiboseQC** | Comprehensive QC | P-site offsets, frame distributions, codon occupancy, biotype QC |
| [03](03-orfquant.md) | **ORFquant** | Splice-aware ORF quantification | Multi-level classification, DPSS periodicity, isoform-aware |
| [04](04-ribotricer.md) | **Ribotricer** | Reference-guided ORF detection | Fourier coherence phase score, cross-correlation P-site |
| [05](05-ribotish.md) | **Ribo-TISH** | De novo TIS/ORF prediction | TIS enrichment p-values, frame Q-values, TI-seq integration |
| [06](06-price.md) | **PRICE (rp-bp)** | Bayesian ORF prediction | Bayes factors (periodic/translated), GP model comparison |
| [02b](02b-ribowaltz.md) | **riboWaltz** | P-site offset & diagnostic QC | Two-step coherence correction P-site, frame distributions, metagene profiles |

## QC Module Design

| # | Document | Content |
|---|----------|---------|
| [07](07-orf-qc-module-design.md) | **ORF QC Module Design** | Complete design: QC dimensions, harmonized metrics, OCS scoring, architecture, implementation plan |

## Quick Reference: QC Metrics by Dimension

### Read-Level Periodicity & P-site
| Metric | RiboCode | RiboseQC | riboWaltz | Ribotricer | Ribo-TISH | PRICE |
|--------|----------|----------|-----------|------------|-----------|-------|
| P-site offset | ✅ | ✅ | ✅ (best algorithm) | ✅ | ✅ | ✅ |
| 3-nt periodicity score | `f0_percent` | frame distributions | frame 0/1/2 % per region | `phase_score` | metagene PDF | `BF` (periodic) |

### ORF-Level Translation Evidence
| Metric | RiboCode | RiboseQC | Ribotricer | Ribo-TISH | PRICE |
|--------|----------|----------|------------|-----------|-------|
| Significance p-value | ✅ `adjusted_pval` | — | ✅ (phase_score p) | ✅ `FisherQvalue` | — |
| Bayes factor | — | — | — | — | ✅ |
| Phase/Frame score | — | ✅ (frame %) | ✅ `phase_score` | ✅ `FrameQvalue` | ✅ (BF) |
| Coverage completeness | ✅ | — | ✅ `valid_codons_ratio` | — | — |
| ORF abundance | ✅ `RPKM` | ✅ `RPKM` | ✅ `read_density` | — | ✅ (profile sum) |

### ORF Classification
| Category | RiboCode | Ribotricer | Ribo-TISH | ORFquant |
|----------|----------|------------|-----------|----------|
| CDS/annotated | ✅ | ✅ | ✅ | ✅ |
| uORF | ✅ | ✅ | ✅ | ✅ |
| dORF | ✅ | ✅ | ✅ | ✅ |
| Overlap types | ✅ | — | — | ✅ |
| Internal | ✅ | — | — | — |
| Novel/ncRNA | ✅ | ✅ | ✅ | ✅ |
