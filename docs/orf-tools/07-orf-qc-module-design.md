# ORF QC Module — Systematic Design (v2)

## 1. Design Rationale

The riboseq pipeline integrates **7 analysis tools** spanning two categories:

| Category | Tools | Role |
|----------|-------|------|
| **Pure QC** | RiboseQC, riboWaltz | P-site offset, periodicity, read distribution — no ORF prediction |
| **ORF Prediction** | RiboCode, ORFquant, Ribotricer, Ribo-TISH, PRICE | ORF detection + per-ORF QC metrics |

Currently:
- QC metrics are scattered across tool-specific output files with incompatible formats
- No cross-tool comparison or harmonization exists
- No aggregate ORF quality score is computed
- Users must manually inspect multiple output files to assess data quality
- **P-site offsets from different tools may disagree** — no consensus mechanism
- The **P-site fallback chain** (RiboseQC → riboWaltz → hardcoded defaults) is implicit, not documented in output

**Goal**: Design a unified ORF QC module that:
1. Collects and harmonizes QC metrics from all 7 tools
2. Computes cross-tool agreement metrics
3. Produces a single, comprehensive ORF QC report
4. Assigns per-ORF confidence scores based on multi-tool evidence
5. Documents the P-site offset consensus and fallback decisions

---

## 2. Tool Taxonomy & Data Availability Matrix

### 2.1 Which tools produce which QC data?

```
                    ┌──────────┬──────────┬───────────┬──────────┬──────────┬──────────┬───────┐
                    │RiboseQC  │riboWaltz │RiboCode   │ORFquant  │Ribotricer│Ribo-TISH │PRICE  │
┌───────────────────┼──────────┼──────────┼───────────┼──────────┼──────────┼──────────┼───────┤
│ P-site offset     │    ✅    │   ✅¹    │    ✅     │    ✗     │    ✅    │    ✅    │  ✅   │
│ per read length   │          │          │           │          │          │          │       │
├───────────────────┼──────────┼──────────┼───────────┼──────────┼──────────┼──────────┼───────┤
│ Read-length       │    ✅    │    ✅    │    ✅     │    ✗     │    ✅    │    ✅    │  ✅   │
│ distribution      │          │          │           │          │          │          │       │
├───────────────────┼──────────┼──────────┼───────────┼──────────┼──────────┼──────────┼───────┤
│ 3-nt periodicity  │    ✅²   │    ✅³   │    ✅⁴    │    ✅⁵   │    ✅⁶   │    ✅⁷   │  ✅⁸  │
│ (read-level)      │          │          │           │          │          │          │       │
├───────────────────┼──────────┼──────────┼───────────┼──────────┼──────────┼──────────┼───────┤
│ 3-nt periodicity  │    ✗     │    ✗     │    ✅⁴    │    ✅⁵   │    ✅⁶   │    ✅⁷   │  ✅⁸  │
│ (ORF-level)       │          │          │           │          │          │          │       │
├───────────────────┼──────────┼──────────┼───────────┼──────────┼──────────┼──────────┼───────┤
│ ORF abundance     │    ✗     │    ✗     │    ✅     │    ✅    │    ✅    │    ✗     │  ✅   │
│ (per ORF)         │          │          │           │          │          │          │       │
├───────────────────┼──────────┼──────────┼───────────┼──────────┼──────────┼──────────┼───────┤
│ Coverage          │    ✗     │    ✗     │    ✅     │    ✗     │    ✅    │    ✗     │  ✗    │
│ completeness      │          │          │           │          │          │          │       │
├───────────────────┼──────────┼──────────┼───────────┼──────────┼──────────┼──────────┼───────┤
│ ORF classification│    ✗     │    ✗     │    ✅     │    ✅    │    ✅    │    ✅    │  ✗    │
├───────────────────┼──────────┼──────────┼───────────┼──────────┼──────────┼──────────┼───────┤
│ Biotype /         │    ✅    │    ✗     │    ✗      │    ✗     │    ✗     │    ✗     │  ✗    │
│ contamination     │          │          │           │          │          │          │       │
├───────────────────┼──────────┼──────────┼───────────┼──────────┼──────────┼──────────┼───────┤
│ Codon occupancy   │    ✅    │    ✅    │    ✗      │    ✗     │    ✗     │    ✗     │  ✗    │
├───────────────────┼──────────┼──────────┼───────────┼──────────┼──────────┼──────────┼───────┤
│ Metagene profile  │    ✅    │    ✅    │    ✅     │    ✗     │    ✅    │    ✗     │  ✅   │
│ (visual)          │    (PDF)  │    (PDF)  │    (PDF)   │          │    (PDF)  │    (PDF)  │  (CSV)│
└───────────────────┴──────────┴──────────┴───────────┴──────────┴──────────┴──────────┴───────┘

¹ Two-step coherence correction — best algorithm among all tools
² frame distributions per transcript × read length (boxplots)
³ frame 0/1/2 % per region (5'UTR/CDS/3'UTR) — aggregated
⁴ Wilcoxon test f0>f1, f0>f2 → pval_combined
⁵ DPSS multitaper spectral analysis
⁶ Fourier coherence phase score
⁷ FrameQvalue (TIS-level frame preference)
⁸ Bayes factor: periodic vs non-periodic GP models
```

### 2.2 Tool Roles in QC Module

| Role | Tools | What they contribute |
|------|-------|---------------------|
| **P-site authority** | riboWaltz (best algorithm), RiboseQC (fallback) | Read-level P-site offsets |
| **Periodicity authority** | RiboseQC (per-transcript distributions), riboWaltz (per-region %) | Read-level 3-nt periodicity |
| **ORF evidence** | RiboCode, ORFquant, Ribotricer, Ribo-TISH, PRICE | Per-ORF translation significance |
| **Classification** | RiboCode, ORFquant, Ribotricer, Ribo-TISH | ORF type assignment |
| **Contamination** | RiboseQC | rRNA/mtDNA/plastid assessment |

---

## 3. QC Dimensions & Harmonized Metrics

### Dimension 1: Read-Level Periodicity & P-site Quality

Assesses whether raw Ribo-seq data shows the expected 3-nt periodic pattern.

| Harmonized Metric | riboWaltz (authority) | RiboseQC | RiboCode | Ribotricer | Ribo-TISH | PRICE |
|-------------------|----------------------|----------|----------|------------|-----------|-------|
| **P-site offset** | `corrected_offset_from_5` in `*_psite_offset.tsv` | `P_sites_calcs` | `predicted_psite` | `*_psite_offsets.txt` | `offdict` in `*.para.py` | `highest_peak_offset` |
| **Periodicity score** | Frame 0% in CDS (from `*_frame_distribution.tsv`) | Per-transcript frame distrib. | `f0_percent` (0–1) | Phase score (0–1) | Metagene PDF (qualitative) | `bayes_factor_mean` |
| **Read length distrib.** | `*_read_length_distribution.pdf` | `rld` counts | `length_counter` | `*_read_length_distribution.pdf` | `*_qual.txt` | `profile_sum` |
| **Optimal extremity** | Auto 5′ vs 3′ selection | — | — | — | — | — |

**P-site harmonization strategy:**

```
1. Primary reference: riboWaltz corrected_offset_from_5
   (two-step coherence correction, most accurate algorithm)
2. Fallback: RiboseQC P_sites_calcs
3. Cross-validation: flag lengths where any tool disagrees with reference by >1 nt
4. Consensus offset = median of all available tools (if ≥2 agree within 1 nt)
5. Pipeline fallback chain for downstream use:
   RiboseQC valid → use RiboseQC
   RiboseQC empty + riboWaltz available → use riboWaltz
   Neither → hardcoded default (28–32 nt → 12)
```

**Periodicity harmonization:**

| Raw Metric | Tool | Normalization to 0–1 |
|------------|------|----------------------|
| `f0_percent` | RiboCode, RiboseQC | Direct (already 0–1) |
| Frame 0% in CDS | riboWaltz | Direct (already 0–1) |
| `phase_score` | Ribotricer | Direct (already 0–1) |
| `bayes_factor_mean` | PRICE | `min(1, log10(BF) / 3)` |
| Metagene PDF | Ribo-TISH | Manual assessment → 0/0.5/1 |

**Aggregate read-level periodicity score:** Mean of all available normalized periodicity scores.

### Dimension 2: ORF-Level Translation Evidence

| Harmonized Metric | RiboCode | ORFquant | Ribotricer | Ribo-TISH | PRICE |
|-------------------|----------|----------|------------|-----------|-------|
| **Significance** | `adjusted_pval` | `p_value` (DPSS) | `phase_score` + p-value | `FisherQvalue` | `bayes_factor` |
| **Frame preference** | `pval_frame0_vs_frame1/2` | DPSS spectral power | `phase_score` | `FrameQvalue` | Implicit in BF |
| **ORF abundance** | `Psites_frame0_RPKM` | P-site coverage | `read_count`, `read_density` | `RiboPvalue` (TIS only) | ORF profile sum |
| **Coverage** | `Psites_coverage_frame0` | — | `valid_codons_ratio` | — | — |
| **ORF length** | `ORF_length` (nt) | GTF span | `length` (nt) | `AALen` × 3 | BED span |

### Dimension 3: ORF Classification

| Harmonized Category | RiboCode | Ribotricer | Ribo-TISH | ORFquant |
|---------------------|----------|------------|-----------|----------|
| **CDS (annotated)** | `annotated` | `annotated` | `CDS` | `CDS` |
| **uORF** | `uORF` | `uORF` | `uORF` | `uORF` |
| **dORF** | `dORF` | `dORF` | `dORF` | `dORF` |
| **Overlap uORF** | `Overlap_uORF` | — | — | `overlap_uORF` |
| **Overlap dORF** | `Overlap_dORF` | — | — | `overlap_dORF` |
| **Internal** | `internal` | — | — | — |
| **Novel** | `novel` | `novel` | `novel` | `ncRNA` |

### Dimension 4: Cross-Tool Agreement

| Metric | Description | Computation |
|--------|-------------|-------------|
| **ORF overlap rate** | Fraction of ORFs detected by ≥N tools | bedtools intersect between each tool pair |
| **Boundary concordance** | Agreement on start/stop (±nt tolerance) | Compare coordinates of overlapping ORFs |
| **Tool yield** | ORF count per tool per sample | Direct count from output |
| **Pairwise Jaccard** | Overlap coefficient | `|A ∩ B| / |A ∪ B|` |
| **Consensus ORF set** | ORFs supported by ≥2 tools | Intersection-based filtering |
| **P-site offset consensus** | Read lengths with consistent offsets | Median ± tolerance across tools |

---

## 4. QC Module Architecture

### 4.1 Pipeline Integration Point

The QC module runs **after** ORF unification (post `UNIFY_ORF_PREDICTIONS`) and **in parallel with** ORF classification. It consumes outputs from all upstream tools.

```
RiboseQC ──┐
riboWaltz ─┤
            ├──→ ORF_QC_MODULE ──→ qc_report.html
RiboCode ──┤                      ├─→ orf_confidence.tsv
ORFquant ──┤                      ├─→ qc_metrics.tsv
Ribotricer─┤                      ├─→ tool_agreement.tsv
Ribo-TISH ─┤                      └─→ psite_harmonized.tsv
PRICE ─────┤
            │
UNIFY_ORF ─┘
```

**Key principle:** QC module does NOT block the pipeline. It runs with `errorStrategy 'ignore'`.

### 4.2 Module: `ORF_QC`

**Input channels (all optional — QC adapts to available tools):**

```
// === Pure QC tools ===
ch_riboseqc_psites       // RiboseQC *_P_sites_calcs
ch_ribowaltz_offsets     // riboWaltz *_psite_offset.tsv
ch_ribowaltz_frames      // riboWaltz *_frame_distribution.tsv

// === ORF Prediction tools ===
ch_ribocode_collapsed    // RiboCode *_collapsed.txt
ch_orfquant_gtf          // ORFquant *_Detected_ORFs.gtf.gz
ch_ribotricer_tsv        // Ribotricer *_translating_ORFs.tsv
ch_ribotish_pred         // Ribo-TISH *_pred.txt
ch_ribotish_qual         // Ribo-TISH *_qual.txt
ch_rpbp_bayes            // PRICE *bayes-factors.bed.gz

// === Post-unification ===
ch_unified_orfs_bed      // Unified ORFs BED12 from UNIFY_ORF_PREDICTIONS
ch_unified_orfs_meta     // Unified ORFs metadata TSV (with 'tools' column)
```

**Outputs:**
```
qc_report.html           // Interactive HTML QC report
qc_metrics.tsv           // Per-sample aggregate QC metrics
orf_confidence.tsv       // Per-ORF confidence scores (OCS + tier)
tool_agreement.tsv       // Cross-tool pairwise agreement matrix
psite_harmonized.tsv     // Harmonized P-site offsets across all tools
sample_flags.json        // Machine-readable quality flags
```

### 4.3 Processing Stages

```
Stage 1: COLLECT
  ├── Detect which tools produced output (handle missing/failed tools)
  ├── Parse each tool's output format into standardized internal schema
  └── Emit: tool_data{} dict with per-tool structured DataFrames

Stage 2: HARMONIZE (Read-Level)
  ├── P-site offsets: extract per-length offsets from all tools → psite_df
  │   ├── Authority order: riboWaltz > RiboseQC > RiboCode > Ribotricer > Ribo-TISH > PRICE
  │   ├── Consensus = median where ≥2 tools agree within 1 nt
  │   └── Flag discordant read lengths (max delta > 1 nt)
  ├── Periodicity: normalize all scores to 0–1 → periodicity_df
  │   ├── riboWaltz: CDS frame-0% directly
  │   ├── RiboseQC: mean f0 proportion across transcripts
  │   ├── RiboCode: f0_percent directly
  │   ├── Ribotricer: phase_score directly
  │   └── PRICE: min(1, log10(BF) / 3)
  └── Read length distribution: merge all → rldist_df

Stage 3: COMPARE (ORF-Level Cross-Tool)
  ├── For each tool's ORF set, compute pairwise overlap:
  │   ├── bedtools intersect -f 0.5 -r → overlapping ORF pairs
  │   ├── Jaccard index per tool pair
  │   └── Boundary concordance (start/stop delta distribution)
  ├── Build consensus ORF set: ORFs with ≥2 tool support
  ├── Compare ORF classifications across tools (agreement matrix)
  └── Emit: agreement_df, consensus_orfs.bed

Stage 4: SCORE (Per-Unified-ORF Confidence)
  ├── For each unified ORF:
  │   ├── S_translation: best available significance from detecting tool(s)
  │   ├── S_agreement: fraction of tools that detect this ORF
  │   ├── S_coverage: coverage completeness (imputed where unavailable)
  │   ├── S_periodicity: ORF-level periodicity score
  │   └── S_readlevel: sample-wide read-level QC pass rate (global modifier)
  ├── OCS = weighted combination → 0–1 score
  ├── Assign tier: High (≥0.7) / Medium (0.4–0.7) / Low (0.2–0.4) / Uncertain (<0.2)
  └── Emit: orf_confidence.tsv

Stage 5: SAMPLE FLAGS
  ├── LOW_PERIODICITY: aggregate periodicity < 0.5
  ├── NO_PERIODICITY: no read length passes threshold
  ├── P_SITE_DISCORDANCE: >30% lengths disagree
  ├── LOW_YIELD: <100 total ORFs
  ├── LOW_CONSENSUS: <10% ORFs detected by ≥2 tools
  ├── TOOL_FAILURE: tools with no output
  └── Emit: sample_flags.json

Stage 6: REPORT
  ├── HTML report with Plotly/ggplot2:
  │   ├── Tab 1: Sample Dashboard (summary metrics + flags)
  │   ├── Tab 2: Read-Level QC (P-site table, periodicity heatmap, length dist)
  │   ├── Tab 3: ORF-Level QC (confidence distribution, length dist, classification pie)
  │   ├── Tab 4: Cross-Tool Comparison (UpSet plot, Jaccard heatmap, boundary concordance)
  │   └── Tab 5: Per-ORF Table (searchable, filterable by confidence tier)
  └── Export all tabular data as TSV
```

---

## 5. ORF Confidence Score (OCS) Specification

### 5.1 Component Definitions

```
OCS = w1·S_translation + w2·S_agreement + w3·S_coverage + w4·S_periodicity + w5·S_readlevel
```

#### S_translation (weight 0.30) — Normalized Translation Significance

Per-tool normalization to [0, 1]:

| Tool | Raw Metric | Normalization | Clamp Range |
|------|-----------|---------------|-------------|
| RiboCode | `adjusted_pval` | `min(1, -log10(p) / 10)` | [0, 1] |
| ORFquant | `p_value` | `min(1, -log10(p) / 10)` | [0, 1] |
| Ribotricer | `phase_score` | Direct use | [0, 1] |
| Ribo-TISH | `FisherQvalue` | `min(1, -log10(q) / 10)` | [0, 1] |
| PRICE | `bayes_factor` | `min(1, log10(BF) / 5)` | [0, 1] |

If an ORF is detected by multiple tools, take the **maximum** normalized score.

#### S_agreement (weight 0.30) — Cross-Tool Support

Based on the number of ORF prediction tools that detect this ORF (out of the N tools that ran):

| Tools detecting | Score |
|-----------------|-------|
| 1 of 5 | 0.20 |
| 2 of 5 | 0.50 |
| 3 of 5 | 0.75 |
| 4 of 5 | 0.90 |
| 5 of 5 | 1.00 |

Scaled by tool availability: `S_agreement = score × (n_detected / n_available)`.

#### S_coverage (weight 0.20) — Coverage Completeness

| Tool | Metric | Normalization |
|------|--------|---------------|
| RiboCode | `Psites_coverage_frame0` | Direct (0–1, already fraction) |
| Ribotricer | `valid_codons_ratio` | Direct (0–1) |
| Other | Imputed from ORF abundance percentile | Percentile rank / 100 |

#### S_periodicity (weight 0.15) — ORF-Level 3-nt Periodicity

| Tool | Metric | Normalization |
|------|--------|---------------|
| RiboCode | `pval_combined` | `min(1, -log10(p) / 10)` |
| Ribotricer | `phase_score` | Direct |
| Ribo-TISH | `FrameQvalue` | `min(1, -log10(q) / 10)` |
| PRICE | `bayes_factor` | `min(1, log10(BF) / 5)` |
| ORFquant | DPSS p-value | `min(1, -log10(p) / 10)` |

#### S_readlevel (weight 0.05) — Sample-Wide Read QC Modifier

A global modifier based on read-level QC:
- Computed as the mean normalized periodicity score across all read lengths (from Dimension 1)
- If >0.7: 1.0 (good data)
- If 0.5–0.7: 0.8 (marginal data)
- If 0.3–0.5: 0.5 (poor data)
- If <0.3: 0.2 (bad data — all ORFs penalized)

### 5.2 Confidence Tiers

| Tier | OCS Range | Interpretation | Recommended Use |
|------|-----------|----------------|-----------------|
| **High** | ≥ 0.7 | Multi-tool support, strong periodicity, good coverage | Downstream analysis, validation studies |
| **Medium** | 0.4–0.7 | Moderate evidence, 1–2 tools, acceptable metrics | Exploratory analysis |
| **Low** | 0.2–0.4 | Weak evidence, single tool, marginal significance | Requires orthogonal validation |
| **Uncertain** | < 0.2 | Below thresholds, likely artifact | Filter out for most analyses |

---

## 6. Output File Specifications

### 6.1 `psite_harmonized.tsv`

```
read_length | ribowaltz_offset | riboseqc_offset | ribocode_offset | ribotricer_offset | ribotish_offset | price_offset | consensus_offset | n_tools | max_delta | flag
28          | 12               | 12              | 12              | 12                | 12              | —            | 12                | 5       | 0         | OK
29          | 12               | 12              | 12              | 13                | 12              | —            | 12                | 5       | 1         | WARN
30          | 13               | 13              | 12              | 13                | 13              | 12           | 13                | 6       | 1         | OK
```

### 6.2 `orf_confidence.tsv`

```
orf_id            | unified_orf_id | ocs  | tier      | s_translation | s_agreement | s_coverage | s_periodicity | s_readlevel | detecting_tools     | orf_type | orf_length_nt
ENST001_100_400_100| ORF00001       | 0.82 | High      | 0.85          | 1.0         | 0.75       | 0.70          | 1.0         | RiboCode,Ribotricer| uORF     | 300
ENST002_50_200_50  | ORF00002       | 0.35 | Low       | 0.20          | 0.2         | —          | —             | 0.5         | RiboCode            | novel    | 150
```

### 6.3 `tool_agreement.tsv` (Pairwise)

```
tool_a     | tool_b     | jaccard | overlap_count | a_only | b_only | boundary_concordance_5p | boundary_concordance_3p
RiboCode   | Ribotricer | 0.45    | 320           | 180    | 210    | 0.72                    | 0.68
RiboCode   | Ribo-TISH  | 0.38    | 260           | 240    | 190    | 0.65                    | 0.70
...
```

### 6.4 `sample_flags.json`

```json
{
  "sample_id": "HEK293_rep1",
  "flags": [
    {"flag": "P_SITE_DISCORDANCE", "severity": "WARNING", "detail": "2/5 read lengths disagree on P-site offset"},
    {"flag": "LOW_CONSENSUS", "severity": "INFO", "detail": "8% ORFs detected by ≥2 tools"}
  ],
  "summary": {
    "total_orfs_unified": 1250,
    "consensus_orfs": 100,
    "tool_count_ran": 5,
    "tool_count_failed": 0,
    "mean_ocs": 0.42,
    "periodicity_pass_rate": 0.67,
    "psite_consensus_rate": 0.60
  }
}
```

---

## 7. Read-Level QC Report Contents

### 7.1 P-site Offset Harmonization Table

riboWaltz serves as the **primary reference** (best algorithm). Each row = one read length.

| Column | Source | Description |
|--------|--------|-------------|
| `read_length` | — | RPF length (nt) |
| `ribowaltz_offset` | `corrected_offset_from_5` | Primary reference |
| `riboseqc_offset` | `P_sites_calcs` | Fallback reference |
| `[tool]_offset` | Respective tool output | Cross-validation |
| `consensus_offset` | Median of available | Final harmonized value |
| `n_tools` | Count | Number of tools providing data |
| `max_delta` | max − min across tools | Discordance metric |
| `flag` | OK / WARN / FAIL | Quality flag |

### 7.2 Periodicity Assessment Table

| Read Length | riboWaltz CDS f0% | RiboseQC f0% | RiboCode f0% | Ribotricer Phase | PRICE log10(BF) | Aggregate Score | Periodic? |
|-------------|-------------------|--------------|--------------|------------------|-----------------|-----------------|-----------|
| 28 | 78% | 75% | 82% | 0.72 | 3.5 | 0.78 | ✅ YES |
| 29 | 71% | 68% | 75% | 0.65 | 2.8 | 0.70 | ✅ YES |
| 30 | 62% | 58% | 68% | 0.55 | 1.2 | 0.58 | ⚠️ MARGINAL |
| 31 | 45% | 42% | 52% | 0.38 | 0.3 | 0.40 | ❌ NO |

- **Periodic**: Aggregate score ≥ 0.65 AND ≥2 tools agree
- **Marginal**: Aggregate score 0.5–0.65 OR only 1 tool available
- **Non-periodic**: Aggregate score < 0.5

---

## 8. Sample-Level QC Summary

### 8.1 Key Aggregate Metrics

| Metric | Description | Derivation |
|--------|-------------|------------|
| `n_tools_ran` | Number of ORF prediction tools that produced output | Count |
| `n_tools_failed` | Tools that failed or produced no ORFs | Count |
| `total_orfs_raw` | Sum of ORF counts across all tools (pre-unification) | Sum |
| `total_orfs_unified` | ORFs after deduplication | From UNIFY_ORF_PREDICTIONS |
| `consensus_orfs` | ORFs detected by ≥2 prediction tools | bedtools intersect |
| `consensus_rate` | consensus_orfs / total_orfs_unified | Ratio |
| `mean_ocs` | Mean OCS across all unified ORFs | Mean |
| `high_confidence_rate` | Fraction of ORFs with OCS ≥ 0.7 | Ratio |
| `periodicity_pass_rate` | Fraction of read lengths classified as "Periodic" | Ratio |
| `psite_consensus_rate` | Fraction of read lengths with max_delta ≤ 1 nt | Ratio |
| `cds_recovery_rate` | Fraction of annotated CDS detected by ≥1 tool | bedtools intersect with reference |
| `novel_orf_rate` | Fraction of ORFs classified as "novel" | From classification |

### 8.2 Quality Flags

| Flag | Condition | Severity | Suggested Action |
|------|-----------|----------|-----------------|
| `NO_PERIODICITY` | 0 read lengths periodic | 🔴 CRITICAL | Do not trust ORF calls; check library quality |
| `LOW_PERIODICITY` | <50% read lengths periodic | 🟡 WARNING | Use high-confidence ORFs only |
| `P_SITE_DISCORDANCE` | >30% lengths disagree on offset | 🟡 WARNING | Manual inspection of metagene profiles |
| `LOW_CONSENSUS` | <10% ORFs detected by ≥2 tools | 🟡 WARNING | Low agreement suggests noisy data |
| `LOW_YIELD` | <100 unified ORFs | 🔵 INFO | May be expected for low-complexity samples |
| `TOOL_FAILURE` | ≥1 prediction tool failed | 🔵 INFO | Review tool logs for cause |
| `HIGH_NOVEL_RATE` | >80% ORFs classified as novel | 🔵 INFO | May indicate poorly annotated genome |
| `P_SITE_FALLBACK_USED` | riboWaltz used as fallback | 🔵 INFO | RiboseQC didn't produce valid P-site data |

---

## 9. Implementation Plan

### Phase 1: Data Schema & Extraction (`bin/extract_orf_qc_metrics.py`)

**Input:** All tool output files (auto-detected from file patterns in a results directory).

**Processing:**
```
1. Auto-detect which tools produced output:
   - Check for file existence by glob pattern
   - Record tool status: OK / EMPTY / FAILED / SKIPPED

2. Parse each tool's output into pandas DataFrames:
   ribocode    → parse _collapsed.txt (TSV)
   riboseqc    → parse _P_sites_calcs (TSV with # comments)
   ribowaltz   → parse _psite_offset.tsv, _frame_distribution.tsv
   ribotricer  → parse _translating_ORFs.tsv
   ribotish    → parse _pred.txt, parse *.para.py (eval)
   orfquant    → parse _Detected_ORFs.gtf.gz (extract attributes)
   price       → parse bayes-factors.bed.gz

3. Validate: check required columns exist, handle empty files
```

**Output:** `tool_data.json` — structured per-tool DataFrames in JSON format.

### Phase 2: Harmonization (`bin/harmonize_orf_qc.py`)

**Processing:**
```
1. P-site harmonization:
   - Extract per-length offsets from all tools
   - Compute consensus (median), max_delta, flag per length
   - Output: psite_harmonized.tsv

2. Periodicity harmonization:
   - Normalize all scores to 0–1
   - Classify each read length as Periodic/Marginal/Non-periodic
   - Compute sample-wide aggregate score

3. Read-length distribution merge
```

### Phase 3: Cross-Tool Comparison (`bin/compare_orf_tools.py`)

**Processing:**
```
1. For each tool with ORF output, convert to BED6 format
2. Pairwise bedtools intersect:
   - -f 0.5 -r for reciprocal overlap ≥50%
   - Compute Jaccard index
   - Compute boundary concordance (delta_start, delta_stop distributions)
3. Build consensus ORF set (≥2 tools)
4. Classification agreement matrix
```

### Phase 4: ORF Confidence Scoring (`bin/score_orf_confidence.py`)

**Processing:**
```
For each unified ORF:
  1. Identify which tools detect it (from unified metadata 'tools' column)
  2. S_translation = max normalized significance across detecting tools
  3. S_agreement = score based on n_detecting / n_available
  4. S_coverage = from tool data (imputed if missing)
  5. S_periodicity = from tool data
  6. S_readlevel = global modifier from Phase 2
  7. OCS = weighted sum
  8. Assign tier
```

### Phase 5: Report Generation (`bin/generate_orf_qc_report.py`)

**Framework:** Python with Plotly (interactive HTML) or R with RMarkdown/flexdashboard.

**Report tabs:**
1. **Dashboard**: Key metrics, flags, tool status summary
2. **Read QC**: P-site harmonization table + periodicity heatmap + length distribution
3. **ORF QC**: OCS distribution histogram, confidence tier pie, length distribution by tier
4. **Cross-Tool**: UpSet plot, Jaccard heatmap, boundary concordance scatter
5. **Detail Table**: Searchable DataTable with all ORFs and their scores

### Phase 6: Nextflow Module (`modules/local/orf_qc/main.nf`)

```nextflow
process ORF_QC {
    tag "$meta.id"
    label 'process_medium'
    errorStrategy 'ignore'  // QC failure must not block pipeline

    input:
    // Accept all channels; QC adapts to what's available
    tuple val(meta), path(riboseqc_psites)      // optional
    tuple val(meta), path(ribowaltz_offsets)     // optional
    tuple val(meta), path(ribowaltz_frames)      // optional
    tuple val(meta), path(ribocode_collapsed)    // optional
    tuple val(meta), path(orfquant_gtf)          // optional
    tuple val(meta), path(ribotricer_tsv)        // optional
    tuple val(meta), path(ribotish_pred)         // optional
    tuple val(meta), path(ribotish_qual)         // optional
    tuple val(meta), path(rpbp_bayes)            // optional
    tuple val(meta), path(unified_bed)           // required (post-unification)
    tuple val(meta), path(unified_meta)          // required

    output:
    tuple val(meta), path("qc_report.html")
    tuple val(meta), path("qc_metrics.tsv")
    tuple val(meta), path("orf_confidence.tsv")
    tuple val(meta), path("tool_agreement.tsv")
    tuple val(meta), path("psite_harmonized.tsv")
    tuple val(meta), path("sample_flags.json")
    path "versions.yml"

    script:
    """
    extract_orf_qc_metrics.py --input-dir . --output tool_data.json
    harmonize_orf_qc.py --input tool_data.json --output-prefix ${meta.id}
    compare_orf_tools.py --input tool_data.json --unified ${unified_bed} --output-prefix ${meta.id}
    score_orf_confidence.py --unified ${unified_meta} --tool-data tool_data.json --read-qc ${meta.id}_periodicity.tsv --output ${meta.id}_orf_confidence.tsv
    generate_orf_qc_report.py --output ${meta.id}_qc_report.html [all TSV inputs]
    """
}
```

### Phase 7: MultiQC Integration (Optional)

Generate a MultiQC-compatible YAML section so ORF QC metrics appear in the existing MultiQC report:
```yaml
id: "orf_qc"
section_name: "ORF Prediction QC"
plot_type: "table"
data:
  sample_name: ["HEK293_rep1"]
  total_orfs_unified: [1250]
  consensus_orfs: [100]
  mean_ocs: [0.42]
  ...
```

---

## 10. Configuration Parameters

```nextflow
// ORF QC Module
params.orf_qc_enabled              = true
params.orf_qc_min_consensus_tools  = 2       // Min tools for consensus ORF
params.orf_qc_confidence_weights   = [0.30, 0.30, 0.20, 0.15, 0.05]  // OCS weights
params.orf_qc_periodicity_min_f0   = 0.6     // Min f0 proportion for "periodic"
params.orf_qc_periodicity_min_agree= 2       // Min tools agreeing for "periodic"
params.orf_qc_offset_max_delta     = 1       // Max nt difference in P-site consensus
params.orf_qc_report_format        = 'html'  // html, json, or both
params.orf_qc_skip                 = false   // Skip entire QC module
```

---

## 11. Key Design Decisions

| # | Decision | Rationale |
|---|----------|-----------|
| 1 | **QC module is non-blocking** | Pipeline must complete even if QC fails; `errorStrategy 'ignore'` |
| 2 | **Post-unification placement** | Operates on unified ORFs to avoid per-tool BED format complexity |
| 3 | **riboWaltz as P-site authority** | Two-step coherence correction is the most accurate algorithm |
| 4 | **Harmonization, not replacement** | Tool-specific metrics preserved alongside harmonized scores |
| 5 | **OCS as primary per-ORF output** | Enables confidence-based filtering in downstream analysis |
| 6 | **Read-level + ORF-level QC** | Both dimensions needed: good reads can still produce bad ORF calls |
| 7 | **Global read-level modifier (5%)** | Sample-wide data quality should slightly influence all ORF scores |
| 8 | **Extensible tool registry** | Adding a new tool = adding a parser + normalization function |
| 9 | **Machine-readable flags** | `sample_flags.json` enables automated pipeline quality gating |
| 10 | **Auto-detect available tools** | QC adapts to whatever completed successfully; no hard dependency |
