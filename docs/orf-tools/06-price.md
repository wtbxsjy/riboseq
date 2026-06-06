# PRICE (rp-bp) — Bayesian ORF Prediction

## Overview

**Rp-Bp** (Ribosome profiling with Bayesian predictions), also known as **PRICE**, is an unsupervised Bayesian approach to predict translated ORFs from ribosome profiles. It uses Gaussian process models to distinguish periodic (translated) from non-periodic (untranslated) ribosome footprint patterns. The pipeline comes with interactive dashboards for QC and ORF discovery.

- **Version in pipeline**: 4.0.1
- **Source**: https://github.com/dieterich-lab/rp-bp
- **Publication**: Malone B. et al., *Nucleic Acids Research*, 2017 (doi: 10.1093/nar/gkw1350)
- **Language**: Python 3
- **Key Dependencies**: numpy, scipy, pandas, pystan/pycmdstan, pbiotools

---

## Analysis Steps

The PRICE pipeline has 6 sequential steps, all run within the `RPBP_PREDICT` module:

### Step 1: `extract-metagene-profiles` — Metagene Profile Extraction

**Inputs:**
- Genome-aligned BAM file
- ORF genomic coordinates file

**Processing:**
- For each read length, extracts aggregate coverage profiles around annotated CDS start/stop codons
- Profiles include leader (upstream) and trailer (downstream) regions

**Output:** `*.metagene-profiles.csv.gz`

### Step 2: `estimate-metagene-profile-bayes-factors` — Periodicity Assessment

**Inputs:**
- Metagene profiles CSV
- **Periodic Gaussian process models** (pre-trained Stan models)
- **Non-periodic Gaussian process models**

**Processing:**
- For each read length, fits both periodic and non-periodic GP models to the metagene profile
- Computes **Bayes factor** = P(data | periodic model) / P(data | non-periodic model)
- A high Bayes factor means strong evidence for 3-nt periodicity

**Output:** `*.metagene-profile-bayes-factors.csv.gz`

**QC metrics per read length:**
| Column | Description |
|--------|-------------|
| `length` | Read length |
| `offset` | Candidate P-site offset |
| `profile_sum` | Total signal at this offset |
| `profile_peak` | Peak signal height |
| `bayes_factor_mean` | **Mean Bayes factor (periodic vs non-periodic)** |
| `bayes_factor_var` | Variance of Bayes factor estimate |

### Step 3: `select-periodic-offsets` — P-site Offset Selection

**Input:** Metagene Bayes factors CSV

**Processing:**
- Groups Bayes factors by read length
- For each length, selects the offset with the **highest periodic Bayes factor peak**
- This gives the most probable P-site position for each read length

**Output:** `*.periodic-offsets.csv.gz`

**Selected offset columns:**
| Column | Description |
|--------|-------------|
| `length` | Read length |
| `highest_peak_offset` | **Selected P-site offset** |
| `highest_peak_peak` | Peak profile height |
| `highest_peak_profile_sum` | Total profile signal |
| `highest_peak_bf_mean` | **Bayes factor at selected offset** |
| `highest_peak_bf_var` | Bayes factor variance |

**Periodicity QC threshold options:**
| Parameter | Purpose |
|-----------|---------|
| `min_metagene_profile_count` | Minimum reads at a length to consider |
| `min_metagene_bf_mean` | Minimum periodic Bayes factor |
| `max_metagene_bf_var` | Maximum Bayes factor variance |
| `min_metagene_bf_likelihood` | Minimum likelihood threshold |

### Step 4: `extract-orf-profiles` — ORF Signal Extraction

**Inputs:**
- BAM file
- ORF coordinates (genomic + exonic)
- Selected read lengths and P-site offsets

**Processing:**
- Shifts reads by P-site offsets
- Extracts per-nucleotide P-site density for each candidate ORF
- Generates sparse matrix of ORF × position profiles

**Output:** `*.profiles.mtx.gz` (sparse matrix)

### Step 5: `estimate-orf-bayes-factors` — Translation Assessment

**Inputs:**
- ORF profile matrix
- **Translated GP models** (pre-trained)
- **Untranslated GP models** (pre-trained)

**Processing:**
- Fits translated and untranslated GP models to each ORF's P-site profile
- Computes **Bayes factor** = P(profile | translated model) / P(profile | untranslated model)
- High Bayes factor → strong evidence of translation
- Optional: chi-square only mode (skips GP model fitting)

**Output:** `*.bayes-factors.bed.gz` (BED format with Bayes factor scores)

**ORF-level QC thresholds:**
| Parameter | Default | Description |
|-----------|---------|-------------|
| `min_bf_mean` | — | Minimum Bayes factor for calling translated |
| `max_bf_var` | — | Maximum acceptable Bayes factor variance |
| `min_bf_likelihood` | — | Minimum likelihood threshold |
| `min_orf_length` | — | Minimum ORF length (nt) |
| `min_profile` | — | Minimum P-site count in ORF |
| `chisq_alpha` | — | Chi-square significance level |

### Step 6: `select-final-prediction-set` — Final ORF Selection

**Inputs:**
- Bayes factors BED
- Candidate ORFs
- Genome FASTA

**Processing:**
- Filters ORFs by Bayes factor thresholds
- For each stop codon: selects longest ORF (`--select-longest-by-stop`)
- Among overlapping ORFs: selects highest Bayes factor (`--select-best-overlapping`)
- Generates DNA and protein sequences for predicted ORFs

**Outputs:**
- `*.predicted-orfs.bed.gz` — Final predicted ORFs
- `*.predicted-orfs.dna.fa` — ORF DNA sequences
- `*.predicted-orfs.protein.fa` — ORF protein sequences

---

## Output Files

| File | Content |
|------|---------|
| `*.metagene-profiles.csv.gz` | Raw metagene profiles per read length |
| `*.metagene-profile-bayes-factors.csv.gz` | Periodic vs non-periodic Bayes factors |
| `*.periodic-offsets.csv.gz` | Selected P-site offsets with Bayes factors |
| `*.profiles.mtx.gz` | ORF × position P-site density matrix |
| `*.bayes-factors.bed.gz` | Per-ORF translated vs untranslated Bayes factors |
| `*.predicted-orfs.bed.gz` | Final predicted ORFs |
| `*.predicted-orfs.dna.fa` | ORF nucleotide sequences |
| `*.predicted-orfs.protein.fa` | ORF amino acid sequences |

---

## QC Metrics Summary

### A. Read-Level QC

| Metric | Source | Type |
|--------|--------|------|
| Read length distribution | `*.metagene-profiles.csv.gz` | Abundance per length |
| **Periodic Bayes factor per (length, offset)** | `*.metagene-profile-bayes-factors.csv.gz` | 3-nt periodicity evidence |
| Profile peak/sum per offset | `*.metagene-profile-bayes-factors.csv.gz` | Signal strength |
| Selected P-site offset per length | `*.periodic-offsets.csv.gz` | P-site position |

### B. ORF-Level QC

| Metric | Source | Type |
|--------|--------|------|
| **Translated vs untranslated Bayes factor** | `*.bayes-factors.bed.gz` | Translation evidence |
| ORF P-site count | `*.profiles.mtx.gz` | Abundance |
| ORF length | BED file | Size |
| Classification (filtered/unfiltered) | BED file | Redundancy resolution |

### C. Bayesian Framework Characteristics

The key distinction of PRICE's QC is its **Bayesian model comparison approach**:

1. **Periodicity assessment**: GP model comparison rather than heuristic thresholds
2. **Translation assessment**: Full Bayesian model comparison at the ORF level
3. **Uncertainty quantification**: Bayes factor variance captures estimation uncertainty
4. **Hierarchical decisions**: Length selection → offset selection → ORF selection

---

## Pipeline Integration

### Module: `RPBP_PREDICT`

**Input channels:**
- `tuple val(meta), path(bam), path(bai)` — Genome-aligned BAM
- `path(orfs_genomic)` — ORF genomic coordinates
- `path(orfs_exons)` — ORF exon coordinates

**Processing note:** The pipeline runs all 6 steps within a single module instance.

**Output channels:**
| Channel | Type | Content |
|---------|------|---------|
| `predictions` | `tuple val(meta), path("*.predicted-orfs.bed.gz")` | Final ORFs |
| `predicted_dna` | `tuple val(meta), path("*.predicted-orfs.dna.fa")` | DNA sequences |
| `predicted_protein` | `tuple val(meta), path("*.predicted-orfs.protein.fa")` | Protein sequences |
| `bayes_factors` | `tuple val(meta), path("*.bayes-factors.bed.gz")` | ORF Bayes factors |

**Parameters:**
| Parameter | Default | Description |
|-----------|---------|-------------|
| `--run_rpbp` | `false` | Enable PRICE (opt-in) |
| `--rpbp_container` | (auto) | Custom container path |

---

## QC Module Design Implications

### Strengths for QC integration:
1. **Rigorous statistical framework**: Bayesian model comparison is more principled than ad-hoc thresholds
2. **Full uncertainty quantification**: Bayes factor variance provides confidence estimates
3. **Two-level periodicity assessment**: Read-length level (metagene periodicity) + ORF level (translation)
4. **Automatic offset selection**: Data-driven with quality metrics
5. **Pre-trained models**: Consistent comparisons across samples/datasets

### Metrics to carry into unified QC module:
| Metric | Source | Type |
|--------|--------|------|
| Periodic BF per read length | `*.metagene-profile-bayes-factors.csv.gz` | Read-level periodicity |
| Selected P-site offset per length | `*.periodic-offsets.csv.gz` | Read-level offset |
| **Translated vs untranslated BF** | `*.bayes-factors.bed.gz` | ORF-level translation evidence |
| Bayes factor variance | `*.bayes-factors.bed.gz` | ORF-level confidence |
| Number of predicted ORFs | Output count | Yield summary |
| ORF length distribution | BED file | Size distribution |

### Limitations:
1. **Computationally expensive**: GP model fitting with MCMC is slow; uses pre-compiled Stan models
2. **Model dependency**: Requires pre-trained periodic/nonperiodic/translated/untranslated models
3. **Opt-in only**: Not enabled by default (`--run_rpbp`)
4. **No explicit ORF classification**: Doesn't classify ORF types (uORF, dORF, etc.) in standard output
5. **Intermediate files not exposed**: Internal step outputs not surfaced through pipeline channels
6. **Complex configuration**: Requires YAML config file for full parameterization
