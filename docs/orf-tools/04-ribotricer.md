# Ribotricer — Reference-Guided ORF Detection

## Overview

**Ribotricer** is a Python-based tool for detecting translating ORFs from Ribo-seq data using a reference-guided approach. It computes a **phase score** based on Fourier coherence between the observed P-site signal and an idealized 3-nt periodic pattern (1-0-0). It also provides metagene analysis for P-site offset determination.

- **Version in pipeline**: Latest (via Bioconda)
- **Source**: https://github.com/smithlabcode/ribotricer
- **Publication**: Choudhary S., Li W., Smith A.D., *BMC Genomics*, 2021
- **Language**: Python 3
- **Key Dependencies**: numpy, scipy, pandas, pysam, quicksect, matplotlib

---

## Analysis Steps

### Step 1: `prepare_orfs` — ORF Index Construction

**Inputs:**
- GTF annotation file
- Genome FASTA file

**Processing:**
- Extracts all annotated CDS regions
- Identifies all possible ORFs from transcript sequences
- Builds interval tree index for rapid lookup

**Output:** Ribotricer index file (TSV)

### Step 2: `detect_orfs` — Translation Detection

**Inputs:**
- Genome-aligned BAM file
- Ribotricer index file

**Sub-steps:**

#### 2a. Protocol Inference (optional)
- Infers strandedness from read distribution on annotated CDS

#### 2b. BAM Splitting
- Splits BAM by read length and strand

#### 2c. Metagene Analysis
- Builds metagene coverage profiles around start/stop codons for each read length
- Computes **phase score at 5′ and 3′ ends** for each read length
- Generates metagene profile plots

#### 2d. P-site Offset Determination
- **Cross-correlation method**: Aligns metagene profiles across read lengths
- Selects the most abundant read length as reference
- Computes cross-correlation lag for all other lengths
- Applies typical offset correction

#### 2e. ORF Coverage Computation
- Shifts reads by P-site offsets and merges across read lengths
- Computes per-nucleotide coverage for each candidate ORF

#### 2f. Phase Score Calculation (Translation Detection)
- For each ORF, computes **phase score** (`phasescore()`)
- Compares to ideal 1-0-0 periodic signal using **Fourier coherence**
- Classifies ORF as "translating" or "nontranslating" based on multiple thresholds

---

## Phase Score Algorithm (`statistics.py`)

The core QC metric is computed as follows:

```
For each reading frame (0, 1, 2):
  1. Extract coverage values starting at that frame
  2. Normalize each codon triplet by its vector magnitude
  3. Compute spectral coherence between normalized signal and [1,0,0] pattern
  4. At frequency 1/3 (period = 3 nt), extract coherence value
  5. Return sqrt(max(coherence across frames)), valid_codons

Phase Score = sqrt(Cxy) at f = 1/3
```

**Statistical significance**: p-value computed using non-central chi-squared distribution:
```
x = 2 * N² * phase_score / (N - 1)
p = ncx2.sf(x, df=2, nc=2/(N-1))
```

---

## Output Files

| File | Content |
|------|---------|
| `*_translating_ORFs.tsv` | Detected ORFs with QC metrics |
| `*_metagene_profiles_5p.tsv` | Metagene profiles aligned at start codon |
| `*_metagene_profiles_3p.tsv` | Metagene profiles aligned at stop codon |
| `*_psite_offsets.txt` | P-site offsets per read length |
| `*_read_length_distribution.pdf` | Read length distribution plot |
| `*_metagene_profiles.pdf` | Metagene profile plots |
| `*_pos.wig` / `*_neg.wig` | P-site shifted WIG tracks |

### `*_translating_ORFs.tsv` Columns

| Column | Description |
|--------|-------------|
| `ORF_ID` | ORF identifier |
| `ORF_type` | Classification (annotated, uORF, dORF, novel, etc.) |
| `status` | `translating` or `nontranslating` |
| **`phase_score`** | **3-nt periodicity score (Fourier coherence)** |
| `read_count` | Total reads in ORF |
| `length` | ORF length (nt) |
| `valid_codons` | Number of codons with sufficient signal |
| `valid_codons_ratio` | Fraction of valid codons |
| `read_density` | Average reads per codon |
| `transcript_id` | Associated transcript |
| `gene_id` / `gene_name` | Associated gene |
| `chrom` / `strand` | Genomic location |
| `start_codon` | Start codon sequence |
| `profile` | Per-nucleotide coverage vector |

### Translation Detection Thresholds (default)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `phase_score_cutoff` | 0.4 | Minimum phase score for "translating" |
| `min_valid_codons` | 10 | Minimum valid codon positions |
| `min_reads_per_codon` | 3.0 | Minimum average reads per codon |
| `min_valid_codons_ratio` | 0.5 | Minimum fraction of valid codons |
| `min_density_over_orf` | 5.0 | Minimum average reads per codon over entire ORF |

---

## Pipeline Integration

### Module in riboseq pipeline

Ribotricer is run as part of the ORF prediction workflow. The pipeline uses the `ribotricer` Bioconda package directly.

**Inputs:**
- Genome-aligned BAM (sORF-filtered)
- Reference GTF + FASTA

**Outputs used downstream:**
- `*_translating_ORFs.tsv` → ORF unification module
- Metagene profiles → QC assessment

**Parameters:**
| Parameter | Default | Description |
|-----------|---------|-------------|
| `--skip_ribotricer` | `false` | Skip Ribotricer |

---

## QC Module Design Implications

### Strengths for QC integration:
1. **Excellent periodicity quantification**: Phase score based on signal processing (Fourier coherence) is statistically rigorous
2. **Per-ORF p-values**: Provides statistical significance for each ORF's periodicity
3. **Multi-threshold classification**: Uses phase score + codon coverage + density thresholds
4. **Data-driven P-site offsets**: Cross-correlation method for offset determination
5. **Metagene profiles**: Rich visualization of read distribution around start/stop codons

### Metrics to carry into unified QC module:
| Metric | Source | Type |
|--------|--------|------|
| **`phase_score`** | `*_translating_ORFs.tsv` | ORF-level periodicity |
| `status` (translating/nontranslating) | `*_translating_ORFs.tsv` | ORF classification |
| `read_count` / `read_density` | `*_translating_ORFs.tsv` | ORF-level abundance |
| `valid_codons` / `valid_codons_ratio` | `*_translating_ORFs.tsv` | Coverage completeness |
| `ORF_type` | `*_translating_ORFs.tsv` | ORF classification |
| P-site offsets per read length | `*_psite_offsets.txt` | Read-level offset |
| Metagene phase scores (5′/3′) | `*_metagene_profiles_*.tsv` | Read-level periodicity |

### Limitations:
1. **Reference-dependent**: Only detects ORFs in annotated transcript regions
2. **No de-novo ORF discovery**: Unlike RiboCode, relies on pre-defined ORF index
3. **Phase score affected by coverage depth**: Low-coverage ORFs may have unstable scores
4. **Single-threshold classification**: ORF is either "translating" or not; no confidence gradation
