# RiboCode — ORF Prediction Tool Analysis

## Overview

**RiboCode** is a Python-based computational algorithm for identifying genome-wide translated ORFs using ribosome-profiling (Ribo-seq) data. It detects ORFs by analyzing the 3-nucleotide (3-nt) periodicity of P-site densities along transcripts.

- **Version in pipeline**: 1.2.11 (Bioconda)
- **Source**: https://github.com/xryanglab/RiboCode
- **Publication**: Xiao Z. et al., *Nucleic Acids Research*, 2018 (doi: 10.1093/nar/gky179)
- **Language**: Python 2/3
- **Key Dependencies**: pysam, Biopython, NumPy, SciPy, statsmodels, matplotlib, HTSeq

---

## Analysis Steps

RiboCode performs three sequential steps:

### Step 1: `prepare_transcripts` — Transcript Annotation Preparation

**Inputs:**
- GTF annotation file (GENCODE/Ensembl format)
- Genome FASTA file

**Processing:**
- Parses GTF to build gene→transcript→exon hierarchy
- Extracts transcript sequences from genome FASTA
- Identifies annotated start codons (CDS start) and stop codons
- Serializes gene/transcript objects to `transcripts.pickle`

**Output:** `RiboCode_annot/` directory containing:
- `transcripts.pickle` — serialized gene/transcript objects
- `transcripts_sequence.fa` — transcript sequences in FASTA format

### Step 2: `metaplots` — P-site offset determination & periodicity analysis

**Inputs:**
- `RiboCode_annot/` from Step 1
- Transcriptome-aligned BAM file(s)

**Processing:**
1. Selects one principal transcript per coding gene (CCDS preferred, then longest)
2. For each read, computes distance from 5′ end to annotated start/stop codons
3. For each read length `L`, builds a metagene density profile (positions −50 to +50 around start/stop codons)
4. **Automatic P-site prediction** (`_predict_psite()` function):
   - Scans offset positions 4 to min(L, 20) from 5′ end
   - For each candidate offset, extracts frame-0, frame-1, frame-2 densities from the downstream 51 nt window
   - Applies two filters:
     - **Frame-0 proportion**: f0/(f0+f1+f2) > 0.65 (default)
     - **Wilcoxon signed-rank test**: f0 > f1 and f0 > f2 with p < 0.001 (default)
   - Selects offset with maximum f0 density passing both filters

**Outputs:**
- `*_pre_config.txt` — P-site configuration with QC metrics per read length
- `*_metaplots.pdf` — Metagene profile plots (5′→start/stop codon distance for each read length)
- `*_readlength_distribution.pdf` — Read length distribution histogram

### Step 3: `RiboCode` — ORF Detection

**Inputs:**
- `RiboCode_annot/` from Step 1
- `config.txt` — P-site configuration (from Step 2 or user-defined)
- Transcriptome BAM file(s)

**Processing:**
1. Counts P-sites at each transcript nucleotide position using the P-site offsets
2. For each transcript:
   - Finds all possible ORFs (start→stop codon combinations) using `orf_finder()`
   - For each stop codon, determines the best start codon via `start_check()`:
     - If `--longest-orf yes`: uses most distal AUG (canonical behavior)
     - If `--longest-orf no`: tests each candidate start codon for significant frame-0 enrichment in the downstream region
   - Extracts frame-0/1/2 P-site arrays for each candidate ORF
   - Performs statistical tests (see QC Metrics below)
   - Filters ORFs with combined p-value ≤ cutoff (default 0.05)
3. Collapses ORFs sharing the same stop codon across transcript isoforms (keeps most upstream start)
4. Applies multiple testing correction (default: Benjamini-Hochberg FDR)

**Outputs:**
- `<name>.txt` — All predicted ORFs (per transcript)
- `<name>_collapsed.txt` — Non-redundant ORFs (per gene, one per stop codon)
- `<name>.gtf` — ORF annotations in GTF format (optional, `-g`)
- `<name>.bed` — ORF annotations in BED format (optional, `-b`)
- `<name>_collapsed.gtf` / `<name>_collapsed.bed` — Collapsed versions
- `*_pie_chart.pdf` — ORF category distribution pie chart

---

## QC Metrics

### A. Read-Level QC (from `metaplots` / `_pre_config.txt`)

| Metric | Column | Description |
|--------|--------|-------------|
| Read length | `read_length` | RPF read length (nt) |
| Read proportion | `proportion` | Fraction of total mapped reads at this length |
| **Predicted P-site offset** | `predicted_psite` | Distance from 5′ end to P-site (nt) |
| Frame-0 P-site sum | `f0_sum` | Total P-site density in frame 0 around start codons |
| Frame-1 P-site sum | `f1_sum` | Total P-site density in frame 1 |
| Frame-2 P-site sum | `f2_sum` | Total P-site density in frame 2 |
| **3-nt periodicity** | `f0_percent` | **f0/(f0+f1+f2)** — key periodicity indicator (threshold: >65%) |
| Frame-0 vs Frame-1 p-value | `pvalue1` | Wilcoxon test: f0 > f1 |
| Frame-0 vs Frame-2 p-value | `pvalue2` | Wilcoxon test: f0 > f2 |
| Combined p-value | `pvalue_combined` | Stouffer's method for combining pvalue1 and pvalue2 |

**Periodicity check algorithm** (`metaplots.py:_predict_psite`):
```
For each read length L:
  For psite_offset in [4, min(L, 20)]:
    f0, f1, f2 = extract_frames(downstream_51nt)
    if f0/(f0+f1+f2) < 0.65: REJECT (insufficient periodicity)
    if wilcoxon(f0 > f1) >= 0.001: REJECT
    if wilcoxon(f0 > f2) >= 0.001: REJECT
    ACCEPT: save offset with max f0
```

### B. ORF-Level QC (from `_collapsed.txt`)

| Metric | Column | Description |
|--------|--------|-------------|
| ORF ID | `ORF_ID` | `{gene_id}_{gstart}_{gstop}_{aa_length}` |
| ORF type | `ORF_type` | See classification below |
| Frame-0 P-site sum | `Psites_sum_frame0` | Total P-sites in frame 0 within ORF |
| Frame-1 P-site sum | `Psites_sum_frame1` | Total P-sites in frame 1 |
| Frame-2 P-site sum | `Psites_sum_frame2` | Total P-sites in frame 2 |
| Frame-0 coverage | `Psites_coverage_frame0` | Fraction of frame-0 codons with ≥1 P-site |
| Frame-1 coverage | `Psites_coverage_frame1` | Fraction of frame-1 codons with ≥1 P-site |
| Frame-2 coverage | `Psites_coverage_frame2` | Fraction of frame-2 codons with ≥1 P-site |
| **Frame-0 RPKM** | `Psites_frame0_RPKM` | Normalized expression (RPKM) of frame-0 P-sites |
| **f0 vs f1 p-value** | `pval_frame0_vs_frame1` | Wilcoxon test: ORF's frame-0 > frame-1 P-site density |
| **f0 vs f2 p-value** | `pval_frame0_vs_frame2` | Wilcoxon test: ORF's frame-0 > frame-2 P-site density |
| Frame dependence | `dependence_frame1_frame2` | MIC or PCC between frame1/frame2 (optional) |
| **Combined p-value** | `pval_combined` | Stouffer's combined p-value |
| **Adjusted p-value** | `adjusted_pval` | Multiple-testing corrected (default: Benjamini-Hochberg FDR) |

### C. ORF Classification (`ORF_type`)

| Type | Definition |
|------|-----------|
| `annotated` | Overlaps annotated CDS, same stop codon |
| `uORF` | Upstream of CDS, no overlap |
| `dORF` | Downstream of CDS, no overlap |
| `Overlap_uORF` | Upstream, overlaps CDS |
| `Overlap_dORF` | Downstream, overlaps CDS |
| `internal` | Inside CDS, different reading frame |
| `novel` | From non-coding genes or non-coding transcripts of coding genes |

---

## Pipeline Integration

### Module: `RIBOCODE_DETECT` (modules/local/ribocode/detect/main.nf)

**Input channels:**
- `tuple val(meta), path(bam), path(bai)` — Transcriptome BAM from STAR/HISAT2
- `path(gtf)` — Reference GTF annotation
- `path(fasta)` — Reference genome FASTA

**Key configuration:**
- `errorStrategy 'ignore'` — Fails gracefully when periodicity is insufficient
- Uses `RiboCode_onestep` (all-in-one command)
- Strandedness mapping: nf-core `forward→yes`, `reverse→reverse`, `unstranded→no`

**Output channels:**
| Channel | Type | Content |
|---------|------|---------|
| `txt` | `tuple val(meta), path("*.txt")` | All predicted ORFs |
| `collapsed` | `tuple val(meta), path("*_collapsed.txt")` | Non-redundant ORFs |
| `gtf` | `tuple val(meta), path("*.gtf.gz")` | ORF GTF annotations |
| `bed` | `tuple val(meta), path("*.bed.gz")` | ORF BED annotations |
| `results` | `tuple val(meta), path("*.*")` | All outputs (including QC PDFs, config) |

**Parameters:**
| Parameter | Default | Description |
|-----------|---------|-------------|
| `--skip_ribocode` | `false` | Skip RiboCode |
| `--extra_ribocode_args` | `null` | Extra CLI arguments |
| `--ribocode_maxForks` | `1` | Max parallel instances |
| `--ribocode_memory` | `36.GB` | Memory allocation |

**Constraints:**
- Requires STAR or HISAT2 alignment (needs transcriptome BAM)
- Automatically skipped in BAM-input mode
- `--run_ribocode` flag in older pipeline versions

---

## QC Module Design Implications

### Strengths for QC integration:
1. **Direct 3-nt periodicity measurement**: `f0_percent` is a clean, interpretable metric
2. **Per-ORF statistical confidence**: Each ORF has individual p-values and FDR
3. **Read-length aware**: Separates QC by read length, appropriate for Ribo-seq
4. **P-site offset determination**: Automatic, data-driven offset detection

### Metrics to carry into unified QC module:
| Metric | Source | Type |
|--------|--------|------|
| `f0_percent` (per read length) | `_pre_config.txt` | Read-level periodicity |
| `predicted_psite` (per read length) | `_pre_config.txt` | P-site offset |
| `Psites_sum_frame0/1/2` | `_collapsed.txt` | ORF-level abundance |
| `pval_combined` | `_collapsed.txt` | ORF-level significance |
| `adjusted_pval` | `_collapsed.txt` | ORF-level FDR |
| `ORF_type` | `_collapsed.txt` | ORF classification |
| `ORF_length` | `_collapsed.txt` | ORF size |
| Number of ORFs detected (total × type) | `_collapsed.txt` | Yield summary |

### Limitations for QC:
1. **No cross-sample comparison metrics** — per-sample only
2. **No metagene profile summary statistic** beyond individual read-length tables
3. **P-site offset may differ from other tools** (RiboseQC, riboWaltz) — needs harmonization
4. **Depends on annotated CDS** for periodicity assessment — non-coding transcriptome ORFs not evaluated for periodicity at read-level
