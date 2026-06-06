# Ribo-TISH — De Novo ORF Prediction

## Overview

**Ribo-TISH** (Translation Initiation Site Hunter) is a Python-based tool for de novo prediction of translated ORFs from Ribo-seq data. It can optionally incorporate TI-seq (translation initiation sequencing) data for improved start codon identification. It consists of two subcommands: `quality` (QC + P-site offset) and `predict` (ORF detection).

- **Version in pipeline**: 0.2.7
- **Source**: https://github.com/zhpn1024/ribotish (Bioconda)
- **Publication**: Zhang P. et al., *Nucleic Acids Research*, 2017 (doi: 10.1093/nar/gkx452)
- **Language**: Python
- **Key Dependencies**: pysam, numpy, scipy, statsmodels

---

## Analysis Steps

### Step 1: `ribotish quality` — P-site Offset & Read Length QC

**Inputs:**
- Genome-aligned BAM file
- Reference GTF annotation

**Processing:**
1. **Read length distribution**: Counts reads per length around annotated start codons
2. **P-site offset determination**: 
   - Computes aggregate profile of 5′ end distances to annotated start codons
   - Identifies the most enriched distance (peak) for each read length
   - Evaluates 3-nt periodicity by comparing frame-0 enrichment
3. **Quality assessment**: Evaluates whether sufficient periodic signal is present

**Outputs:**
- `*_qual.txt` — Read length distribution table (read_length, count, proportion)
- `*_qual.pdf` — Quality plots (read length distribution + metagene profiles)
- `*.para.py` — **P-site offset parameter file** (Python dict: `offdict = {28: 12, 29: 12, ...}`)

### Step 2: `ribotish predict` — ORF Detection

**Inputs:**
- Ribo-seq BAM file(s)
- Optional: TI-seq BAM file(s)
- Reference FASTA + GTF
- P-site parameter file (`*.para.py`) from quality step
- Candidate ORF definitions (optional)

**Processing:**
1. **ORF candidate enumeration**: Identifies all possible ATG→Stop ORFs across transcripts
2. **P-site shifting**: Applies per-read-length offsets to convert read positions to P-sites
3. **TIS scoring**: For each candidate ORF start codon:
   - **RiboPvalue**: Tests whether P-site density at start codon exceeds background
   - **TISPvalue** (if TI-seq): Tests whether TI-seq reads enrich at start codon
   - **FisherPvalue**: Combined test integrating Ribo-seq + TI-seq signals
4. **Frame assessment**: Evaluates reading frame preference
5. **Significance filtering**: Multiple testing correction (Q-values) for all p-values

**Outputs:**
- `*_pred.txt` — Final predicted ORFs
- `*_all.txt` — All candidate ORFs with scores
- `*_transprofile.py` — Translation profile data

---

## Output Files

### `*_pred.txt` Columns (Predicted ORFs)

| Column | Description |
|--------|-------------|
| `Gid` | Gene ID |
| `Tid` | Transcript ID |
| `Symbol` | Gene symbol |
| `GeneType` | Gene biotype |
| `GenomePos` | Genomic position (chr:start-end:strand) |
| `StartCodon` | Start codon sequence |
| `Strand` | Strand (+/-) |
| `AALen` | Amino acid length |
| `TisType` | TIS classification (e.g., `CDS`, `uORF`, `dORF`, `novel`) |
| `TISGroup` | TIS grouping |
| `TISCounts` | TI-seq read count at TIS |
| **`TISPvalue`** | TI-seq p-value for TIS enrichment |
| **`RiboPvalue`** | Ribo-seq p-value for TIS enrichment |
| `RiboPStatus` | Ribo-seq P-site status |
| **`FisherPvalue`** | Combined Fisher's p-value |
| `TISQvalue` | Adjusted TI-seq q-value |
| `RiboQvalue` | Adjusted Ribo-seq q-value |
| **`FrameQvalue`** | **Adjusted frame preference q-value** |
| `FisherQvalue` | Adjusted combined q-value |

### `*_qual.txt` Columns (Quality)

| Column | Description |
|--------|-------------|
| `read_length` | RPF read length (nt) |
| `count` | Number of reads at this length |
| `proportion` | Fraction of total reads |

### `*.para.py` Content

```python
offdict = {28: 12, 29: 12, 30: 12, 31: 12, 32: 12}
```
Per-read-length P-site offsets (default values shown).

---

## QC Metrics

### A. Read-Level QC (from `quality`)

| Metric | Source | Description |
|--------|--------|-------------|
| Read length distribution | `*_qual.txt` | Counts/proportions per read length |
| **P-site offset per read length** | `*.para.py` | Determined from aggregate start codon profiles |
| Metagene profile quality | `*_qual.pdf` | Visual assessment of periodicity |
| Data sufficiency flag | log output | Whether periodic signal was detected |

### B. ORF-Level QC (from `predict`)

| Metric | Column | Description |
|--------|--------|-------------|
| **RiboPvalue** | `*_pred.txt` | Ribo-seq P-site enrichment at start codon |
| **FrameQvalue** | `*_pred.txt` | Reading frame preference (3-nt periodicity) significance |
| **FisherPvalue** | `*_pred.txt` | Combined Ribo-seq + TI-seq evidence |
| `TISQvalue` / `RiboQvalue` | `*_pred.txt` | Multiple-testing corrected q-values |
| `AALen` | `*_pred.txt` | ORF length in amino acids |
| `TisType` | `*_pred.txt` | ORF classification relative to CDS |

---

## Pipeline Integration

### Module: `RIBOTISH_QUALITY`

**Input:** `tuple val(meta), path(bam), path(bai)` + `tuple val(meta2), path(gtf)`

**Output:**
- `distribution` — `*.txt` (read length distribution)
- `pdf` — `*.pdf` (quality plots)
- `offset` — `*.para.py` (P-site offsets)

### Module: `RIBOTISH_PREDICT`

**Inputs:** Ribo-seq BAM, TI-seq BAM (optional), FASTA, GTF, candidate ORFs, para files

**Output:**
- `predictions` — `*_pred.txt` (final ORF predictions)
- `all` — `*_all.txt` (all candidates)
- `transprofile` — `*_transprofile.py` (profile data)

**Parameters:**
| Parameter | Default | Description |
|-----------|---------|-------------|
| `--skip_ribotish` | `false` | Skip Ribo-TISH |

**Key configuration:**
- Graceful failure for low-signal samples (creates placeholder files)
- Default P-site offsets `{28:12, 29:12, 30:12, 31:12, 32:12}` used when quality step fails

---

## QC Module Design Implications

### Strengths for QC integration:
1. **Separate quality module**: Dedicated QC step before prediction
2. **Multiple p-value types**: RiboPvalue, TISPvalue, FisherPvalue, FrameQvalue — different evidence dimensions
3. **Multiple testing correction**: Q-values for all p-value types
4. **TI-seq integration**: Can incorporate additional experimental evidence
5. **TIS classification**: Annotates ORF type relative to CDS

### Metrics to carry into unified QC module:
| Metric | Source | Type |
|--------|--------|------|
| P-site offsets per read length | `*.para.py` | Read-level offset |
| Read length distribution | `*_qual.txt` | Read-level abundance |
| **`RiboPvalue`** | `*_pred.txt` | ORF-level TIS enrichment |
| **`FrameQvalue`** | `*_pred.txt` | ORF-level frame preference |
| `FisherPvalue` / `FisherQvalue` | `*_pred.txt` | ORF-level combined evidence |
| `TisType` | `*_pred.txt` | ORF classification |
| `AALen` | `*_pred.txt` | ORF size |

### Limitations:
1. **TIS-focused QC**: Periodicity assessment tied to start codon enrichment, not full ORF body
2. **No explicit phase score**: `FrameQvalue` captures frame preference but less directly than Ribotricer's phase score
3. **Quality step can fail**: Low-depth data may not yield valid offsets; pipeline falls back to hardcoded defaults
4. **No ORF-level coverage metrics**: Doesn't report per-ORF read counts or density
5. **TI-seq dependent for best performance**: Ribo-seq-only mode has reduced sensitivity
