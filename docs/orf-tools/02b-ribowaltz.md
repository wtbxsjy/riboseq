# riboWaltz — P-site Offset & Diagnostic Analysis

## Overview

**riboWaltz** is an R/Bioconductor package specialized in calculating optimal P-site offsets, performing diagnostic analyses, and generating visualizations for ribosome profiling data. In the riboseq pipeline, it serves as both a complementary QC tool to RiboseQC and a **P-site offset fallback** when RiboseQC data is unavailable.

- **Version in pipeline**: 2.0+ (from GitHub, patched for Bioc 3.20+)
- **Source**: https://github.com/LabTranslationalArchitectomics/riboWaltz
- **Publication**: Lauria F. et al., *PLoS Computational Biology*, 2018 (doi: 10.1371/journal.pcbi.1006169)
- **Language**: R
- **Key Dependencies**: data.table, ggplot2, GenomicFeatures, GenomicAlignments, txdbmaker (Bioc 3.20+)

---

## Analysis Steps

riboWaltz performs a sequential analysis pipeline:

### Step 1: Annotation Loading

**Inputs:**
- GTF annotation file (GENCODE/Ensembl)
- Alternatively: pre-built cached annotation RDS

**Processing:**
- Parses GTF to extract per-transcript lengths: total (`l_tr`), 5′ UTR (`l_utr5`), CDS (`l_cds`), 3′ UTR (`l_utr3`)
- Uses `txdbmaker::makeTxDbFromGFF()` (Bioc 3.20+ compatible, replacing defunct `GenomicFeatures::makeTxDbFromGFF`)
- **Strips transcript ID version suffixes** (e.g., `.11`) so BAM reads match GTF transcripts

**Output:** `annotation` data table:
```
| transcript | l_tr | l_utr5 | l_cds | l_utr3 |
```

### Step 2: BAM Loading

**Inputs:**
- Transcriptome-aligned BAM file
- Annotation data table

**Processing:**
- Reads alignments via `GenomicAlignments::readGAlignments()`
- **Strips transcript ID versions from BAM seqlevels** for matching with annotation
- Filters to annotated transcripts only
- Computes `cds_start` and `cds_stop` for each read from annotation

### Step 3: Read Length Filtering

**Default:** `c(28, 29, 30)` (configurable via `--ribowaltz_read_lengths`)

**Methods:**
- **Custom mode**: User-specified read lengths (pipeline default)
- **Periodicity mode**: Automatically selects read lengths with periodicity above threshold

### Step 4: P-site Offset Calculation (`psite()`)

**This is riboWaltz's core algorithm — a two-step correction procedure:**

#### Phase 1: Temporary Offset Determination
1. Groups reads by length around the reference codon (start codon, `start = TRUE`)
2. For each read length bin, computes 5′ and 3′ end density profiles
3. Identifies the global maximum in each profile → **temporary 5′ and 3′ offsets**
4. The `flanking` parameter (default 6 nt) excludes reads too close to the reference codon

#### Phase 2: Offset Correction
1. Automatically selects the optimal read extremity (5′ or 3′) by comparing profile patterns — selects the extremity whose distance to the reference codon is most **stable across read lengths**
2. Identifies the **most frequent temporary offset** among the predominant (highest-signal) read length bins
3. For smaller bins: adjusts the temporary offset to the local maximum whose distance from the reference codon is **closest to the reference offset**
4. This maximizes offset coherence across read length populations

**Output — `psite_offset` data table columns:**

| Column | Description |
|--------|-------------|
| `length` | Read length (nt) |
| `total_percentage` | Percentage of total reads at this length |
| `start_percentage` | Percentage of reads at this length aligning on the start codon |
| `around_start` | Whether offset is computed around start codon (T/F) |
| `offset_from_5` | Distance from 5′ end to P-site (before correction) |
| `offset_from_3` | Distance from 3′ end to P-site (before correction) |
| **`corrected_offset_from_5`** | **Final P-site offset from 5′ end (after correction)** |
| `corrected_offset_from_3` | Final P-site offset from 3′ end (after correction) |
| `sample` | Sample identifier |

### Step 5: P-site Assignment (`psite_info()`)

Assigns P-site positions to each read:
- `psite` — P-site position relative to transcript start (1-based)
- `psite_from_start` — Distance from annotated start codon
- `psite_from_stop` — Distance from annotated stop codon
- `psite_region` — Transcript region (5utr, cds, 3utr)

### Step 6: QC Analyses

| Analysis | Function | Output |
|----------|----------|--------|
| **CDS Coverage** | `cds_coverage()` | `*_cds_coverage.tsv` — P-site counts per transcript CDS |
| **Frame Distribution** | `frame_psite()` | `*_frame_distribution.tsv` — % P-sites in frames 0/1/2 per region |
| **Region Distribution** | (derived from `psite_info`) | `*_region_distribution.tsv` — P-site counts per transcript region |
| **Read Length Distribution** | `rlength_distr()` | `*_read_length_distribution.pdf` |
| **Metagene Profiles** | `metaprofile_psite()` | `*_metaprofile.pdf` |
| **Meta Heatmaps** | `metaheatmap_psite()` | `*_metaheatmap.pdf` |
| **P-site Offset Plots** | `psite(plot=TRUE)` | Per-read-length offset metaprofiles in `*_ribowaltz_plots/` |

---

## QC Metrics

### A. Read-Level QC (P-site Offset)

| Metric | Source Column | Description |
|--------|--------------|-------------|
| **Corrected P-site offset** | `corrected_offset_from_5` in `*_psite_offset.tsv` | Final P-site position from 5′ end (primary output) |
| Read proportion | `total_percentage` | Contribution of each read length |
| Start codon enrichment | `start_percentage` | Fraction of reads at start codon |
| Optimal extremity | `best_offset.txt` | Whether 5′ or 3′ end was used for correction |

**Key advantage over RiboseQC/RiboCode:** riboWaltz's two-step correction algorithm maximizes P-site offset coherence across read lengths, producing more accurate offsets when multiple read length populations are present.

### B. Transcript-Level QC

| Metric | Source | Description |
|--------|--------|-------------|
| **Frame 0/1/2 distribution** | `*_frame_distribution.tsv` | **3-nt periodicity**: % of P-sites in each reading frame, per transcript region (5′UTR/CDS/3′UTR) |
| CDS coverage | `*_cds_coverage.tsv` | P-site counts per transcript CDS |
| Region distribution | `*_region_distribution.tsv` | P-sites per region (5utr/cds/3utr) |

### C. Visual Diagnostics

| Plot | Content | QC Purpose |
|------|---------|------------|
| Read length distribution | Bar plot of reads per length | Identify major RPF populations |
| P-site offset metaprofiles | Per-length 5′/3′ end density around start codon | Verify offset determination accuracy |
| Metagene profile | P-site density across CDS (±UTR flanks) | Visualize 3-nt periodicity pattern |
| Meta heatmap | Heatmap of P-site density around start/stop | Compare periodicity across samples |
| Frame distribution barplot | % P-sites in frames 0/1/2 by region | Quantify periodicity enrichment in CDS |

---

## Pipeline Integration

### Module: `RIBOWALTZ`

**Runs before RiboseQC** — its P-site offsets serve as fallback for `EXTRACT_RL_CUTOFF` → `PREPARE_FOR_ORFQUANT_CORRECTED` → ORFquant chain.

**Input channels:**
- `tuple val(meta), path(bam), path(bai)` — Transcriptome BAM (STAR/HISAT2 alignment mode) or genome BAM (BAM-input mode)
- `path(gtf)` — Reference GTF annotation
- `path(fasta)` — Reference genome FASTA

**Output channels:**

| Channel | Type | Content |
|---------|------|---------|
| `psite_offset` | `tuple val(meta), path("*_psite_offset.tsv")` | P-site offsets (tab-separated) |
| `psite_offset_txt` | `tuple val(meta), path("*_psite_offset.txt")` | Best offset + extremity summary |
| `cds_coverage` | `tuple val(meta), path("*_cds_coverage.tsv")` | CDS P-site counts |
| `codon_usage` | `tuple val(meta), path("*_codon_usage.tsv")` | Codon usage (optional) |
| `frame_distribution` | `tuple val(meta), path("*_frame_distribution.tsv")` | Frame 0/1/2 proportions |
| `region_distribution` | `tuple val(meta), path("*_region_distribution.tsv")` | P-sites per region |
| `plots` | `tuple val(meta), path("*_ribowaltz_plots")` | All diagnostic plots |

**P-site offset fallback chain:**

```
1. RiboseQC P_sites_calcs valid → use RiboseQC offset (unchanged)
2. RiboseQC data empty + riboWaltz available → use riboWaltz corrected_offset_from_5 per read length
3. Neither available → hardcoded defaults (28–32 nt → offset 12)
```

**Parameters:**

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--skip_ribowaltz` | `false` | Skip riboWaltz |
| `--ribowaltz_read_lengths` | `[28, 29, 30]` | Read lengths to analyze |
| `--extra_ribowaltz_args` | `null` | Extra arguments |
| `--ribowaltz_container` | (auto) | Custom container path |

**Constraints:**
- Uses **transcriptome BAM** in alignment mode (STAR with `--quantMode TranscriptomeSAM`)
- Falls back to **genome BAM** in BAM-input mode
- Bioc 3.20+ compatibility: uses `txdbmaker::makeTxDbFromGFF()` instead of defunct `GenomicFeatures::makeTxDbFromGFF()`
- Transcript ID version stripping for 100% BAM-GTF match rate
- Annotation caching (RDS) to avoid rebuilding TxDb per sample (~3 min → <1 sec)

---

## QC Module Design Implications

### Strengths for QC integration:
1. **Best P-site offset algorithm**: Two-step correction maximizes inter-length coherence — more accurate than simple peak-picking used by other tools
2. **Automatic extremity selection**: Data-driven choice of 5′ vs 3′ end for offset determination
3. **Comprehensive frame analysis**: Per-region frame distributions (5′UTR/CDS/3′UTR) with quantitative proportions
4. **CDS coverage quantification**: Per-transcript P-site counts usable for expression analysis
5. **Rich visualization**: Multiple plot types for periodicity diagnosis
6. **Fallback role**: Ensures P-site offsets are available even when RiboseQC fails

### Metrics to carry into unified QC module:

| Metric | Source | Type |
|--------|--------|------|
| `corrected_offset_from_5` per read length | `*_psite_offset.tsv` | Read-level P-site offset |
| `total_percentage` / `start_percentage` | `*_psite_offset.tsv` | Read abundance + start enrichment |
| Frame 0/1/2 proportions (CDS) | `*_frame_distribution.tsv` | **3-nt periodicity** (per region) |
| CDS P-site counts | `*_cds_coverage.tsv` | Transcript-level expression |
| P-site region distribution | `*_region_distribution.tsv` | Enrichment in CDS vs UTR |
| Optimal extremity | `*_psite_offset.txt` | 5′ vs 3′ end preference |

### Comparison with RiboseQC:

| Feature | riboWaltz | RiboseQC |
|---------|-----------|----------|
| P-site algorithm | Two-step coherence correction | Peak-based with frame distribution |
| Requires | Transcriptome BAM | Genome BAM |
| Frame analysis | Per-region proportions | Per-transcript distributions |
| Codon usage | Yes (needs FASTA) | Yes (built-in) |
| Biotype analysis | No | Yes |
| HTML report | No (PDF plots only) | Yes |
| ORFquant input | No | Yes (`*_for_ORFquant`) |

### Limitations:
1. **Transcriptome BAM required**: In BAM-input mode, only genome BAM available — may reduce accuracy
2. **No biotype/contamination analysis**: Cannot assess rRNA/mtDNA/plastid contamination
3. **No HTML report**: Only PDF plots; less convenient for multi-sample comparison
4. **No ORF-level metrics**: Operates at read/transcript level, not individual ORF level
5. **Per-read-length offset only**: Does not integrate offsets into ORF-level scores
