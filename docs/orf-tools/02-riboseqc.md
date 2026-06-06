# RiboseQC — Comprehensive Ribo-seq QC Tool

## Overview

**RiboseQC** (Ribo-seQC) is an R/Bioconductor package for comprehensive quality control analysis of Ribo-seq data. It is the central QC engine in the pipeline, running both before and after sORF filtering. It also generates the `*_for_ORFquant` input files required by ORFquant.

- **Version in pipeline**: 1.1
- **Source**: https://github.com/lcalviell/Ribo-seQC
- **Publication**: Calviello L. et al. (in preparation)
- **Language**: R (Bioconductor)
- **Key Dependencies**: GenomicFeatures, GenomicAlignments, Rsamtools, BSgenome, rtracklayer, ggplot2

---

## Analysis Steps

RiboseQC performs a single comprehensive analysis via `RiboseQC_analysis()`:

### Step 1: Annotation Loading

**Inputs:**
- `*_Rannot` file (from `prepare_annotation_files`)
- BAM file (genome-aligned)
- Optional: genome FASTA (FaFile)

### Step 2: Read Statistics Collection

Processes BAM in chunks (default 5M reads/chunk) to minimize RAM:

1. **Read length distribution** (`rld`, `rld_unq`): Per compartment (nuclear, mitochondrial, chloroplast)
2. **Biotype mapping** (`positions`, `positions_unq`): 5′ end counts on CDS, 5′UTRs, 3′UTRs, ncRNAs, introns, intergenic
3. **Gene-level counts** (`counts_cds_genes`, `counts_all_genes`): RPKM and TPM values for CDS and all genes
4. **Read length × compartment summary** (`reads_summary`)

### Step 3: Five-Prime Profile Construction

For each read length and compartment:
- **Binned profiles** (`five_prime_bins`): Signal over 50 5′UTR bins + 100 CDS bins + 50 3′UTR bins, using one representative transcript per gene
- **Subcodon profiles** (`five_prime_subcodon`): High-resolution signal around TSS, start codon, mid-ORF, stop codon, and TES

### Step 4: P-site Cutoff Calculation (`calc_cutoffs_from_profiles`)

**This is where 3-nt periodicity is assessed:**

For each read length on each compartment:
1. Builds frame-0/1/2 coverage distributions across all CDS transcripts
2. Computes per-transcript frame proportions (f0/(f0+f1+f2))
3. Generates frame distribution boxplots
4. Calculates `analysis_frame_cutoff` statistics (frame proportions per read length per transcript)

**Read length selection methods** (`choose_readlengths`):
- `max_coverage` (default): Selects read lengths with highest total coverage
- `max_inframe`: Selects read lengths with highest frame-0 proportion
- `all`: Uses all read lengths regardless of periodicity

### Step 5: P-site Position Calculation

For each selected read length:
1. Extracts P-site positions from 5′ end mapping around start codons
2. Identifies the most common offset (peak) for each read length
3. Generates `summary_P_sites` with per-read-length P-site offsets and usage statistics

**Three P-site variants computed:**
- `P_sites_all`: From all mapped reads
- `P_sites_all_uniq`: From uniquely mapping reads (MAPQ > 50)
- `P_sites_uniq_mm`: Uniquely mapping reads with mismatches

### Step 6: P-site Profile Construction

- **Binned P-site profiles** (`P_sites_bins`): P-site signal over binned transcript regions
- **Subcodon P-site profiles** (`P_sites_subcodon`): P-site signal at high resolution around TSS/start/stop/TES
- **Codon occupancy** (`Codon_counts`, `P_sites_percodon`, `P_sites_percodon_ratio`): P-site counts per codon at ORF start/middle/end

### Step 7: Additional Analyses

- **Top mapping sequences** (`sequence_analysis`): Top 50 mapping locations with DNA sequence, read counts, feature annotation
- **Junction statistics** (`junctions`): Read mapping on annotated splice junctions
- **Coverage bedgraphs**: Normalized (sum to 1M) bedgraph files for plus/minus strands

---

## Output Files

| File | Content |
|------|---------|
| `*_results_RiboseQC` | Core results for HTML report (RData) |
| `*_results_RiboseQC_all` | Complete analysis results (RData) |
| `*_for_ORFquant` | P-site positions for ORFquant input (RData) |
| `*_P_sites_calcs` | **P-site offset table** per read length with stats |
| `*_coverage_plus.bedgraph` | Normalized coverage, plus strand |
| `*_coverage_minus.bedgraph` | Normalized coverage, minus strand |
| `*_P_sites_plus.bedgraph` | P-site positions, plus strand |
| `*_P_sites_minus.bedgraph` | P-site positions, minus strand |
| `*_junctions` | Junction mapping statistics |
| `*_ggribo.tsv` | P-sites in ggRibo format (count, chr, pos, strand) |
| `*_RiboseQC_report.html` | Interactive HTML QC report (optional) |

---

## QC Metrics

### A. Read-Level QC

| Metric | Source | Description |
|--------|--------|-------------|
| Read length distribution | `rld` | Counts per read length per compartment |
| Biotype distribution | `positions` | 5′ end distribution across CDS/UTR/ncRNA/intron/intergenic |
| Read counts per gene | `counts_cds_genes`, `counts_all_genes` | RPKM + TPM per gene |
| Mapping statistics | `reads_summary` | Reads per biotype × read length × organelle |

### B. Periodicity QC (P-site Analysis)

| Metric | Source | Description |
|--------|--------|-------------|
| **Frame-0/1/2 proportions** | `analysis_frame_cutoff` | Per-transcript frame distribution for each read length |
| **P-site offset per read length** | `P_sites_calcs` / `summary_P_sites` | Most probable P-site distance from 5′ end |
| **P-site usage fraction** | `summary_P_sites` | Fraction of total reads used after read length selection |
| **Selected read lengths** | `results_choice` | Read lengths passing coverage/periodicity filters |
| Codon occupancy | `P_sites_percodon_ratio` | P-site enrichment per codon position |

### C. Coverage QC

| Metric | Source | Description |
|--------|--------|-------------|
| CDS coverage (5′ end) | `five_prime_bins` | Metagene profile across CDS |
| CDS coverage (P-site) | `P_sites_bins` | P-site metagene profile across CDS |
| Strand-specific coverage | `*_coverage_*.bedgraph` | Genome-wide coverage tracks |

---

## Pipeline Integration

### Module: `RIBOSEQC_ANALYSIS`

**Input channels:**
- `tuple val(meta), path(bam), path(bai)` — Genome-aligned BAM
- `path(annotation)` — `*_Rannot` annotation file
- `path(fasta)` — Genome FASTA

**Key configuration:**
- Runs twice: pre-filter (unfiltered BAM) and post-filter (sORF-filtered BAM)
- `fast_mode = TRUE` (default): Uses top 500 genes for faster profiling
- Read length selection: `max_coverage` method
- Error handling: Creates placeholder files for low-signal samples

**Output channels:**
| Channel | Type | Content |
|---------|------|---------|
| `results` | `tuple val(meta), path("*_results_RiboseQC")` | Core results |
| `results_all` | `tuple val(meta), path("*_results_RiboseQC_all")` | Full results (optional) |
| `orfquant` | `tuple val(meta), path("*_for_ORFquant")` | ORFquant input |
| `coverage` | `tuple val(meta), path("*_coverage_*.bedgraph")` | Coverage tracks |
| `psites_bedgraph` | `tuple val(meta), path("*_P_sites_*.bedgraph")` | P-site tracks |
| `psites_calcs` | `tuple val(meta), path("*_P_sites_calcs")` | P-site offsets |
| `junctions` | `tuple val(meta), path("*_junctions")` | Junction stats |
| `ggribo` | `tuple val(meta), path("*_ggribo.tsv")` | ggRibo format P-sites |

**Parameters:**
| Parameter | Default | Description |
|-----------|---------|-------------|
| `--skip_riboseqc` | `false` | Skip RiboseQC entirely |
| `--extra_riboseqc_args` | `null` | Extra arguments |

---

## QC Module Design Implications

### Strengths for QC integration:
1. **Most comprehensive periodicity analysis**: Frame distributions per transcript × read length, boxplot distributions
2. **Codon occupancy analysis**: Position-specific P-site enrichment at start/middle/end of ORFs
3. **Multi-compartment aware**: Separates nuclear/mitochondrial/plastid signals
4. **P-site offset with usage fraction**: Knows what proportion of total reads are used
5. **Generates ORFquant input**: Bridges QC and ORF prediction

### Metrics to carry into unified QC module:
| Metric | Source | Type |
|--------|--------|------|
| `P_sites_calcs` table | `*_P_sites_calcs` | Per-read-length P-site offsets |
| Frame distribution (f0/f1/f2) | `analysis_frame_cutoff` | 3-nt periodicity per transcript |
| Selected read lengths | `results_choice` | QC-passing read lengths |
| Read proportion used | `summary_P_sites` | P-site usage efficiency |
| Biotype distribution | `positions` | Contamination assessment |
| CDS coverage profile | `P_sites_bins` | Metagene profile shape |
| RPKM/TPM | `counts_cds_genes` | Expression quantification |

### Limitations:
1. **No per-ORF QC**: Reports aggregate profiles, not individual ORF metrics
2. **No ORF classification**: Does not classify ORF types (uORF, dORF, etc.)
3. **P-site offset may differ from other tools**: Needs harmonization with RiboCode/riboWaltz offsets
4. **Low-signal failure**: Can fail on low-depth samples (pipeline handles gracefully)
5. **Version**: Pipeline uses v1.1 (older Bioconductor 3.6/R 3.6), newer versions available
