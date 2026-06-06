# ORFquant — Splice-Aware ORF Quantification

## Overview

**ORFquant** is an R/Bioconductor package for detecting and quantifying ORF translation on complex transcriptomes using Ribo-seq data. It performs splice-aware, single-ORF-level quantification considering multiple transcript isoforms per gene. It depends on RiboseQC for input preparation.

- **Version in pipeline**: 1.02 (patched)
- **Source**: https://github.com/ohlerlab/ORFquant
- **Publication**: Calviello L. et al., *Nature Structural & Molecular Biology*, 2020 (doi: 10.1038/s41594-019-0299-6)
- **Language**: R (Bioconductor)
- **Key Dependencies**: GenomicFeatures, GenomicAlignments, BSgenome, rtracklayer, multitaper (DPSS), Biostrings

---

## Analysis Steps

ORFquant operates in three main phases via `run_ORFquant()`:

### Phase 1: Transcript Filtering & ORF Detection (`detect_translated_orfs`)

**Inputs:**
- `*_for_ORFquant` file (P-site positions from RiboseQC)
- `*_Rannot` annotation file
- n_cores for parallelization

**Processing:**
1. **Transcript filtering**: Removes low-expression transcripts and those with insufficient P-site signal
2. **P-site aggregation**: Combines P-site positions across replicates/chunks
3. **De-novo ORF finding** (`get_orfs()`): For each transcript:
   - Translates sequence in all 3 reading frames
   - Identifies all ATG→Stop codon pairs
   - Considers alternative start codons if specified

### Phase 2: ORF P-site Signal Quantification

For each detected ORF:
1. Extracts P-site density profile along the ORF
2. Computes **DPSS (Slepian) multitaper spectral analysis** for statistical testing
3. Calculates periodicity/translation significance via `calc_orf_pval()`
4. Filters ORFs by significance threshold

### Phase 3: ORF Classification & Annotation

Classifies each ORF at three levels:
- **`ORF_category_Gen`**: Genomic-level classification (relative to annotated CDS)
- **`ORF_category_Tx`**: Transcript-level classification (per isoform)
- **`ORF_category_Tx_compatible`**: Best classification across all compatible transcripts

---

## Output Files

| File | Content |
|------|---------|
| `*_final_ORFquant_results` | Complete ORFquant results (RData) |
| `*_Detected_ORFs.gtf.gz` | ORF annotations in GTF format |
| `*_Protein_sequences.fasta` | Predicted ORF protein sequences |
| `*_tmp_ORFquant_results` | Intermediate results (optional) |
| `*_ORFquant_plots_RData` | Plot data (optional) |
| `*_plots/` | Visualization plots (optional) |

---

## QC Metrics

### A. ORF-Level QC

From `*_Detected_ORFs.gtf` attributes and `*_final_ORFquant_results`:

| Metric | Description |
|--------|-------------|
| `ORF_category_Gen` | Genomic classification (e.g., `CDS`, `uORF`, `dORF`, `overlap_uORF`, `ncRNA`) |
| `ORF_category_Tx` | Transcript-level classification |
| `ORF_category_Tx_compatible` | Best classification across isoforms |
| `ORF_id` | Unique ORF identifier |
| `gene_id` / `gene_name` | Associated gene |
| `transcript_id` | Associated transcript |
| `p_value` | Statistical significance of translation (from DPSS multitaper test) |

### B. Signal-Level QC (internal)

| Metric | Description |
|--------|-------------|
| P-site coverage per ORF | Total P-sites within ORF boundaries |
| DPSS spectral power | Measures 3-nt periodicity in P-site signal |
| Translation score | Combined periodicity + abundance metric |

---

## Pipeline Integration

### Module: `ORFQUANT_RUN`

**Input channels:**
- `tuple val(meta), path(for_orfquant)` — `*_for_ORFquant` from RiboseQC
- `path(annotation)` — `*_Rannot` annotation
- `path(fasta)` — Genome FASTA
- `path(orfquant_pkg)` — Optional pre-downloaded ORFquant R package

**Key configuration:**
- Requires RiboseQC output (`*_for_ORFquant`); auto-skipped if RiboseQC is skipped
- **Custom patched container**: Uses patched ORFquant to resolve BiocGenerics namespace conflicts
- DPSS caching for performance optimization
- Error handling: graceful failure for low-signal samples

**Output channels:**
| Channel | Type | Content |
|---------|------|---------|
| `results` | `tuple val(meta), path("*_final_ORFquant_results")` | Complete results |
| `gtf` | `tuple val(meta), path("*_Detected_ORFs.gtf.gz")` | ORF GTF |
| `proteins` | `tuple val(meta), path("*_Protein_sequences.fasta")` | Protein sequences |
| `tmp_results` | `tuple val(meta), path("*_tmp_ORFquant_results")` | Intermediate (optional) |
| `plots_data` | `tuple val(meta), path("*_ORFquant_plots_RData")` | Plot data (optional) |
| `plots_dir` | `tuple val(meta), path("*_plots")` | Visualization plots (optional) |

**Parameters:**
| Parameter | Default | Description |
|-----------|---------|-------------|
| `--skip_orfquant` | `false` | Skip ORFquant |
| `--orfquant_container` | (auto) | Custom ORFquant container path |
| `--orfquant_pkg` | (auto) | Pre-downloaded ORFquant package |

---

## QC Module Design Implications

### Strengths for QC integration:
1. **Splice-aware**: Correctly handles multi-exon ORFs across transcript isoforms
2. **Multi-level classification**: Genomic + transcript + best-isoform annotation
3. **Statistical testing**: DPSS multitaper provides rigorous periodicity evaluation
4. **Direct P-site integration**: Uses RiboseQC P-sites, ensuring consistency

### Metrics to carry into unified QC module:
| Metric | Source | Type |
|--------|--------|------|
| `ORF_category_Gen` / `_Tx` / `_Tx_compatible` | GTF attributes | Classification |
| `ORF_id` | GTF | ORF identifier |
| `gene_id` / `transcript_id` | GTF | Genomic context |
| Number of ORFs detected | Output count | Yield summary |
| ORF length distribution | GTF | Size distribution |

### Limitations:
1. **Heavy dependency on RiboseQC**: Cannot run independently
2. **No per-read-length metrics**: Aggregates all read lengths internally
3. **No P-site offset output**: Relies entirely on RiboseQC for offset calculation
4. **Computationally expensive**: DPSS multitaper analysis is CPU-intensive
5. **Binary RData results**: Main results are in RData format, not tabular
