# GENCODE ORF Format Converters

## Overview

This directory contains standalone Python scripts that convert ORF predictions from various Ribo-seq tools to gencode-riboseqORFs compatible format.

## Available Converters

### 1. ribotish_to_gencode.py ✅

**Status**: Completed and tested

Converts Ribo-TISH predictions to GENCODE format.

```bash
python3 bin/ribotish_to_gencode.py \
    --predict ribotish_predict.txt \
    --fasta genome.fa \
    --study_id SAMPLE1 \
    --output_prefix output \
    --min_length 16
```

**Test data**: `test_data/ribotish_to_gencode/`

### 2. ribotricer_to_gencode.py ✅

**Status**: Completed and tested

Converts Ribotricer predictions to GENCODE format.

```bash
python3 bin/ribotricer_to_gencode.py \
    --tsv sample_translating_ORFs.tsv \
    --fasta genome.fa \
    --study_id SAMPLE1 \
    --output_prefix output \
    --min_length 16 \
    --min_phase_score 0.5
```

**Test data**: `test_data/ribotricer_to_gencode/`

### 3. ribocode_to_gencode.py ⏳

**Status**: Planned

Converts RiboCode predictions to GENCODE format.

**Input**: RiboCode GTF + FASTA output

### 4. rpbp_to_gencode.py ⏳

**Status**: Planned

Converts rp-bp predictions to GENCODE format.

**Input**: rp-bp BED + FASTA output

### 5. orfquant_to_gencode.py ⏳

**Status**: Planned

Converts ORFquant predictions to GENCODE format.

**Input**: ORFquant GTF/BED output

## Output Format

All converters produce two files:

### 1. FASTA File (`{prefix}.gencode.fa`)

```
>{GENE}_{START}_{LENGTH}aa--{STUDY_ID}
SEQUENCE*
```

- Header format: `{ORF_NAME}--{STUDY_ID}`
- Sequence must end with stop codon `*`

### 2. BED File (`{prefix}.gencode.bed`)

```
chr  start  end  ORF_NAME  STUDY_ID  strand
```

- **Important**: Uses 1-based coordinates (gencode-riboseqORFs requirement)
- Tab-delimited, 6 columns
- No header line

## Common Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--study_id` | Sample/study identifier | Required |
| `--output_prefix` | Output file prefix | Required |
| `--fasta` | Genome FASTA for sequence extraction | Required |
| `--min_length` | Minimum ORF length (aa) | 16 |

## Dependencies

### Required
- Python 3.9+
- Biopython 1.79+

### Optional (for real sequence extraction)
- pyfaidx 0.7.2+

Install:
```bash
pip install biopython pyfaidx
```

## Testing

Each converter has test data in `test_data/{converter_name}/`:

```bash
# Test Ribo-TISH converter
cd test_data/ribotish_to_gencode
python3 ../../bin/ribotish_to_gencode.py \
    --predict ribotish_predict.txt \
    --fasta test_genome.fa \
    --study_id TEST \
    --output_prefix test_output

# Test Ribotricer converter
cd test_data/ribotricer_to_gencode
python3 ../../bin/ribotricer_to_gencode.py \
    --tsv ribotricer_translating_ORFs.tsv \
    --fasta test_genome.fa \
    --study_id TEST \
    --output_prefix test_output
```

## Coordinate Systems

**Critical Note**: Different tools use different coordinate systems.

| Tool | Input Coordinates | Output BED | Notes |
|------|------------------|------------|-------|
| Ribo-TISH | 1-based (GenomePos) | 1-based | Already in correct format |
| Ribotricer | 1-based (genomic) | 1-based | start_codon is 1-based |
| RiboCode | 0-based (BED) | 1-based | **Need conversion** (+1 to start) |
| rp-bp | 0-based (BED) | 1-based | **Need conversion** (+1 to start) |
| ORFquant | GTF (1-based) | 1-based | Already in correct format |

## Integration Workflow

```
┌─────────────────┐
│  ORF Predictor  │
│ (Ribo-TISH/etc) │
└────────┬────────┘
         │
         ▼
┌─────────────────┐
│    Converter    │
│  (this script)  │
└────────┬────────┘
         │
         ▼
┌─────────────────┐
│  GENCODE Format │
│  (FA + BED)     │
└────────┬────────┘
         │
         ▼
┌─────────────────┐
│ gencode-        │
│ riboseqORFs     │
└────────┬────────┘
         │
         ▼
┌─────────────────┐
│ Unified ORF     │
│ Annotations     │
└─────────────────┘
```

## Next Steps

After running converters, use gencode-riboseqORFs:

```bash
# 1. Download Ensembl annotation
bash scripts/retrieve_ensembl_data.sh 110 GRCh38

# 2. Merge all converter outputs
cat sample1.gencode.fa sample2.gencode.fa > all_orfs.fa
cat sample1.gencode.bed sample2.gencode.bed > all_orfs.bed

# 3. Run gencode-riboseqORFs mapper
python3 ORF_mapper_to_GENCODE_v1.1.py \
    -d Ens110/ \
    -f all_orfs.fa \
    -b all_orfs.bed \
    -o project_name \
    -l 16 \
    -c 0.9 \
    -m longest_string
```

## Development Status

| Component | Status | Priority |
|-----------|--------|----------|
| ribotish_to_gencode.py | ✅ Complete | High |
| ribotricer_to_gencode.py | ✅ Complete | High |
| Test data | ✅ Complete | High |
| Documentation | ✅ Complete | High |
| ribocode_to_gencode.py | ⏳ Pending | Medium |
| rpbp_to_gencode.py | ⏳ Pending | Medium |
| orfquant_to_gencode.py | ⏳ Pending | Medium |
| Nextflow module integration | ⏳ Pending | High |

## See Also

- [GENCODE Integration Project Summary](../GENCODE_INTEGRATION_PROJECT_SUMMARY.md)
- [Implementation Plan](../docs/GENCODE_INTEGRATION_IMPLEMENTATION_PLAN.md)
- [gencode-riboseqORFs GitHub](https://github.com/jorruior/gencode-riboseqORFs)

---

**Last Updated**: 2026-01-17
**Maintainer**: nf-core/riboseq team
