# ribotricer_to_gencode.py - Test Documentation

## Overview

The `ribotricer_to_gencode.py` script converts Ribotricer ORF predictions to gencode-riboseqORFs compatible format.

## Script Location

```bash
bin/ribotricer_to_gencode.py
```

## Usage

```bash
python3 bin/ribotricer_to_gencode.py \
    --tsv <sample_translating_ORFs.tsv> \
    --fasta <genome.fa> \
    --study_id <STUDY_ID> \
    --output_prefix <output_prefix> \
    [--min_length 16] \
    [--min_phase_score 0.5]
```

## Parameters

| Parameter | Required | Default | Description |
|-----------|----------|---------|-------------|
| `--tsv` | Yes | - | Ribotricer translating_ORFs.tsv output file |
| `--fasta` | Yes | - | Genome FASTA file for sequence extraction |
| `--study_id` | Yes | - | Study identifier (e.g., sample name, dataset ID) |
| `--output_prefix` | Yes | - | Prefix for output files |
| `--min_length` | No | 16 | Minimum ORF length in amino acids |
| `--min_phase_score` | No | 0.5 | Minimum phase score for quality filtering |

## Input Format: Ribotricer TSV

The input file `*_translating_ORFs.tsv` contains 18 columns:

### Required Columns

1. **ORF_ID** - Unique identifier (e.g., `ENST00000456328.2_1`)
2. **ORF_type** - ORF classification (see categories below)
3. **status** - Translation status (`translating` or `non-translating`)
4. **phase_score** - Quality metric for ribosome periodicity (0-1)
5. **read_count** - Total reads covering the ORF
6. **length** - ORF length in nucleotides
7. **valid_codons** - Number of valid codons
8. **valid_codons_ratio** - Ratio of valid codons
9. **read_density** - Read coverage density
10. **transcript_id** - Transcript identifier (e.g., `ENST00000456328.2`)
11. **transcript_type** - Transcript biotype
12. **gene_id** - Gene identifier (e.g., `ENSG00000223972.5`)
13. **gene_name** - Gene symbol (e.g., `DDX11L1`)
14. **gene_type** - Gene biotype
15. **chrom** - Chromosome name
16. **strand** - Strand orientation (`+` or `-`)
17. **start_codon** - Genomic position of start codon (1-based)
18. **profile** - Coverage profile array

### ORF Type Categories

Ribotricer classifies ORFs into 8 types:

- **annotated**: CDS annotated in GTF
- **super_uORF**: Upstream ORF, not overlapping any CDS
- **super_dORF**: Downstream ORF, not overlapping any CDS
- **uORF**: Upstream ORF, not overlapping main CDS
- **dORF**: Downstream ORF, not overlapping main CDS
- **overlap_uORF**: Upstream ORF overlapping main CDS
- **overlap_dORF**: Downstream ORF overlapping main CDS
- **novel**: ORF in non-coding genes/transcripts

### Example Input

```tsv
ORF_ID	ORF_type	status	phase_score	read_count	length	valid_codons	valid_codons_ratio	read_density	transcript_id	transcript_type	gene_id	gene_name	gene_type	chrom	strand	start_codon	profile
ENST00000456328.2_1	annotated	translating	0.95	1500	333	111	1.0	4.5	ENST00000456328.2	protein_coding	ENSG00000223972.5	DDX11L1	protein_coding	chr1	+	100000	[15,18,20,...]
ENST00000515242.2_1	uORF	translating	0.88	800	150	50	0.98	5.3	ENST00000515242.2	protein_coding	ENSG00000227232.5	WASH7P	protein_coding	chr1	+	200000	[8,10,12,...]
```

## Output Formats

### 1. FASTA Output (`<prefix>.gencode.fa`)

Format:
```
>{GENE_NAME}_{START}_{LENGTH}aa--{STUDY_ID}
SEQUENCE*
```

Example:
```
>DDX11L1_100000_111aa--TEST_RIBOTRICER
MAAGTLQSQLQNLQSQLQNLQ*
```

### 2. BED Output (`<prefix>.gencode.bed`)

Format (1-based coordinates):
```
chr	start	end	ORF_NAME	STUDY_ID	strand
```

Example:
```
chr1	100000	100333	DDX11L1_100000_111aa	TEST_RIBOTRICER	+
chr1	200000	200150	WASH7P_200000_50aa	TEST_RIBOTRICER	+
```

## Filtering Logic

The script applies two quality filters:

### 1. Length Filter (`--min_length`)

Only ORFs with `length_aa >= min_length` are included.

**Default**: 16 amino acids

**Calculation**:
```python
length_aa = length_nt / 3
```

### 2. Phase Score Filter (`--min_phase_score`)

Only ORFs with `phase_score >= min_phase_score` are included.

**Default**: 0.5

**What is phase score?**
- Metric for ribosome 3-nt periodicity (quality of translation signal)
- Range: 0 (poor) to 1 (perfect)
- Higher score = stronger translation signal

**Recommended values**:
- `0.5`: Permissive (includes moderate-quality ORFs)
- `0.7`: Moderate (balanced)
- `0.9`: Stringent (high-confidence only)

## Example Test Run

```bash
# Navigate to test directory
cd test_data/ribotricer_to_gencode

# Run the script
python3 ../../bin/ribotricer_to_gencode.py \
    --tsv ribotricer_translating_ORFs.tsv \
    --fasta test_genome.fa \
    --study_id TEST_RIBOTRICER \
    --output_prefix ribotricer_output \
    --min_length 16 \
    --min_phase_score 0.5

# Expected output:
# Parsing Ribotricer TSV file: ribotricer_translating_ORFs.tsv
# Found 8 ORFs (>= 16 aa, phase_score >= 0.5)
# Extracting sequences from genome: test_genome.fa
# Writing gencode-riboseqORFs format files
# ✅ Successfully converted 8 ORFs to gencode format
#    Output: ribotricer_output.gencode.fa
#            ribotricer_output.gencode.bed
```

## Verify Output

```bash
# Check FASTA format
head -4 ribotricer_output.gencode.fa
# >DDX11L1_100000_111aa--TEST_RIBOTRICER
# MMMMMM...*

# Check BED format
head -2 ribotricer_output.gencode.bed
# chr1	100000	100333	DDX11L1_100000_111aa	TEST_RIBOTRICER	+

# Count ORFs
grep -c ">" ribotricer_output.gencode.fa
# 8
```

## Coordinate System

**Important**: Ribotricer uses genomic coordinates.

| Field | Coordinate System | Description |
|-------|-------------------|-------------|
| `start_codon` | 1-based genomic | Position of start codon |
| `length` | Nucleotides | ORF length in nt |
| Output BED | 1-based (gencode) | Required by gencode-riboseqORFs |

**Coordinate Calculation**:

For **positive strand** (`+`):
```python
genomic_start = start_codon
genomic_end = start_codon + length
```

For **negative strand** (`-`):
```python
genomic_end = start_codon
genomic_start = start_codon - length
```

## Test Data Description

### ribotricer_translating_ORFs.tsv

Contains 10 test ORFs:

| ORF | Gene | Type | Strand | Length (aa) | Phase Score | Status |
|-----|------|------|--------|-------------|-------------|--------|
| 1 | DDX11L1 | annotated | + | 111 | 0.95 | translating |
| 2 | WASH7P | uORF | + | 50 | 0.88 | translating |
| 3 | MIR6859-1 | dORF | - | 85 | 0.76 | translating |
| 4 | FAM138A | overlap_uORF | + | 200 | 0.92 | translating |
| 5 | LINC00115 | novel | - | 60 | 0.71 | translating |
| 6 | MIR1302-2HG | super_uORF | + | 30 | 0.65 | translating |
| 7 | TSPAN6 | annotated | - | 150 | 0.97 | translating |
| 8 | KDM5D | uORF | + | 20 | 0.55 | translating |
| 9 | DDX11L1 | overlap_dORF | + | 40 | 0.35 | non-translating |
| 10 | WASH7P | uORF | + | 15 | 0.42 | non-translating |

**With default filters** (`min_length=16`, `min_phase_score=0.5`):
- **8 ORFs pass** (ORFs 1-8)
- **2 ORFs filtered out** (ORFs 9-10: low phase score and short length)

## Troubleshooting

### Issue: "No ORFs passed the filters"

**Possible causes**:
1. All ORFs below `--min_length`
2. All ORFs below `--min_phase_score`
3. Only `non-translating` ORFs in file

**Solution**:
```bash
# Check ORF lengths
awk -F'\t' 'NR>1 {print $6/3}' ribotricer_translating_ORFs.tsv | sort -n

# Check phase scores
awk -F'\t' 'NR>1 {print $4}' ribotricer_translating_ORFs.tsv | sort -n

# Relax filters
python3 bin/ribotricer_to_gencode.py ... --min_length 10 --min_phase_score 0.3
```

### Issue: "Required column 'XXX' not found in TSV"

**Cause**: Input file is not a valid Ribotricer output

**Solution**: Verify file format
```bash
# Check header
head -1 ribotricer_translating_ORFs.tsv

# Expected columns (18 total):
# ORF_ID, ORF_type, status, phase_score, read_count, length, ...
```

### Issue: "Error parsing line X"

**Cause**: Malformed TSV row (wrong number of columns, invalid values)

**Solution**: Check specific line
```bash
# View problematic line
sed -n 'Xp' ribotricer_translating_ORFs.tsv
```

### Issue: "pyfaidx not available, using placeholder sequences"

**Cause**: pyfaidx not installed

**Impact**: Placeholder polyMethionine sequences used instead of real sequences

**Solution**:
```bash
pip install pyfaidx biopython
```

## Quality Filtering Recommendations

### Conservative (High Confidence)

```bash
--min_length 20 --min_phase_score 0.8
```

Best for:
- Publication-quality datasets
- Functional validation studies
- Low false positive rate needed

### Balanced (Default)

```bash
--min_length 16 --min_phase_score 0.5
```

Best for:
- General analysis
- Exploratory studies
- Balance between sensitivity and specificity

### Permissive (High Sensitivity)

```bash
--min_length 10 --min_phase_score 0.3
```

Best for:
- Discovery of novel ORFs
- Pooling across multiple datasets
- Prioritizing sensitivity over specificity

## Integration with gencode-riboseqORFs

After generating output files:

```bash
# Download Ensembl annotation
bash scripts/retrieve_ensembl_data.sh 110 GRCh38

# Run ORF mapper
python3 ORF_mapper_to_GENCODE_v1.1.py \
    -d Ens110/ \
    -f ribotricer_output.gencode.fa \
    -b ribotricer_output.gencode.bed \
    -o ribotricer_annotated \
    -l 16 \
    -c 0.9 \
    -m longest_string
```

## Performance

- **Runtime**: ~1-5 seconds for typical datasets (<1000 ORFs)
- **Memory**: <500 MB
- **Bottleneck**: Sequence extraction (use pyfaidx for efficiency)

## Comparison with Ribo-TISH Converter

| Feature | ribotish_to_gencode.py | ribotricer_to_gencode.py |
|---------|------------------------|--------------------------|
| Input format | GenomePos string | TSV with coordinates |
| Quality metric | P-value | Phase score |
| ORF types | TisGroup (6 types) | ORF_type (8 types) |
| Filtering | Length only | Length + phase score |
| Coordinate parsing | Regex on string | Direct column access |

## See Also

- [Ribotricer GitHub](https://github.com/smithlabcode/ribotricer)
- [Ribotricer Publication](https://doi.org/10.1093/nar/gkz878)
- [GENCODE Integration Plan](../../docs/GENCODE_INTEGRATION_IMPLEMENTATION_PLAN.md)
- [gencode-riboseqORFs](https://github.com/jorruior/gencode-riboseqORFs)

---

**Last Updated**: 2026-01-17
**Tested Python Versions**: 3.9, 3.10, 3.11
**Dependencies**: Python 3.9+, Biopython 1.79+, pyfaidx 0.7+ (optional)

**Sources:**
- [nf-core/riboseq documentation](https://nf-co.re)
- [Ribotricer GitHub repository](https://github.com/smithlabcode/ribotricer)
