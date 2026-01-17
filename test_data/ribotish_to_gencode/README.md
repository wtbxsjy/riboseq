# ribotish_to_gencode.py - Test Documentation

## Overview

The `ribotish_to_gencode.py` script converts Ribo-TISH ORF predictions to gencode-riboseqORFs compatible format.

## Script Location

```bash
bin/ribotish_to_gencode.py
```

## Usage

```bash
python3 bin/ribotish_to_gencode.py \
    --predict <ribotish_predict.txt> \
    --fasta <genome.fa> \
    --study_id <STUDY_ID> \
    --output_prefix <output_prefix> \
    [--min_length 16] \
    [--gtf annotation.gtf] \
    [--quality ribotish_quality.txt]
```

## Parameters

| Parameter | Required | Default | Description |
|-----------|----------|---------|-------------|
| `--predict` | Yes | - | Ribo-TISH predict output file |
| `--fasta` | Yes | - | Genome FASTA file for sequence extraction |
| `--study_id` | Yes | - | Study identifier (e.g., sample name, dataset ID) |
| `--output_prefix` | Yes | - | Prefix for output files |
| `--min_length` | No | 16 | Minimum ORF length in amino acids |
| `--gtf` | No | - | GTF annotation file (optional, for future use) |
| `--quality` | No | - | Ribo-TISH quality file (optional, for future use) |

## Input Formats

### Ribo-TISH Predict File

Tab-delimited file with the following required columns:

- `GenomePos`: Genomic position in format `chr:start-end:strand`
- `Tid`: Transcript ID (e.g., ENST00000456328.2)
- `TisType`: Translation initiation site type (e.g., ATG, CTG)
- `TisGroup`: ORF classification (e.g., annotated, uORF, dORF)
- `TisLen`: ORF length in nucleotides

Example:
```
GenomePos	Tid	TisType	TisGroup	TisLen	TisCount	TisPvalue	...
chr1:100000-100333:+	ENST00000456328.2	ATG	annotated	333	150	1.2e-10	...
chr1:200000-200150:+	ENST00000515242.2	ATG	uORF	150	80	2.5e-06	...
```

### Genome FASTA

Standard FASTA format:
```
>chr1
ATGGCGGCGGGCACGCTG...
>chr2
ATGAAGCCGCCTGCGGCA...
```

## Output Formats

### 1. FASTA Output (`<prefix>.gencode.fa`)

Format required by gencode-riboseqORFs:
```
>{ORF_NAME}--{STUDY_ID}
SEQUENCE*
```

Example:
```
>ENST00000456328_100000_111aa--TEST_STUDY
MAAGTLQSQLQNLQSQLQNLQSQLQNLQ*
```

**Naming Convention**:
- `ORF_NAME`: `{GENE_ID}_{GENOMIC_START}_{LENGTH}aa`
  - `GENE_ID`: Transcript ID without version (e.g., ENST00000456328)
  - `GENOMIC_START`: 1-based genomic start position
  - `LENGTH`: ORF length in amino acids
- All sequences **must** end with stop codon `*`

### 2. BED Output (`<prefix>.gencode.bed`)

**Important**: Uses **1-based coordinates** (gencode-riboseqORFs requirement)

Format:
```
chr	start	end	ORF_NAME	STUDY_ID	strand
```

Example:
```
chr1	100000	100333	ENST00000456328_100000_111aa	TEST_STUDY	+
chr1	200000	200150	ENST00000515242_200000_50aa	TEST_STUDY	+
```

## Example Test Run

```bash
# Navigate to test directory
cd test_data/ribotish_to_gencode

# Run the script
python3 ../../bin/ribotish_to_gencode.py \
    --predict ribotish_predict.txt \
    --fasta test_genome.fa \
    --study_id TEST_STUDY \
    --output_prefix test_output \
    --min_length 16

# Expected output:
# Parsing Ribo-TISH predict file: ribotish_predict.txt
# Found 8 ORFs (>= 16 aa)
# Extracting sequences from genome: test_genome.fa
# Writing gencode-riboseqORFs format files
# ✅ Successfully converted 8 ORFs to gencode format
#    Output: test_output.gencode.fa
#            test_output.gencode.bed
```

## Verify Output

```bash
# Check FASTA format
head -4 test_output.gencode.fa
# >ENST00000456328_100000_111aa--TEST_STUDY
# MMMMMM...*

# Check BED format (1-based coordinates)
head -2 test_output.gencode.bed
# chr1	100000	100333	ENST00000456328_100000_111aa	TEST_STUDY	+

# Count ORFs
grep -c ">" test_output.gencode.fa
# 8
wc -l test_output.gencode.bed
# 8
```

## Coordinate System Conversion

**Important Note**: Ribo-TISH outputs are already in **1-based coordinates**, so the script maintains them as-is:

| Format | Coordinate System | Example |
|--------|-------------------|---------|
| Ribo-TISH GenomePos | 1-based (closed interval) | chr1:100000-100333 |
| Output BED | 1-based (gencode requirement) | chr1 100000 100333 |
| Standard BED | 0-based (half-open) | chr1 99999 100333 |

⚠️ **This is different from standard BED format!** gencode-riboseqORFs requires 1-based coordinates in BED files.

## Sequence Extraction

The script attempts to extract actual ORF sequences from the genome FASTA:

1. **With pyfaidx** (recommended):
   ```bash
   pip install pyfaidx biopython
   ```
   - Extracts real sequences from genome
   - Handles reverse complement for `-` strand
   - Translates to protein sequence

2. **Without pyfaidx** (fallback):
   - Uses placeholder sequences (polyMethionine + stop)
   - Example: `MMMMMMMMMM*` for 10 aa ORF

## Test Data Description

### ribotish_predict.txt

Contains 8 test ORFs representing different categories:

| ORF | Type | Strand | Length (aa) | Description |
|-----|------|--------|-------------|-------------|
| 1 | annotated | + | 111 | Canonical annotated ORF |
| 2 | uORF | + | 50 | Upstream ORF |
| 3 | dORF | - | 85 | Downstream ORF |
| 4 | uoORF | + | 200 | Upstream overlapping ORF |
| 5 | intORF | - | 60 | Internal out-of-frame ORF |
| 6 | lncRNA_ORF | + | 30 | lncRNA-derived ORF |
| 7 | annotated | - | 150 | Annotated ORF on X chromosome |
| 8 | uORF | + | 20 | Short uORF on Y chromosome |

### test_genome.fa

Simplified genome FASTA with:
- chr1, chr2, chr3: Autosomes
- chrX, chrY: Sex chromosomes
- Both strands represented

## Troubleshooting

### Issue: "No ORFs passed the filters"

**Cause**: All ORFs are shorter than `--min_length`

**Solution**: Reduce `--min_length` or check input file format

```bash
# Check ORF lengths in predict file
awk -F'\t' 'NR>1 {print $5/3}' ribotish_predict.txt | sort -n
```

### Issue: "Could not parse genome position"

**Cause**: Invalid GenomePos format

**Expected format**: `chr:start-end:strand` (e.g., `chr1:100000-100333:+`)

**Solution**: Verify Ribo-TISH output format

### Issue: "pyfaidx not available, using placeholder sequences"

**Cause**: pyfaidx not installed

**Solution**: Install pyfaidx for real sequence extraction:
```bash
pip install pyfaidx biopython
```

### Issue: Translation errors

**Cause**: Sequence length not divisible by 3, or non-ATCG characters

**Solution**: Check genome FASTA quality and ORF coordinates

## Integration with gencode-riboseqORFs

After generating the output files, use them with gencode-riboseqORFs:

```bash
# Download Ensembl annotation
bash scripts/retrieve_ensembl_data.sh 110 GRCh38

# Run ORF mapper
python3 ORF_mapper_to_GENCODE_v1.1.py \
    -d Ens110/ \
    -f test_output.gencode.fa \
    -b test_output.gencode.bed \
    -o output_name \
    -l 16 \
    -c 0.9 \
    -m longest_string
```

## Performance

- **Runtime**: ~1-5 seconds for typical datasets (\u003c1000 ORFs)
- **Memory**: \u003c500 MB (genome FASTA loaded by pyfaidx is memory-mapped)
- **Bottleneck**: Sequence extraction (use pyfaidx for efficiency)

## Future Enhancements

- [ ] Use Ribo-TISH quality metrics for filtering
- [ ] Extract sequences from transcriptome instead of genome
- [ ] Add GTF-based validation of ORF coordinates
- [ ] Support for splice-aware ORF sequence extraction
- [ ] Parallel processing for large datasets

## See Also

- [GENCODE Integration Implementation Plan](../../docs/GENCODE_INTEGRATION_IMPLEMENTATION_PLAN.md)
- [gencode-riboseqORFs GitHub](https://github.com/jorruior/gencode-riboseqORFs)
- [Ribo-TISH Documentation](https://github.com/zhpn1024/ribotish)

---

**Last Updated**: 2026-01-17
**Tested Python Versions**: 3.9, 3.10, 3.11
**Dependencies**: Python 3.9+, Biopython 1.79+, pyfaidx 0.7+ (optional)
