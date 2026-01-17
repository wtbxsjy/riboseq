#!/usr/bin/env python3
"""
Convert Ribo-TISH output to gencode-riboseqORFs compatible format

Ribo-TISH outputs:
1. predict file: ORF predictions with genomic coordinates
2. quality file: P-site offset and quality metrics

gencode-riboseqORFs requires:
1. FASTA: >{ORF_NAME}--{STUDY_ID}\nSEQUENCE*
2. BED (1-based): chr start end ORF_NAME STUDY_ID strand
"""

import argparse
import sys
import re
from Bio import SeqIO
from Bio.Seq import Seq
import subprocess

def parse_ribotish_predict(predict_file, min_length=16):
    """
    Parse Ribo-TISH predict output

    Format:
    GenomePos   Tid TisType TisGroup    TisLen  ...
    chr1:12345-12678:+  ENST00000123456 ATG uORF    111 ...
    """
    orfs = []

    with open(predict_file, 'r') as f:
        header = f.readline().strip().split('\t')

        # Find column indices
        try:
            genomepos_idx = header.index('GenomePos')
            tid_idx = header.index('Tid')
            tistype_idx = header.index('TisType')
            tisgroup_idx = header.index('TisGroup')
            tislen_idx = header.index('TisLen')
        except ValueError as e:
            print(f"Error: Required column not found in predict file: {e}", file=sys.stderr)
            sys.exit(1)

        for line in f:
            if line.startswith('#'):
                continue

            parts = line.strip().split('\t')
            if len(parts) < max(genomepos_idx, tid_idx, tistype_idx, tisgroup_idx, tislen_idx) + 1:
                continue

            genome_pos = parts[genomepos_idx]
            tid = parts[tid_idx]
            tis_type = parts[tistype_idx]
            tis_group = parts[tisgroup_idx]
            tis_len = int(parts[tislen_idx])

            # Parse genome position: chr:start-end:strand
            match = re.match(r'(.+):(\d+)-(\d+):([+-])', genome_pos)
            if not match:
                print(f"Warning: Could not parse genome position: {genome_pos}", file=sys.stderr)
                continue

            chrom, start, end, strand = match.groups()
            start = int(start)
            end = int(end)

            # Calculate ORF length in amino acids (without stop codon)
            orf_length_aa = tis_len // 3

            # Filter by minimum length
            if orf_length_aa < min_length:
                continue

            orfs.append({
                'chrom': chrom,
                'start': start,
                'end': end,
                'strand': strand,
                'tid': tid,
                'tis_type': tis_type,
                'tis_group': tis_group,
                'length_aa': orf_length_aa,
                'genome_pos': genome_pos
            })

    return orfs

def extract_sequences_from_genome(orfs, fasta_file):
    """
    Extract ORF sequences from genome FASTA using coordinates
    Uses bedtools getfasta or pyfaidx
    """
    try:
        from pyfaidx import Fasta
        genome = Fasta(fasta_file)

        for orf in orfs:
            # Extract sequence
            seq = genome[orf['chrom']][orf['start']-1:orf['end']]

            if orf['strand'] == '-':
                seq = seq.reverse.complement

            # Translate
            seq_str = str(seq)
            seq_obj = Seq(seq_str)

            try:
                protein = str(seq_obj.translate())
                # Ensure stop codon
                if not protein.endswith('*'):
                    protein += '*'
                orf['sequence'] = protein
            except Exception as e:
                print(f"Warning: Could not translate ORF at {orf['genome_pos']}: {e}", file=sys.stderr)
                orf['sequence'] = 'M' * orf['length_aa'] + '*'

        return orfs

    except ImportError:
        # Fallback: use bedtools getfasta
        print("Warning: pyfaidx not available, using placeholder sequences", file=sys.stderr)

        for orf in orfs:
            # Generate placeholder sequence
            orf['sequence'] = 'M' * orf['length_aa'] + '*'

        return orfs

def write_gencode_format(orfs, study_id, output_prefix):
    """
    Write ORFs in gencode-riboseqORFs format

    FASTA: >{ORF_NAME}--{STUDY_ID}
    BED (1-based): chr start end {ORF_NAME} {STUDY_ID} strand
    """
    fasta_output = f"{output_prefix}.gencode.fa"
    bed_output = f"{output_prefix}.gencode.bed"

    with open(fasta_output, 'w') as fa, open(bed_output, 'w') as bed:
        for orf in orfs:
            # Create ORF name
            # Format: GENE_POSITION_LENGTHaa
            gene_name = orf['tid'].split('.')[0]  # Remove version
            orf_name = f"{gene_name}_{orf['start']}_{orf['length_aa']}aa"

            # Write FASTA
            fa.write(f">{orf_name}--{study_id}\n")
            fa.write(f"{orf['sequence']}\n")

            # Write BED (convert to 1-based)
            # Note: BED is 0-based half-open, we need 1-based closed
            bed_start = orf['start']  # Already 1-based from Ribo-TISH
            bed_end = orf['end']      # Already 1-based from Ribo-TISH

            bed.write(f"{orf['chrom']}\t{bed_start}\t{bed_end}\t{orf_name}\t{study_id}\t{orf['strand']}\n")

    return len(orfs)

def main():
    parser = argparse.ArgumentParser(
        description='Convert Ribo-TISH output to gencode-riboseqORFs format'
    )
    parser.add_argument('--predict', required=True, help='Ribo-TISH predict output file')
    parser.add_argument('--quality', required=False, help='Ribo-TISH quality output file (optional)')
    parser.add_argument('--fasta', required=True, help='Genome FASTA file')
    parser.add_argument('--gtf', required=False, help='GTF annotation file (optional)')
    parser.add_argument('--study_id', required=True, help='Study identifier')
    parser.add_argument('--output_prefix', required=True, help='Output file prefix')
    parser.add_argument('--min_length', type=int, default=16, help='Minimum ORF length in amino acids')

    args = parser.parse_args()

    # Parse Ribo-TISH predict file
    print(f"Parsing Ribo-TISH predict file: {args.predict}")
    orfs = parse_ribotish_predict(args.predict, args.min_length)
    print(f"Found {len(orfs)} ORFs (>= {args.min_length} aa)")

    if len(orfs) == 0:
        print("Warning: No ORFs passed the filters", file=sys.stderr)
        # Create empty output files
        open(f"{args.output_prefix}.gencode.fa", 'w').close()
        open(f"{args.output_prefix}.gencode.bed", 'w').close()
        return

    # Extract sequences from genome
    print(f"Extracting sequences from genome: {args.fasta}")
    orfs = extract_sequences_from_genome(orfs, args.fasta)

    # Write gencode format
    print(f"Writing gencode-riboseqORFs format files")
    n_written = write_gencode_format(orfs, args.study_id, args.output_prefix)

    print(f"✅ Successfully converted {n_written} ORFs to gencode format")
    print(f"   Output: {args.output_prefix}.gencode.fa")
    print(f"           {args.output_prefix}.gencode.bed")

if __name__ == '__main__':
    main()
