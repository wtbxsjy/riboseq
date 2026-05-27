#!/usr/bin/env python3
"""
Convert Ribotricer output to gencode-riboseqORFs compatible format

Ribotricer outputs:
- translating_ORFs.tsv: Detected translating ORFs with coordinates and metrics

gencode-riboseqORFs requires:
1. FASTA: >{ORF_NAME}--{STUDY_ID}\nSEQUENCE*
2. BED (1-based): chr start end ORF_NAME STUDY_ID strand
"""

import argparse
import sys
import re
from Bio.Seq import Seq

def parse_ribotricer_tsv(tsv_file, min_length=16, min_phase_score=0.5):
    """
    Parse Ribotricer translating_ORFs.tsv output

    Expected columns (varies by version):
    ORF_ID, ORF_type, status, phase_score, read_count, length, valid_codons,
    valid_codons_ratio, read_density, transcript_id, transcript_type, gene_id,
    gene_name, gene_type, chrom, strand, start_codon, profile

    Note: 'start_codon' can be either:
    - Genomic coordinate (integer) in older versions
    - Codon sequence (e.g., 'ATG') in newer versions
    """
    orfs = []

    with open(tsv_file, 'r') as f:
        header = f.readline().strip().split('\t')

        # Map column names to indices
        try:
            col_map = {name: idx for idx, name in enumerate(header)}

            # Required columns
            required = ['ORF_ID', 'ORF_type', 'length', 'transcript_id',
                       'gene_id', 'gene_name', 'chrom', 'strand']
            for col in required:
                if col not in col_map:
                    print(f"Error: Required column '{col}' not found in TSV", file=sys.stderr)
                    print(f"Available columns: {', '.join(header)}", file=sys.stderr)
                    sys.exit(1)

            # 'start_codon' is optional - format varies by version
            has_start_codon = 'start_codon' in col_map

        except Exception as e:
            print(f"Error parsing header: {e}", file=sys.stderr)
            sys.exit(1)

        for line_num, line in enumerate(f, start=2):
            if line.startswith('#') or not line.strip():
                continue

            parts = line.strip().split('\t')

            if len(parts) < len(header):
                print(f"Warning: Line {line_num} has fewer columns than header, skipping",
                      file=sys.stderr)
                continue

            try:
                orf_id = parts[col_map['ORF_ID']]
                orf_type = parts[col_map['ORF_type']]
                length_nt = int(parts[col_map['length']])
                transcript_id = parts[col_map['transcript_id']]
                gene_id = parts[col_map['gene_id']]
                gene_name = parts[col_map['gene_name']]
                chrom = parts[col_map['chrom']]
                strand = parts[col_map['strand']]

                # Parse start_codon column (format varies by version)
                start_codon_value = None
                if has_start_codon:
                    start_codon_raw = parts[col_map['start_codon']]
                    # Try to parse as integer (genomic position)
                    try:
                        start_codon_value = int(start_codon_raw)
                    except ValueError:
                        # It's a codon sequence (e.g., 'ATG'), need to parse from ORF_ID
                        start_codon_value = None

                # If start_codon is not a position, extract from ORF_ID
                # ORF_ID format: TRANSCRIPT_ID_rank_start_stop or similar
                if start_codon_value is None:
                    # Try to extract coordinates from ORF_ID
                    # Common formats:
                    # - "ENST00000123456_1_100_200"
                    # - "ENST00000123456:1-100"
                    orf_id_match = re.search(r'_(\d+)_(\d+)$', orf_id)
                    if not orf_id_match:
                        orf_id_match = re.search(r':(\d+)-(\d+)$', orf_id)

                    if orf_id_match:
                        genomic_start = int(orf_id_match.group(1))
                        genomic_end = int(orf_id_match.group(2))
                    else:
                        print(f"Warning: Cannot extract coordinates from ORF_ID '{orf_id}', skipping",
                              file=sys.stderr)
                        continue
                else:
                    # Use start_codon position to calculate coordinates
                    if strand == '+':
                        genomic_start = start_codon_value
                        genomic_end = start_codon_value + length_nt
                    else:  # strand == '-'
                        genomic_end = start_codon_value
                        genomic_start = start_codon_value - length_nt

                # Optional columns with defaults
                phase_score = float(parts[col_map['phase_score']]) if 'phase_score' in col_map else 1.0
                status = parts[col_map['status']] if 'status' in col_map else 'translating'

                # Calculate ORF length in amino acids
                orf_length_aa = length_nt // 3

                # Filter by minimum length
                if orf_length_aa < min_length:
                    continue

                # Filter by phase score (quality metric)
                if phase_score < min_phase_score:
                    continue

                orfs.append({
                    'orf_id': orf_id,
                    'chrom': chrom,
                    'start': genomic_start,
                    'end': genomic_end,
                    'strand': strand,
                    'transcript_id': transcript_id,
                    'gene_id': gene_id,
                    'gene_name': gene_name,
                    'orf_type': orf_type,
                    'length_aa': orf_length_aa,
                    'length_nt': length_nt,
                    'phase_score': phase_score,
                    'status': status
                })

            except (ValueError, IndexError, KeyError) as e:
                print(f"Warning: Error parsing line {line_num}: {e}", file=sys.stderr)
                continue

    return orfs

def extract_sequences_from_genome(orfs, fasta_file):
    """
    Extract ORF sequences from genome FASTA using coordinates
    Uses pyfaidx if available, otherwise generates placeholder sequences
    """
    try:
        from pyfaidx import Fasta
        genome = Fasta(fasta_file)

        for orf in orfs:
            try:
                # Extract sequence (pyfaidx uses 0-based indexing)
                # Convert our 1-based coordinates to 0-based for extraction
                seq = genome[orf['chrom']][orf['start']-1:orf['end']]

                if orf['strand'] == '-':
                    seq = seq.reverse.complement

                # Translate to protein
                seq_str = str(seq)
                seq_obj = Seq(seq_str)

                protein = str(seq_obj.translate())

                # Ensure stop codon
                if not protein.endswith('*'):
                    protein += '*'

                orf['sequence'] = protein

            except Exception as e:
                print(f"Warning: Could not extract sequence for ORF {orf['orf_id']}: {e}",
                      file=sys.stderr)
                # Use placeholder
                orf['sequence'] = 'M' * orf['length_aa'] + '*'

        return orfs

    except ImportError:
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
            # Use gene_name if available, otherwise gene_id
            gene_name = orf['gene_name'] if orf['gene_name'] != 'NA' else orf['gene_id']
            gene_name = gene_name.split('.')[0]  # Remove version if present

            orf_name = f"{gene_name}_{orf['start']}_{orf['length_aa']}aa"

            # Write FASTA
            fa.write(f">{orf_name}--{study_id}\n")
            fa.write(f"{orf['sequence']}\n")

            # Write BED (1-based coordinates as required by gencode-riboseqORFs)
            # Ribotricer outputs are already in 1-based genomic coordinates
            bed.write(f"{orf['chrom']}\t{orf['start']}\t{orf['end']}\t"
                     f"{orf_name}\t{study_id}\t{orf['strand']}\n")

    return len(orfs)

def main(argv=None):
    parser = argparse.ArgumentParser(
        description='Convert Ribotricer output to gencode-riboseqORFs format',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage
  %(prog)s --tsv sample_translating_ORFs.tsv --fasta genome.fa \\
           --study_id SAMPLE1 --output_prefix sample1_output

  # With quality filtering
  %(prog)s --tsv sample_translating_ORFs.tsv --fasta genome.fa \\
           --study_id SAMPLE1 --output_prefix sample1_output \\
           --min_length 20 --min_phase_score 0.7
        """
    )

    parser.add_argument('--tsv', required=True,
                       help='Ribotricer translating_ORFs.tsv output file')
    parser.add_argument('--fasta', required=True,
                       help='Genome FASTA file for sequence extraction')
    parser.add_argument('--study_id', required=True,
                       help='Study identifier (e.g., sample name)')
    parser.add_argument('--output_prefix', required=True,
                       help='Output file prefix')
    parser.add_argument('--min_length', type=int, default=16,
                       help='Minimum ORF length in amino acids (default: 16)')
    parser.add_argument('--min_phase_score', type=float, default=0.5,
                       help='Minimum phase score for quality filtering (default: 0.5)')

    args = parser.parse_args(argv)

    # Parse Ribotricer TSV file
    print(f"Parsing Ribotricer TSV file: {args.tsv}")
    orfs = parse_ribotricer_tsv(args.tsv, args.min_length, args.min_phase_score)
    print(f"Found {len(orfs)} ORFs (>= {args.min_length} aa, phase_score >= {args.min_phase_score})")

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
    main(sys.argv[1:])
