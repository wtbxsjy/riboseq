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

        # Find column indices (support different Ribo-TISH versions)
        try:
            genomepos_idx = header.index('GenomePos')
            tid_idx = header.index('Tid')

            # Support both TisType and TISType
            if 'TisType' in header:
                tistype_idx = header.index('TisType')
            elif 'TISType' in header:
                tistype_idx = header.index('TISType')
            else:
                raise ValueError("'TisType' or 'TISType'")

            # Support both TisGroup and TISGroup
            if 'TisGroup' in header:
                tisgroup_idx = header.index('TisGroup')
            elif 'TISGroup' in header:
                tisgroup_idx = header.index('TISGroup')
            else:
                raise ValueError("'TisGroup' or 'TISGroup'")

            # Support both TisLen and AALen
            if 'TisLen' in header:
                tislen_idx = header.index('TisLen')
                length_in_nt = True  # TisLen is in nucleotides
            elif 'AALen' in header:
                tislen_idx = header.index('AALen')
                length_in_nt = False  # AALen is already in amino acids
            else:
                raise ValueError("'TisLen' or 'AALen'")

        except ValueError as e:
            print(f"Error: Required column not found in predict file: {e}", file=sys.stderr)
            print(f"Available columns: {', '.join(header)}", file=sys.stderr)
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
            tis_len_raw = parts[tislen_idx]

            # Handle 'None' or empty values
            if tis_len_raw in ['None', '', 'NA', 'N/A']:
                continue

            try:
                tis_len = int(tis_len_raw)
            except ValueError:
                print(f"Warning: Invalid length value '{tis_len_raw}' at {genome_pos}", file=sys.stderr)
                continue

            # Parse genome position: chr:start-end:strand
            match = re.match(r'(.+):(\d+)-(\d+):([+-])', genome_pos)
            if not match:
                print(f"Warning: Could not parse genome position: {genome_pos}", file=sys.stderr)
                continue

            chrom, start, end, strand = match.groups()
            start = int(start)
            end = int(end)

            # Calculate ORF length in amino acids
            if length_in_nt:
                # TisLen is in nucleotides, convert to amino acids
                orf_length_aa = tis_len // 3
            else:
                # AALen is already in amino acids
                orf_length_aa = tis_len

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

        # Track extraction statistics
        adjusted_count = 0
        extraction_details = []

        for orf in orfs:
            # Calculate expected nucleotide length (ORF length * 3 + 3 for stop codon)
            expected_nt_length = orf['length_aa'] * 3 + 3
            
            # Extract sequence from genome
            # Use expected length instead of genomic coordinates to avoid over-extraction
            if orf['strand'] == '+':
                # For positive strand, extract from start
                seq_end = orf['start'] - 1 + expected_nt_length
                seq = genome[orf['chrom']][orf['start']-1:seq_end]
            else:
                # For negative strand, extract from end backwards
                seq_start = orf['end'] - expected_nt_length
                seq = genome[orf['chrom']][seq_start:orf['end']]
                seq = seq.reverse.complement

            seq_str = str(seq)
            extracted_length = len(seq_str)
            
            # Verify extracted length is reasonable
            if extracted_length < expected_nt_length - 3:
                print(f"Warning: Extracted sequence too short at {orf['genome_pos']}: "
                      f"expected ~{expected_nt_length}nt, got {extracted_length}nt", file=sys.stderr)
            elif extracted_length > expected_nt_length + 10:
                print(f"Warning: Extracted sequence too long at {orf['genome_pos']}: "
                      f"expected ~{expected_nt_length}nt, got {extracted_length}nt", file=sys.stderr)

            # Trim to multiple of 3
            remainder = len(seq_str) % 3
            if remainder != 0:
                seq_str = seq_str[:-remainder]

            seq_obj = Seq(seq_str)

            try:
                # Translate and look for the first stop codon
                full_protein = str(seq_obj.translate())
                
                # Find first stop codon position
                first_stop = full_protein.find('*')
                
                if first_stop == -1:
                    # No stop codon found, add one
                    protein = full_protein + '*'
                    print(f"Warning: No stop codon in sequence at {orf['genome_pos']}, added one", file=sys.stderr)
                else:
                    # Use sequence up to and including first stop codon
                    protein = full_protein[:first_stop+1]
                
                # Calculate actual length (excluding stop codon)
                actual_length = len(protein.rstrip('*'))
                
                # Check if length matches expected
                if actual_length != orf['length_aa']:
                    extraction_details.append(
                        f"{orf['genome_pos']}: declared {orf['length_aa']}aa, "
                        f"extracted {extracted_length}nt, translated to {actual_length}aa"
                    )
                    adjusted_count += 1
                    orf['length_aa'] = actual_length
                
                orf['sequence'] = protein
                    
            except Exception as e:
                print(f"Warning: Could not translate ORF at {orf['genome_pos']}: {e}", file=sys.stderr)
                orf['sequence'] = 'M' * orf['length_aa'] + '*'

        # Write detailed extraction log to file
        if extraction_details:
            with open('sequence_extraction.log', 'w') as log:
                log.write("# ORF sequence extraction adjustments\n")
                log.write(f"# Total adjusted: {adjusted_count} out of {len(orfs)} ORFs\n")
                log.write(f"# These ORFs had length differences between Ribo-TISH report and actual translation\n\n")
                for detail in extraction_details:
                    log.write(detail + '\n')
            print(f"Adjusted {adjusted_count}/{len(orfs)} ORF lengths based on actual sequences (see sequence_extraction.log)")

        return orfs

    except ImportError:
        # Fallback: use bedtools getfasta
        print("Warning: pyfaidx not available, using bedtools getfasta", file=sys.stderr)

        import tempfile
        import os

        # Track extraction statistics
        adjusted_count = 0
        extraction_details = []

        # Create temporary BED file for bedtools
        with tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as tmp_bed:
            bed_file = tmp_bed.name
            for i, orf in enumerate(orfs):
                # Calculate expected nucleotide length
                expected_nt_length = orf['length_aa'] * 3 + 3
                
                # BED format: 0-based, half-open
                # Extract based on expected length, not full genomic coordinates
                if orf['strand'] == '+':
                    bed_start = orf['start'] - 1
                    bed_end = bed_start + expected_nt_length
                else:
                    bed_end = orf['end']
                    bed_start = bed_end - expected_nt_length
                
                tmp_bed.write(f"{orf['chrom']}\t{bed_start}\t{bed_end}\torf_{i}\t0\t{orf['strand']}\n")

        # Run bedtools getfasta
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as tmp_fa:
            fa_file = tmp_fa.name

        try:
            cmd = ['bedtools', 'getfasta', '-fi', fasta_file, '-bed', bed_file, '-fo', fa_file, '-s', '-name']
            subprocess.run(cmd, check=True, capture_output=True)

            # Parse extracted sequences
            seq_dict = {}
            for record in SeqIO.parse(fa_file, 'fasta'):
                orf_id = record.id.split('(')[0]  # Remove strand info
                seq_dict[orf_id] = str(record.seq)

            # Translate sequences
            for i, orf in enumerate(orfs):
                orf_id = f"orf_{i}"
                if orf_id in seq_dict:
                    seq_str = seq_dict[orf_id]
                    extracted_length = len(seq_str)

                    # Trim to multiple of 3
                    remainder = len(seq_str) % 3
                    if remainder != 0:
                        seq_str = seq_str[:-remainder]

                    seq_obj = Seq(seq_str)
                    try:
                        # Translate and find first stop codon
                        full_protein = str(seq_obj.translate())
                        
                        first_stop = full_protein.find('*')
                        if first_stop == -1:
                            protein = full_protein + '*'
                            print(f"Warning: No stop codon at {orf['genome_pos']}, added one", file=sys.stderr)
                        else:
                            protein = full_protein[:first_stop+1]
                        
                        # Calculate actual length
                        actual_length = len(protein.rstrip('*'))
                        
                        if actual_length != orf['length_aa']:
                            extraction_details.append(
                                f"{orf['genome_pos']}: declared {orf['length_aa']}aa, "
                                f"extracted {extracted_length}nt, translated to {actual_length}aa"
                            )
                            adjusted_count += 1
                            orf['length_aa'] = actual_length
                        
                        orf['sequence'] = protein
                            
                    except Exception as e:
                        print(f"Warning: Could not translate ORF at {orf['genome_pos']}: {e}", file=sys.stderr)
                        orf['sequence'] = 'M' * orf['length_aa'] + '*'
                else:
                    print(f"Warning: No sequence extracted for {orf['genome_pos']}", file=sys.stderr)
                    orf['sequence'] = 'M' * orf['length_aa'] + '*'

            # Cleanup
            os.unlink(bed_file)
            os.unlink(fa_file)

            # Write detailed extraction log to file
            if extraction_details:
                with open('sequence_extraction.log', 'w') as log:
                    log.write("# ORF sequence extraction adjustments (bedtools mode)\n")
                    log.write(f"# Total adjusted: {adjusted_count} out of {len(orfs)} ORFs\n")
                    log.write(f"# These ORFs had length differences between Ribo-TISH report and actual translation\n\n")
                    for detail in extraction_details:
                        log.write(detail + '\n')
                print(f"Adjusted {adjusted_count}/{len(orfs)} ORF lengths based on actual sequences (see sequence_extraction.log)")

        except subprocess.CalledProcessError as e:
            print(f"Error running bedtools getfasta: {e.stderr.decode()}", file=sys.stderr)
            print("Falling back to placeholder sequences", file=sys.stderr)
            for orf in orfs:
                orf['sequence'] = 'M' * orf['length_aa'] + '*'
        except Exception as e:
            print(f"Error during sequence extraction: {e}", file=sys.stderr)
            print("Falling back to placeholder sequences", file=sys.stderr)
            for orf in orfs:
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
            # Create ORF name with full details (gene + coordinates + length)
            # Format: GENE_START_LENGTHaa
            gene_name = orf['tid'].split('.')[0]  # Remove version
            orf_name = f"{gene_name}_{orf['start']}_{orf['length_aa']}aa"

            # Write FASTA
            fa.write(f">{orf_name}--{study_id}\n")
            fa.write(f"{orf['sequence']}\n")

            # Write BED in gencode-riboseqORFs format
            # Format: chr start end ORF_NAME STUDY_ID strand
            # Use the same ORF_NAME for both FASTA and BED to enable matching
            bed_start = orf['start']  # Already 1-based from Ribo-TISH
            bed_end = orf['end']      # Already 1-based from Ribo-TISH

            # Normalize chromosome name (remove 'chr' prefix if present to match reference)
            chrom_normalized = orf['chrom'].replace('chr', '')

            bed.write(f"{chrom_normalized}\t{bed_start}\t{bed_end}\t{orf_name}\t{study_id}\t{orf['strand']}\n")

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
