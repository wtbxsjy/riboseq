"""Sequence extraction utilities."""

import sys

try:
    from Bio.Seq import Seq
    HAS_BIOPYTHON = True
except ImportError:
    Seq = None
    HAS_BIOPYTHON = False

try:
    from pyfaidx import Fasta
    HAS_PYFAIDX = True
except ImportError:
    Fasta = None
    HAS_PYFAIDX = False

from orfont.core.utils import CODON_TABLE, translate_sequence


def extract_sequence_from_genome(fasta_path, chrom, start_1based, end_1based, strand):
    """Extract nucleotide sequence from a genome FASTA using pyfaidx."""
    if not HAS_PYFAIDX:
        return None
    try:
        genome = Fasta(fasta_path)
        # pyfaidx is 1-based, inclusive
        seq = genome[chrom][start_1based - 1:end_1based]
        nt_seq = seq.seq.upper()
        if strand == '-':
            bio_seq = Seq(nt_seq)
            nt_seq = str(bio_seq.reverse_complement())
        return nt_seq
    except Exception as e:
        print(f"Warning: sequence extraction failed: {e}", file=sys.stderr)
        return None


def extract_sequence_bedtools(fasta_path, chrom, start_1based, end_1based, strand):
    """Fallback: extract sequence using bedtools getfasta."""
    import subprocess
    import tempfile
    try:
        with tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as tmp:
            tmp.write(f"{chrom}\t{start_1based - 1}\t{end_1based}\t.\t.\t{strand}\n")
            tmp_path = tmp.name
        result = subprocess.run(
            ['bedtools', 'getfasta', '-fi', fasta_path, '-bed', tmp_path,
             '-fo', '/dev/stdout', '-s', '-name'],
            capture_output=True, text=True, check=True)
        lines = result.stdout.strip().split('\n')
        return ''.join(lines[1:]).upper() if len(lines) > 1 else ""
    except Exception as e:
        print(f"Warning: bedtools extraction failed: {e}", file=sys.stderr)
        return None
    finally:
        import os
        if os.path.exists(tmp_path):
            os.remove(tmp_path)
