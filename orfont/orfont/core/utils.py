"""Shared constants and utility functions — single canonical definition.

All modules (orfont and external pipeline scripts) import from here.
"""

from typing import List, Optional

# ---------------------------------------------------------------------------
# Standard codon table (DNA → amino acid)
# ---------------------------------------------------------------------------

CODON_TABLE = {
    'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
    'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
    'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
    'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
    'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
    'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*',
    'TGC': 'C', 'TGT': 'C', 'TGA': '*', 'TGG': 'W',
}


def translate_sequence(nt_seq):
    """Translate a nucleotide sequence to amino acids.

    Returns empty string for empty input. Unknown codons → 'X'. N-containing codons → 'X'.
    """
    if not nt_seq:
        return ""
    nt_seq = nt_seq.upper().replace('U', 'T')
    protein = []
    for i in range(0, len(nt_seq) - 2, 3):
        codon = nt_seq[i:i + 3]
        if 'N' in codon:
            protein.append('X')
        else:
            protein.append(CODON_TABLE.get(codon, 'X'))
    return ''.join(protein)


def chrom_aliases(chrom: Optional[str]) -> List[str]:
    """Generate common chromosome-name aliases (chrM/MT/M etc.)."""
    if not chrom:
        return []
    aliases = [chrom]
    if chrom.startswith('chr'):
        base = chrom[3:]
        aliases.append(base)
        if base == 'M':
            aliases.append('MT')
        elif base == 'MT':
            aliases.append('M')
    else:
        aliases.append(f"chr{chrom}")
        if chrom == 'MT':
            aliases.extend(['chrM', 'M'])
        elif chrom == 'M':
            aliases.extend(['chrM', 'MT'])
    seen = set()
    result = []
    for a in aliases:
        if a and a not in seen:
            seen.add(a)
            result.append(a)
    return result
