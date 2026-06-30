"""Convert Ribotricer translating_ORFs.tsv to GENCODE ORF mapper input format.

Calls ribotricer_to_gencode.py directly (no subprocess overhead).
"""

import os
import logging

from orfont.core.scripts_bridge import call_convert_ribotricer

logger = logging.getLogger(__name__)


def convert(tsv_file, fasta_file, study_id, output_dir='.',
            min_length=16, min_phase_score=0.5):
    """Convert Ribotricer TSV to GENCODE format — direct Python call.

    Args:
        tsv_file: Ribotricer *_translating_ORFs.tsv
        fasta_file: genome FASTA
        study_id: sample/study identifier
        output_dir: output directory
        min_length: minimum ORF length in AA
        min_phase_score: minimum phase score filter

    Returns:
        dict: {'fasta': path, 'bed': path}
    """
    os.makedirs(output_dir, exist_ok=True)
    prefix = os.path.join(output_dir, 'ribotricer_gencode')

    argv = ['--tsv', tsv_file,
            '--fasta', fasta_file,
            '--study_id', study_id,
            '--output_prefix', prefix,
            '--min_length', str(min_length),
            '--min_phase_score', str(min_phase_score)]

    logger.info(f"Converting Ribotricer (direct): {argv}")
    call_convert_ribotricer(argv)

    return {
        'fasta': f'{prefix}.gencode.fa',
        'bed': f'{prefix}.gencode.bed',
    }
