"""Convert Ribo-TISH predict output to GENCODE ORF mapper input format.

Calls ribotish_to_gencode.py directly (no subprocess overhead).
"""

import os
import logging

from orfont.core.scripts_bridge import call_convert_ribotish

logger = logging.getLogger(__name__)


def convert(predict_file, fasta_file, study_id, output_dir='.',
            quality_file=None, gtf_file=None, min_length=16):
    """Convert Ribo-TISH predict to GENCODE format — direct Python call.

    Args:
        predict_file: Ribo-TISH *_pred.txt output
        fasta_file: genome FASTA
        study_id: sample/study identifier
        output_dir: output directory
        quality_file: optional Ribo-TISH quality file
        gtf_file: optional GTF annotation
        min_length: minimum ORF length in AA

    Returns:
        dict: {'fasta': path, 'bed': path}
    """
    os.makedirs(output_dir, exist_ok=True)
    prefix = os.path.join(output_dir, 'ribotish_gencode')

    argv = ['--predict', predict_file,
            '--fasta', fasta_file,
            '--study_id', study_id,
            '--output_prefix', prefix,
            '--min_length', str(min_length)]

    if quality_file:
        argv.extend(['--quality', quality_file])
    if gtf_file:
        argv.extend(['--gtf', gtf_file])

    logger.info(f"Converting Ribo-TISH (direct): {argv}")
    call_convert_ribotish(argv)

    return {
        'fasta': f'{prefix}.gencode.fa',
        'bed': f'{prefix}.gencode.bed',
    }
