"""Classification wrapper — dispatches to GENCODE, ORFquant, or ORF-type.

Calls classification scripts directly (no subprocess for Python scripts).
"""

import os
import sys
import logging

from orfont.core.scripts_bridge import (
    call_classify_gencode,
    call_classify_orfquant,
    call_classify_orftype,
)

logger = logging.getLogger(__name__)


def classify_gencode(bed_path, metadata_path, ensembl_dir, output_dir='.',
                     gencode_impl='original', cpus=1, extra_args=None):
    """Run GENCODE/Ensembl ORF classification — direct Python call.

    Args:
        bed_path: unified BED12 file
        metadata_path: unified metadata TSV
        ensembl_dir: Ensembl annotation directory (5 standard files)
        output_dir: output directory
        gencode_impl: 'original', 'fast', or 'indexed_fast'
        cpus: number of CPUs
        extra_args: additional CLI args
    """
    os.makedirs(output_dir, exist_ok=True)

    input_base = os.path.splitext(bed_path)[0]
    argv = ['--mode', 'gencode',
            '--input', input_base,
            '--output_dir', output_dir,
            '--ensembl_dir', ensembl_dir,
            '--gencode_impl', gencode_impl,
            '--cpus', str(cpus)]

    if extra_args:
        argv.extend(extra_args.split())

    logger.info(f"Running GENCODE classification (direct): {argv}")
    call_classify_gencode(argv)
    return os.path.join(output_dir, 'gencode_results.orfs.out')


def classify_orfquant(gtf_path, metadata_path, ref_gtf, output_dir='.',
                      cpus=1, parallel=True, extra_args=None):
    """Run ORFquant-based classification with mirai parallelization.

    Args:
        gtf_path: unified ORF GTF file
        metadata_path: unified metadata TSV
        ref_gtf: reference annotation GTF
        output_dir: output directory
        cpus: number of CPU cores for mirai daemons
        parallel: enable mirai-based parallel classification
        extra_args: additional CLI args (unused for now)
    """
    os.makedirs(output_dir, exist_ok=True)

    output_path = os.path.join(output_dir, 'orfquant_classification.tsv')
    output_prefix = os.path.join(output_dir, 'orfquant_results')

    logger.info(f"Running ORFquant classification (Rscript, parallel={parallel}, cpus={cpus})")
    call_classify_orfquant(
        gtf_path=gtf_path,
        annotation_path=ref_gtf,
        output_path=output_path,
        metadata_path=metadata_path,
        output_prefix=output_prefix,
        cpus=cpus,
        parallel=parallel)

    return output_path


def classify_orftype(metadata_path, ref_gtf, output_dir='.',
                     cpus=1, extra_args=None):
    """Run ORF-type classification (gene-level) — direct Python call."""
    os.makedirs(output_dir, exist_ok=True)

    argv = ['--input', metadata_path,
            '--gtf', ref_gtf,
            '--output', os.path.join(output_dir, 'orftype_classification.tsv'),
            '--output_prefix', os.path.join(output_dir, 'orftype_results'),
            '--cpus', str(cpus)]

    if extra_args:
        argv.extend(extra_args.split())

    logger.info(f"Running ORF-type classification (direct): {argv}")
    call_classify_orftype(argv)
    return os.path.join(output_dir, 'orftype_classification.tsv')
