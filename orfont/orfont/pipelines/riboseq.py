"""Ribo-seq pipeline adapter for ORF unification and classification.

Produces the standardised output expected by nf-core/riboseq:
  - Unified BED12, GTF, metadata TSV
  - GENCODE classification results
  - ORFquant classification TSV
  - ORF-type classification TSV
"""

import os
import logging

from orfont.unification.builder import unify
from orfont.classification.wrapper import classify_gencode, classify_orfquant, classify_orftype

logger = logging.getLogger(__name__)


def run(ribotish_files=None, ribotricer_files=None, ribocode_files=None,
        orfquant_files=None, gtf_path=None, fasta_path=None,
        ensembl_dir=None, output_dir='./orf_results',
        unify_prefix='unified_orfs',
        run_gencode=True, run_orfquant=True, run_orftype=True,
        frame_merge=True, frame_merge_min_overlap=0.9,
        seq_cluster=False, bedgraph_dir=None, sample_list=None,
        gencode_impl='original', cpus=1, extra_args=None):
    """Run the complete riboseq ORF pipeline: unify + 3-way classify.

    Args:
        ribotish_files: list of Ribo-TISH *_pred.txt files
        ribotricer_files: list of Ribotricer *_translating_ORFs.tsv files
        ribocode_files: list of RiboCode *_collapsed.gtf files
        orfquant_files: list of ORFquant *_Detected_ORFs.gtf files
        gtf_path: reference GTF annotation
        fasta_path: genome FASTA
        ensembl_dir: Ensembl annotation directory (required for GENCODE mode)
        output_dir: root output directory
        unify_prefix: prefix for unified output files
        run_gencode: enable GENCODE classification
        run_orfquant: enable ORFquant classification
        run_orftype: enable ORF-type classification
        frame_merge: enable frame-aware merge
        seq_cluster: enable sequence clustering
        gencode_impl: GENCODE implementation variant
        cpus: number of CPUs
        extra_args: additional CLI args

    Returns:
        dict: paths to all output files
    """
    os.makedirs(output_dir, exist_ok=True)
    results = {}

    # Stage 1: Unification
    logger.info("Stage 1/2: ORF Unification")
    unified = unify(
        ribotish_files=ribotish_files,
        ribotricer_files=ribotricer_files,
        ribocode_files=ribocode_files,
        orfquant_files=orfquant_files,
        gtf_path=gtf_path,
        fasta_path=fasta_path,
        output_dir=output_dir,
        prefix=unify_prefix,
        frame_merge=frame_merge,
        frame_merge_min_overlap=frame_merge_min_overlap,
        seq_cluster=seq_cluster,
        bedgraph_dir=bedgraph_dir,
        sample_list=sample_list,
        extra_args=extra_args,
    )
    results['unify'] = unified

    # Stage 2: Classification (3 parallel modes)
    logger.info("Stage 2/2: ORF Classification")
    classify_outdir = os.path.join(output_dir, 'classification')

    if run_gencode and ensembl_dir:
        logger.info("  GENCODE classification ...")
        results['gencode'] = classify_gencode(
            bed_path=unified['bed'],
            metadata_path=unified['metadata'],
            ensembl_dir=ensembl_dir,
            output_dir=classify_outdir,
            gencode_impl=gencode_impl,
            cpus=cpus,
            extra_args=extra_args,
        )

    if run_orfquant:
        logger.info("  ORFquant classification ...")
        results['orfquant'] = classify_orfquant(
            gtf_path=unified['gtf'],
            metadata_path=unified['metadata'],
            ref_gtf=gtf_path,
            output_dir=classify_outdir,
            extra_args=extra_args,
        )

    if run_orftype:
        logger.info("  ORF-type classification ...")
        results['orftype'] = classify_orftype(
            metadata_path=unified['metadata'],
            ref_gtf=gtf_path,
            output_dir=classify_outdir,
            cpus=cpus,
            extra_args=extra_args,
        )

    logger.info(f"All results in: {output_dir}")
    return results
