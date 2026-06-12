#!/usr/bin/env python3
"""
Lightweight reformatter for UNIFY's pre-computed per-ORF expression files.

Reads expression_summary.tsv and expression_rpkm_tpm.tsv produced by
unify_orf_predictions.py, optionally applies OCS confidence filtering,
and writes final expression output files.

Replaces the awk-based quantify_orf_expression.py + calc_orf_rpkm_tpm.py
with a sub-second pandas read/filter/write pass.
"""
import argparse
import sys
from pathlib import Path

import pandas as pd


def main():
    parser = argparse.ArgumentParser(
        description="Format UNIFY expression outputs with optional OCS filtering"
    )
    parser.add_argument(
        "--expression-summary", required=True,
        help="UNIFY expression_summary.tsv"
    )
    parser.add_argument(
        "--expression-rpkm-tpm", required=True,
        help="UNIFY expression_rpkm_tpm.tsv"
    )
    parser.add_argument(
        "--orf-confidence", default=None,
        help="ORF confidence TSV from ORF_QC (optional, for OCS filtering)"
    )
    parser.add_argument(
        "--min-ocs", type=float, default=0.0,
        help="Minimum OCS threshold (default: 0.0)"
    )
    parser.add_argument(
        "--max-orfs", type=int, default=0,
        help="Keep at most N ORFs by OCS rank (0 = keep all)"
    )
    parser.add_argument(
        "--output-summary", required=True,
        help="Filtered expression summary output"
    )
    parser.add_argument(
        "--output-rpkm-tpm", required=True,
        help="Filtered RPKM/TPM output"
    )
    args = parser.parse_args()

    # Load UNIFY pre-computed expression files
    summary = pd.read_csv(args.expression_summary, sep='\t')
    rpkm_tpm = pd.read_csv(args.expression_rpkm_tpm, sep='\t')

    print(f"Loaded {len(summary)} ORFs from expression summary")
    print(f"Loaded {len(rpkm_tpm)} ORFs from RPKM/TPM")

    # Optional OCS filtering
    if args.orf_confidence and Path(args.orf_confidence).exists():
        conf = pd.read_csv(args.orf_confidence, sep='\t')
        if 'orf_id' in conf.columns and 'ocs' in conf.columns:
            if args.min_ocs > 0:
                conf = conf[conf['ocs'] >= args.min_ocs]
            if args.max_orfs > 0:
                conf = conf.nlargest(args.max_orfs, 'ocs')
            keep_ids = set(conf['orf_id'])
            summary = summary[summary['orf_id'].isin(keep_ids)]
            rpkm_tpm = rpkm_tpm[rpkm_tpm['orf_id'].isin(keep_ids)]
            print(f"After OCS filter: {len(summary)} ORFs (min_ocs={args.min_ocs}, max_orfs={args.max_orfs})")
        else:
            print("Warning: OCS file missing 'orf_id' or 'ocs' column, skipping filter")
    else:
        if args.orf_confidence:
            print(f"Note: OCS file not found ({args.orf_confidence}), skipping confidence filter")

    # Write filtered outputs
    summary.to_csv(args.output_summary, sep='\t', index=False)
    rpkm_tpm.to_csv(args.output_rpkm_tpm, sep='\t', index=False)
    print(f"Wrote {len(summary)} rows to {args.output_summary}")
    print(f"Wrote {len(rpkm_tpm)} rows to {args.output_rpkm_tpm}")


if __name__ == "__main__":
    main()
