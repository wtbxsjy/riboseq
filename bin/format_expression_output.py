#!/usr/bin/env python3
"""
Lightweight reformatter for UNIFY's pre-computed per-ORF expression files.

Reads expression_summary.tsv and expression_rpkm_tpm.tsv produced by
unify_orf_predictions.py, optionally applies OCS confidence filtering,
and writes final expression output files.

Uses only Python stdlib (zero external dependencies) so it runs instantly
in any Python 3.7+ container without pip install overhead.
"""
import argparse
import csv
import sys
from pathlib import Path


def read_tsv_columns(path):
    """Read a TSV file and return (column_names, list_of_dicts)."""
    with open(path, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        rows = list(reader)
    return reader.fieldnames, rows


def write_tsv(path, fieldnames, rows):
    """Write a list of dicts as a TSV file."""
    with open(path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t',
                                extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)


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

    # Load UNIFY pre-computed expression files (stdlib csv, zero-cost)
    summary_cols, summary_rows = read_tsv_columns(args.expression_summary)
    rpkm_cols, rpkm_rows = read_tsv_columns(args.expression_rpkm_tpm)

    print(f"Loaded {len(summary_rows)} ORFs from expression summary")
    print(f"Loaded {len(rpkm_rows)} ORFs from RPKM/TPM")

    # Optional OCS filtering
    if args.orf_confidence and Path(args.orf_confidence).exists():
        ocs_cols, ocs_rows = read_tsv_columns(args.orf_confidence)
        if 'orf_id' in ocs_cols and 'ocs' in ocs_cols:
            # Filter by min_ocs
            if args.min_ocs > 0:
                ocs_rows = [r for r in ocs_rows if float(r.get('ocs', 0)) >= args.min_ocs]
            # Sort by OCS descending and take top N
            ocs_rows.sort(key=lambda r: float(r.get('ocs', 0)), reverse=True)
            if args.max_orfs > 0:
                ocs_rows = ocs_rows[:args.max_orfs]
            keep_ids = {r['orf_id'] for r in ocs_rows}
            summary_rows = [r for r in summary_rows if r['orf_id'] in keep_ids]
            rpkm_rows = [r for r in rpkm_rows if r['orf_id'] in keep_ids]
            print(f"After OCS filter: {len(summary_rows)} ORFs "
                  f"(min_ocs={args.min_ocs}, max_orfs={args.max_orfs})")
        else:
            print("Warning: OCS file missing 'orf_id' or 'ocs' column, skipping filter")
    elif args.orf_confidence:
        print(f"Note: OCS file not found ({args.orf_confidence}), skipping confidence filter")

    # Write filtered outputs
    write_tsv(args.output_summary, summary_cols, summary_rows)
    write_tsv(args.output_rpkm_tpm, rpkm_cols, rpkm_rows)
    print(f"Wrote {len(summary_rows)} rows to {args.output_summary}")
    print(f"Wrote {len(rpkm_rows)} rows to {args.output_rpkm_tpm}")


if __name__ == "__main__":
    main()
