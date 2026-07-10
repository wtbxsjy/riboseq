#!/usr/bin/env python3
"""
Deduce the optimal PRICE read-length filter range from upstream QC outputs.

Primary source: riboWaltz *_psite_offset.tsv
  - The 'length' column lists read lengths with detectable 3-nt periodicity.
  - Take the min:max across all input files as the filter range.
  - riboWaltz runs before PRICE in the pipeline, so data is available.

Fallback source: RiboseQC *_P_sites_calcs
  - Read lengths where all==TRUE, max_inframe==TRUE, and max_coverage==TRUE.

Output format (stdout): "min:max" (e.g. "28:30") or "" (no filter).
Exit code 0 always (even with empty input — absence of data is not an error).

Usage:
    deduce_price_read_filter.py ribowaltz/*_psite_offset.tsv
    deduce_price_read_filter.py riboseqc/*_P_sites_calcs --source riboseqc
"""

import argparse
import csv
import sys
from pathlib import Path


def from_ribowaltz(files: list[str]) -> str:
    """Extract read-length range from riboWaltz psite_offset TSV files."""
    lengths: set[int] = set()
    for f in files:
        try:
            with open(f) as fh:
                reader = csv.DictReader(fh, delimiter='\t')
                for row in reader:
                    try:
                        lengths.add(int(row['length']))
                    except (KeyError, ValueError):
                        continue
        except (OSError, IOError):
            continue
    if not lengths:
        return ''
    return f'{min(lengths)}:{max(lengths)}'


def from_riboseqc(files: list[str]) -> str:
    """Extract read-length range from RiboseQC P_sites_calcs files."""
    lengths: set[int] = set()
    for f in files:
        try:
            with open(f) as fh:
                reader = csv.DictReader(fh, delimiter='\t')
                for row in reader:
                    # Only include lengths that passed all QC checks
                    if (row.get('all', '').upper() == 'TRUE'
                            and row.get('max_inframe', '').upper() == 'TRUE'
                            and row.get('max_coverage', '').upper() == 'TRUE'):
                        try:
                            lengths.add(int(row['read_length']))
                        except (KeyError, ValueError):
                            continue
        except (OSError, IOError):
            continue
    if not lengths:
        return ''
    return f'{min(lengths)}:{max(lengths)}'


def main():
    parser = argparse.ArgumentParser(
        description='Deduce optimal PRICE read-length filter range from QC outputs.')
    parser.add_argument(
        'files', nargs='+', help='QC output files (riboWaltz psite_offset.tsv or RiboseQC P_sites_calcs)')
    parser.add_argument(
        '--source', choices=['ribowaltz', 'riboseqc'], default='ribowaltz',
        help='Source of the input files (default: ribowaltz)')
    args = parser.parse_args()

    if args.source == 'riboseqc':
        result = from_riboseqc(args.files)
    else:
        result = from_ribowaltz(args.files)

    if result:
        print(result)
    # else: print nothing (empty filter = let PRICE use all read lengths)


if __name__ == '__main__':
    main()
