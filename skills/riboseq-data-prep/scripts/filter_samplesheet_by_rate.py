#!/usr/bin/env python3
"""Filter a riboseq samplesheet by per-sample alignment rates (HQ filtering).

Historical workflow (maize: 100 -> 79 samples): salmon unique alignment
rates per sample were collected into a TSV, and samples with rate < 20%
were removed before running the pipeline (low-quality samples break
RIBOSEQC downstream).

Rate file format: whitespace/TSV lines, default layout is
`<rate> ... <sample>` (rate in column 1, sample in column 3).

Usage:
    python3 filter_samplesheet_by_rate.py \
        --rates /tmp/alignment_rates.txt \
        --samplesheet samplesheet.csv \
        --output samplesheet_hq.csv \
        --min-rate 20 \
        --rate-col 1 --sample-col 3
"""

import argparse
import csv
import sys


def parse_rate_file(path, rate_col, sample_col):
    """Return {sample: rate}. Columns are 1-based. Skips a header line if
    the rate column does not parse as a float."""
    rates = {}
    with open(path) as fh:
        for raw in fh:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < max(rate_col, sample_col):
                continue
            sample = parts[sample_col - 1]
            try:
                rate = float(parts[rate_col - 1])
            except ValueError:
                continue  # header line
            rates[sample] = rate
    return rates


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--rates", required=True, help="per-sample alignment rates file")
    ap.add_argument("--samplesheet", required=True, help="samplesheet.csv to filter")
    ap.add_argument("--output", required=True, help="filtered samplesheet (header preserved)")
    ap.add_argument("--min-rate", type=float, default=20.0,
                    help="keep samples with rate >= this (default: 20)")
    ap.add_argument("--rate-col", type=int, default=1, help="1-based rate column (default: 1)")
    ap.add_argument("--sample-col", type=int, default=3, help="1-based sample column (default: 3)")
    ap.add_argument("--sample-col-csv", type=int, default=1,
                    help="1-based sample column in the samplesheet (default: 1)")
    args = ap.parse_args()

    rates = parse_rate_file(args.rates, args.rate_col, args.sample_col)
    if not rates:
        sys.exit("ERROR: no sample rates parsed from %s" % args.rates)

    kept, removed = [], []
    with open(args.samplesheet, newline="") as fh:
        reader = csv.reader(fh)
        header = next(reader)
        for row in reader:
            sample = row[args.sample_col_csv - 1]
            rate = rates.get(sample)
            if rate is None:
                print("WARN: %s has no rate entry - keeping" % sample, file=sys.stderr)
                kept.append(row)
            elif rate >= args.min_rate:
                kept.append(row)
            else:
                removed.append((sample, rate))

    with open(args.output, "w", newline="") as fh:
        writer = csv.writer(fh, lineterminator="\n")
        writer.writerow(header)
        writer.writerows(kept)

    print("kept:   %d samples (rate >= %.1f)" % (len(kept), args.min_rate))
    print("removed: %d samples" % len(removed))
    for sample, rate in sorted(removed, key=lambda x: x[1]):
        print("  - %s (%.2f)" % (sample, rate))
    print("output: %s" % args.output)


if __name__ == "__main__":
    main()
