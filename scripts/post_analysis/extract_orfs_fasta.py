#!/usr/bin/env python3
"""Generate `<prefix>_orfs.fa` — raw spliced ORF nucleotide sequences with coordinate headers.

Mirrors the header convention of the original PRJEB26593 QC outputs:

    >{orf_id}::{chrom}:{start0}-{end}({strand})
    ATG...

Coordinates are BED-style (0-based start) — i.e. the metadata's 1-based `start`
minus 1; `end` is used as-is. Sequence comes from the unify metadata `sequence`
column (spliced, sense-strand CDS as produced by unify_orf_predictions.py).

Companion outputs (strict-CDS NA/AA, trailing stop stripped and trimmed to
3*len(AA)) come from extract_cds_fasta.py / mouse_extract_cds_fasta.py.

Usage:
  python3 extract_orfs_fasta.py --metadata M.tsv --ids pass.tsv --out prefix_orfs.fa
"""

import argparse

import pandas as pd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--metadata", required=True, help="unified_orfs.metadata.tsv")
    ap.add_argument("--ids", required=True,
                    help="pass-list TSV containing an orf_id column")
    ap.add_argument("--out", required=True, help="output FASTA path")
    a = ap.parse_args()

    ids = set(pd.read_csv(a.ids, sep="\t", usecols=["orf_id"], dtype=str)["orf_id"])
    print(f"ids: {len(ids):,}", flush=True)

    usecols = ["orf_id", "chrom", "strand", "start", "end", "sequence"]
    n = n_empty = 0
    seen = set()
    with open(a.out, "w") as fo:
        for ch in pd.read_csv(a.metadata, sep="\t", usecols=usecols,
                              dtype={"orf_id": str}, chunksize=500_000,
                              low_memory=False):
            m = ch["orf_id"].isin(ids)
            if not m.any():
                continue
            for r in ch[m].itertuples(index=False):
                s = r.sequence
                if not isinstance(s, str) or not s:
                    n_empty += 1
                    continue
                fo.write(f">{r.orf_id}::{r.chrom}:{r.start - 1}-{r.end}({r.strand})\n"
                         f"{s.upper()}\n")
                seen.add(r.orf_id)
                n += 1

    missing = len(ids) - len(seen)
    print(f"written: {n:,}  empty_sequence: {n_empty}  "
          f"ids_not_found_in_metadata: {missing}", flush=True)
    if missing:
        print("WARNING: some pass-list ids had no metadata row", flush=True)
    print(f"→ {a.out}", flush=True)


if __name__ == "__main__":
    main()
