#!/usr/bin/env python3
"""Extract strict-CDS NA/AA FASTA for QC-passed ORF sets from the unify metadata.

Sources:
  - unified_orfs.metadata.tsv columns `sequence` (spliced, sense-strand CDS as
    produced by unify_orf_predictions.py) and `aa_sequence` (translation).
  - Per-set pass TSVs (e.g. {prefix}_stage2_passed.tsv /
    {prefix}_filtered_len50aa.tsv) for the orf_id list and the
    `biotype_final` classification written by stage2_psite_filter.py; the
    biotype goes into the FASTA description.

Normalisation to a uniform strict-CDS convention:
  - NA  = start codon .. last coding codon (trailing TAA/TAG/TGA stripped)
  - AA  = translation without the trailing '*'
  - trailing boundary noise (1/2/4/5 extra nt after the CDS) is trimmed so that
    len(NA) == 3 * len(AA) for EVERY ORF.

⚠️ The stop-codon strip MUST be gated on `len(s) % 3 == 0`. Raw sequences carry
1-2 nt of boundary noise, so `len(s) % 3 != 0` is the norm; for a 37 nt / 12 aa
ORF the last three characters ("TAA") are not a codon at all and stripping them
truncates the CDS by 3 nt. This silently produced 23% length-inconsistent output
on PRJEB26593 (100% of it PRICE-sourced) — see CLAUDE.md gotcha 31.
`n_len_mismatch_na_vs_aa` in the stats output must be 0; treat anything else as
a bug, not a warning.

Companion output: extract_orfs_fasta.py writes the RAW spliced sequences with
coordinate headers (`>{orf_id}::{chrom}:{start0}-{end}({strand})`).

Usage:
  python3 extract_cds_fasta.py \
      --metadata <result>/orf_unification/unified_orfs.metadata.tsv \
      --out-dir  post_analysis/<project> --prefix <project> \
      --set stage2=<project>_stage2_passed.tsv \
      --set len50aa=<project>_filtered_len50aa.tsv
"""

import argparse
import json
import os

import pandas as pd

STOP_CODONS = {"TAA", "TAG", "TGA"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--metadata", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--prefix", required=True)
    ap.add_argument("--set", action="append", required=True, metavar="NAME=TSV",
                    help="pass-list TSV (needs orf_id + biotype_final columns); "
                         "repeatable. NAME becomes the output filename infix.")
    ap.add_argument("--biotype-col", default="biotype_final")
    a = ap.parse_args()

    sets = {}
    for spec in a.set:
        name, _, tsv = spec.partition("=")
        if not tsv:
            raise SystemExit(f"--set must be NAME=TSV, got {spec!r}")
        sets[name] = tsv

    bio = {}
    for name, tsv in sets.items():
        df = pd.read_csv(tsv, sep="\t", usecols=["orf_id", a.biotype_col], dtype=str)
        bio[name] = dict(zip(df["orf_id"], df[a.biotype_col]))
        print(f"{name}: {len(bio[name]):,} ids loaded", flush=True)

    union = set().union(*[set(m.keys()) for m in bio.values()])
    print(f"union: {len(union):,}", flush=True)

    # Stream the (multi-GB) metadata, keeping only the ids we need
    seq, aa = {}, {}
    usecols = ["orf_id", "sequence", "aa_sequence"]
    for ch in pd.read_csv(a.metadata, sep="\t", usecols=usecols,
                          dtype={"orf_id": str}, chunksize=2_000_000,
                          low_memory=False):
        m = ch["orf_id"].isin(union)
        if not m.any():
            continue
        sub = ch[m]
        seq.update(dict(zip(sub["orf_id"],
                            sub["sequence"].fillna("").astype(str))))
        aa.update(dict(zip(sub["orf_id"],
                           sub["aa_sequence"].fillna("").astype(str))))
    print(f"metadata rows collected: {len(seq):,}", flush=True)

    stats = {}
    for name in sets:
        na_out = f"{a.out_dir}/{a.prefix}_{name}_cds_na.fa"
        aa_out = f"{a.out_dir}/{a.prefix}_{name}_cds_aa.fa"
        n = n_stop = n_trim = n_empty = n_bad = 0
        with open(na_out, "w") as fna, open(aa_out, "w") as faa:
            for orf_id, bt in bio[name].items():
                s = seq.get(orf_id, "").upper()
                p = aa.get(orf_id, "").upper()
                # codon-aligned only — see module docstring / gotcha 31
                if len(s) % 3 == 0 and len(s) >= 3 and s[-3:] in STOP_CODONS:
                    s = s[:-3]
                    n_stop += 1
                if p.endswith("*"):
                    p = p[:-1]
                if p and len(s) > 3 * len(p):
                    s = s[:3 * len(p)]
                    n_trim += 1
                if not s or not p:
                    n_empty += 1
                    continue
                if len(s) != 3 * len(p):
                    n_bad += 1
                fna.write(f">{orf_id} biotype={bt}\n{s}\n")
                faa.write(f">{orf_id} biotype={bt}\n{p}\n")
                n += 1
        stats[name] = {
            "n_written": n, "n_stop_codon_stripped": n_stop,
            "n_trailing_nt_trimmed": n_trim, "n_empty_sequence": n_empty,
            "n_len_mismatch_na_vs_aa": n_bad,
            "na_file": os.path.basename(na_out),
            "aa_file": os.path.basename(aa_out),
        }
        print(f"{name}: written={n:,} stop_stripped={n_stop:,} "
              f"trailing_trimmed={n_trim:,} empty={n_empty} "
              f"len_mismatch={n_bad}"
              f"{'  <-- BUG, must be 0' if n_bad else ''}", flush=True)

    with open(f"{a.out_dir}/{a.prefix}_cds_fasta_stats.json", "w") as f:
        json.dump(stats, f, indent=2)
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
