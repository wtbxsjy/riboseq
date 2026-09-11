#!/usr/bin/env python3
"""Post-analysis Step 1: CDS exclusion + Stage-1 expression filter + pN distribution.

Generic version of the per-project scripts used for PRJEB26593 / GSE120762 /
rice / maize. Runs three things:

  1. CDS exclusion — drops ORFs overlapping a CDS (`is_cds_overlap` from the
     orf-type classifier) or mapped to a CDS by the GENCODE mapper
     (`orf_biotype == 'CDS'` in `gencode_results.orfs.out.gz`).
  2. Stage-1 candidate set — `{sample}_reads > 9` in >= 1 sample. pN is
     reported but NOT gated on: it is structurally >= 1 (see CLAUDE.md
     gotcha 28), so a pN threshold is vacuous.
  3. pN distribution + quantiles + plot, for the record.

Outputs (in --out-dir):
  step1_nonCDS.tsv.gz          non-CDS ORFs with length/biotype/type/psites/reads/pN
  step1_stats.json             CDS exclusion counts
  stage1_threshold_table.tsv   candidate pN thresholds -> passing ORF counts
  pn_values_summary.json       pN quantiles
  pn_distribution.png          histogram of per-sample pN (reads > 9)
  {prefix}_stage1_passed.tsv   ORFs with reads > 9 in >= 1 sample

Usage:
  python3 stage1_expression_filter.py \
      --orftype  <result>/orf_classification/orf_type/orftype_classification.tsv \
      --gencode  <result>/orf_classification/gencode/gencode_results.orfs.out.gz \
      --expression <result>/orf_unification/unified_orfs_expression_summary.tsv \
      --out-dir  post_analysis/<project> --prefix <project>
"""

import argparse
import json
import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

READS_MIN = 9          # per-sample read floor for the Stage-1 set
PN_SWEEP = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8]
QUANTILES = [0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--orftype", required=True)
    ap.add_argument("--gencode", required=True)
    ap.add_argument("--expression", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--prefix", required=True, help="project id, e.g. mouse_GSE120762")
    ap.add_argument("--title", default=None, help="plot title (defaults to prefix)")
    a = ap.parse_args()
    title = a.title or a.prefix
    os.makedirs(a.out_dir, exist_ok=True)

    # ── Step 1: CDS exclusion ─────────────────────────────────────────────
    print("Loading orftype_classification.tsv (chunked) ...", flush=True)
    usecols = ["orf_id", "length_aa", "is_cds_overlap", "orf_type_category",
               "total_psites", "total_reads", "pN"]
    chunks = pd.read_csv(a.orftype, sep="\t", usecols=usecols,
                         dtype={"orf_id": str}, chunksize=500_000)
    ort = pd.concat(list(chunks), ignore_index=True)
    print(f"  orftype rows: {len(ort):,}", flush=True)
    ort["is_cds_overlap"] = pd.to_numeric(ort["is_cds_overlap"],
                                          errors="coerce").fillna(0).astype(int)

    print("Loading gencode_results.orfs.out.gz ...", flush=True)
    genc = pd.read_csv(a.gencode, sep="\t", dtype=str,
                       usecols=["orf_biotype", "all_orf_names", "phaseI_id"])
    n_gencode = len(genc)
    print(f"  gencode rows: {n_gencode:,}", flush=True)

    # `all_orf_names` holds comma-separated unified ORF ids
    cds_ids = set()
    for v in genc.loc[genc["orf_biotype"] == "CDS", "all_orf_names"].dropna():
        for part in str(v).split(","):
            cds_ids.add(part.strip())
    print(f"  gencode CDS-mapped unified ids: {len(cds_ids):,}", flush=True)

    biotype_map = {}
    gnon = genc[genc["orf_biotype"] != "CDS"]
    for uid, bt in zip(gnon["all_orf_names"], gnon["orf_biotype"]):
        if pd.isna(uid):
            continue
        for part in str(uid).split(","):
            biotype_map.setdefault(part.strip(), bt)
    del genc, gnon

    ids = ort["orf_id"]
    excl_ov = ort["is_cds_overlap"] == 1
    excl_gen = ids.isin(cds_ids)
    keep = (~excl_ov) & (~excl_gen)
    n_overlap_excl = int(excl_ov.sum())
    n_gencode_excl = int((excl_gen & ~excl_ov).sum())
    n_keep = int(keep.sum())
    print(f"  excluded: is_cds_overlap={n_overlap_excl:,} | "
          f"gencode-CDS additional={n_gencode_excl:,} | keep={n_keep:,}", flush=True)

    step1_stats = {
        "ortype_rows": len(ort),
        "excluded_is_cds_overlap": n_overlap_excl,
        "excluded_gencode_CDS_additional": n_gencode_excl,
        "kept_nonCDS": n_keep,
        "gencode_rows_total": n_gencode,
    }

    ort = ort[keep].copy()
    ort["gencode_biotype"] = ids[keep].map(lambda x: biotype_map.get(x, "")).values
    ort.to_csv(f"{a.out_dir}/step1_nonCDS.tsv.gz", sep="\t", index=False,
               compression="gzip")
    with open(f"{a.out_dir}/step1_stats.json", "w") as f:
        json.dump(step1_stats, f, indent=2)
    print("Wrote step1_nonCDS.tsv.gz + step1_stats.json", flush=True)

    # ── Stage 1: reads/pN sweep on the non-CDS set ────────────────────────
    print("Loading expression summary ...", flush=True)
    header = pd.read_csv(a.expression, sep="\t", nrows=0).columns.tolist()
    read_cols = [c for c in header if c.endswith("_reads") and c != "total_reads"]
    pn_cols = [c for c in header if c.endswith("_pN")]
    samples = [c[: -len("_reads")] for c in read_cols]
    print(f"  samples={len(samples)}: {samples}", flush=True)

    expr = pd.read_csv(a.expression, sep="\t", usecols=["orf_id"] + read_cols + pn_cols,
                       dtype={"orf_id": str})
    keep_ids = set(ort["orf_id"])
    expr = expr[expr["orf_id"].isin(keep_ids)]
    print(f"  expression rows on non-CDS set: {len(expr):,}", flush=True)

    expr_ids = expr["orf_id"].to_numpy()
    reads = expr[read_cols].to_numpy(dtype=float)
    pns = expr[pn_cols].to_numpy(dtype=float)
    del expr

    rows = []
    for T in PN_SWEEP:
        passm = (reads > READS_MIN) & (pns > T)
        n_pass = int(passm.any(axis=1).sum())
        rows.append((T, n_pass, round(100.0 * n_pass / n_keep, 2)))
    tbl = pd.DataFrame(rows, columns=["pN_threshold", "n_orf_pass_any_sample",
                                      "pct_of_nonCDS"])
    tbl.to_csv(f"{a.out_dir}/stage1_threshold_table.tsv", sep="\t", index=False)
    print("stage1_threshold_table.tsv:")
    print(tbl.to_string(index=False), flush=True)

    m9 = reads > READS_MIN
    out = pd.DataFrame({"orf_id": expr_ids,
                        "n_samples_reads_gt9": m9.sum(axis=1)})
    out = out[out["n_samples_reads_gt9"] > 0]
    out.to_csv(f"{a.out_dir}/{a.prefix}_stage1_passed.tsv", sep="\t", index=False)
    print(f"Wrote {a.prefix}_stage1_passed.tsv: {len(out):,} ORFs "
          f"(reads > {READS_MIN} in >=1 sample, non-CDS set)", flush=True)

    # ── pN distribution ───────────────────────────────────────────────────
    vals = pns[(reads > READS_MIN) & (pns > 0)]
    print(f"  per-sample pN values with reads>{READS_MIN}: {vals.size:,}", flush=True)
    qv = np.quantile(vals, QUANTILES)
    summary = {f"q{int(q * 100)}": round(float(v), 4)
               for q, v in zip(QUANTILES, qv)}
    summary.update({"n": int(vals.size), "mean": round(float(vals.mean()), 4),
                    "max": round(float(vals.max()), 4)})
    with open(f"{a.out_dir}/pn_values_summary.json", "w") as f:
        json.dump(summary, f, indent=2)
    print(json.dumps(summary, indent=2), flush=True)

    fig, ax = plt.subplots(figsize=(9, 5))
    logv = np.log10(vals)
    ax.hist(logv, bins=80, color="#2E75B6", alpha=0.85)
    for q, v in zip(QUANTILES, qv):
        ax.axvline(np.log10(v), color="#C00000", ls="--", lw=1)
        ax.text(np.log10(v), ax.get_ylim()[1] * 0.98, f"q{int(q * 100)}",
                rotation=90, va="top", ha="right", fontsize=8, color="#C00000")
    ax.set_xlabel(f"log10(pN)  (per-sample, reads > {READS_MIN})")
    ax.set_ylabel("count (per-sample ORF records)")
    ax.set_title(f"{title} per-sample pN distribution — {len(samples)} Ribo-seq "
                 f"samples, reads > {READS_MIN}\nred dashed: q10/25/50/75/90/95/99")
    fig.tight_layout()
    fig.savefig(f"{a.out_dir}/pn_distribution.png", dpi=150)
    print(f"Wrote {a.out_dir}/pn_distribution.png", flush=True)
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
