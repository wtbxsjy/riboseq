#!/usr/bin/env python3
"""Post-analysis Step 2: P-site purity filter + length cut + top-N% per biotype.

Generic version of the per-project scripts (PRJEB26593 / GSE120762).

Stage 2 (per sample, pass if >= 1 sample):
    {sample}_p_site_GSE > 9
    AND pct_s >= 0.5,  where pct_s = {sample}_p_site_GSE / {sample}_reads
    AND global_p_site_pos > 2 (ORF-level, weighted avg P-site position from ORF start)

Why the expression-summary reads as the denominator:
    compute_psite_purity.py pairs raw P-site counts with the RiboseQC COVERAGE
    bedgraphs, which are RPM-normalized (genome-wide integral == 1e6) -> its
    p_site_pct column is scale-invalid. The expression summary's {sample}_reads
    are RAW counts on the same scale as the P-site bedgraphs, so the ratio is
    meaningful (empirical q50 = 0.49 on PRJEB26593, see CLAUDE.md gotcha 28).

    Also note compute_psite_purity.py's per-sample {sample}_p_site_pos column
    divides the value-sum by psite -> always 1.0. The correct ORF-level metric
    is `global_p_site_pos`, which is what the position criterion uses.

Step 4: length_aa <= 50
Step 5: top-N% per biotype (gencode orf_biotype, fallback orf_type_category),
        ranked by purity total_psites.

Outputs (in --out-dir):
  {prefix}_stage2_passed.tsv      Stage-2 survivors (full columns)
  {prefix}_filtered_len50aa.tsv   Stage-2 survivors with length_aa <= 50
  {prefix}_top10_per_biotype.tsv  top-N% per biotype of the length-cut set
  {prefix}_filter_stats.json      counts + per-biotype breakdowns

Usage:
  python3 stage2_psite_filter.py \
      --purity  post_analysis/<project>/<prefix>_psite_purity.tsv \
      --stage1  post_analysis/<project>/<prefix>_stage1_passed.tsv \
      --noncds  post_analysis/<project>/step1_nonCDS.tsv.gz \
      --expression <result>/orf_unification/unified_orfs_expression_summary.tsv \
      --out-dir post_analysis/<project> --prefix <project>
"""

import argparse
import json

import numpy as np
import pandas as pd

PSITE_MIN = 9        # per-sample P-site count floor
PCT_MIN = 0.5        # per-sample p_site_GSE / reads floor
POS_MIN = 2          # ORF-level global_p_site_pos floor
LEN_AA_MAX = 50      # Step-4 length cut
TOP_FRAC = 0.10      # Step-5 fraction per biotype


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--purity", required=True)
    ap.add_argument("--stage1", required=True)
    ap.add_argument("--noncds", required=True)
    ap.add_argument("--expression", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--prefix", required=True)
    a = ap.parse_args()

    st1 = pd.read_csv(a.stage1, sep="\t", dtype={"orf_id": str})
    keep = set(st1["orf_id"])
    print(f"stage1 set: {len(keep):,}", flush=True)

    noncds = pd.read_csv(
        a.noncds, compression="gzip", sep="\t", low_memory=False,
        usecols=["orf_id", "length_aa", "orf_type_category", "gencode_biotype",
                 "total_psites", "total_reads", "pN"],
        dtype={"orf_id": str})
    noncds = noncds.rename(columns={"total_psites": "total_psites_ortype",
                                    "total_reads": "total_reads_ortype",
                                    "pN": "pN_ortype"})

    print("Loading expression summary reads ...", flush=True)
    hdr = pd.read_csv(a.expression, sep="\t", nrows=0).columns.tolist()
    rc = [c for c in hdr if c.endswith("_reads") and c != "total_reads"]
    expr = pd.read_csv(a.expression, sep="\t", usecols=["orf_id"] + rc,
                       dtype={"orf_id": str})
    expr = expr[expr["orf_id"].isin(keep)]
    print(f"  expression rows on stage1 set: {len(expr):,}", flush=True)

    print("Loading purity matrix (chunked, stage1 subset) ...", flush=True)
    pheader = pd.read_csv(a.purity, sep="\t", nrows=0).columns.tolist()
    sample_cols = sorted({c[: -len("_p_site_GSE")] for c in pheader
                          if c.endswith("_p_site_GSE")
                          and not c.endswith("_not_p_site_GSE")})
    usecols = (["orf_id", "total_psites", "global_p_site_pos",
                "n_samples_with_psites"]
               + [f"{s}_p_site_GSE" for s in sample_cols])
    chunks = pd.read_csv(a.purity, sep="\t", dtype={"orf_id": str},
                         usecols=usecols, chunksize=500_000, low_memory=False)
    frames = []
    for ch in chunks:
        ch = ch[ch["orf_id"].isin(keep)]
        if len(ch):
            frames.append(ch)
    pur = pd.concat(frames, ignore_index=True)
    print(f"  purity rows on stage1 set: {len(pur):,}", flush=True)

    rc_ordered = [f"{s}_reads" for s in sample_cols]
    missing = [c for c in rc_ordered if c not in expr.columns]
    if missing:
        raise SystemExit(f"missing expression read columns: {missing}")

    m = pur.merge(expr, on="orf_id", how="left")
    psite = m[[f"{s}_p_site_GSE" for s in sample_cols]].to_numpy(float)
    reads = m[rc_ordered].to_numpy(float)

    pct = np.where(reads > 0, psite / np.maximum(reads, 1e-9), 0.0)
    pass_s = (psite > PSITE_MIN) & (pct >= PCT_MIN)
    n_pass = pass_s.sum(axis=1)

    pass_cols = [pd.Series(np.where(pass_s[:, i], s, ""), index=m.index)
                 for i, s in enumerate(sample_cols)]
    m["n_samples_pass_stage2"] = n_pass
    m["pass_samples_stage2"] = (pd.concat(pass_cols, axis=1)
                                  .agg(",".join, axis=1).str.strip(","))

    psite9 = int((psite > PSITE_MIN).any(axis=1).sum())
    pctpass = int((n_pass > 0).sum())
    pos_pass = int(((m["global_p_site_pos"] > POS_MIN) & (n_pass > 0)).sum())
    print(f"  any-sample psite>{PSITE_MIN}: {psite9:,}", flush=True)
    print(f"  psite>{PSITE_MIN} AND pct>={PCT_MIN}: {pctpass:,}", flush=True)
    print(f"  + global_p_site_pos>{POS_MIN}: {pos_pass:,}", flush=True)

    stage2 = m[(n_pass > 0) & (m["global_p_site_pos"] > POS_MIN)].copy()
    print(f"  stage2 passed: {len(stage2):,} "
          f"({100.0 * len(stage2) / len(keep):.1f}% of stage1)", flush=True)

    stage2 = stage2.merge(noncds, on="orf_id", how="left")
    stage2["biotype_final"] = (stage2["gencode_biotype"]
                               .fillna(stage2["orf_type_category"]))

    len50 = stage2[stage2["length_aa"] <= LEN_AA_MAX].copy()
    print(f"  length <= {LEN_AA_MAX} aa: {len(len50):,}", flush=True)

    top_parts, top_counts = [], {}
    for bt, grp in len50.groupby("biotype_final", sort=False):
        grp = grp.sort_values("total_psites", ascending=False)
        k = max(1, int(np.ceil(TOP_FRAC * len(grp))))
        top_parts.append(grp.head(k))
        top_counts[bt] = k
    top10 = (pd.concat(top_parts, ignore_index=True)
             if top_parts else len50.iloc[0:0].copy())
    print(f"  top {int(TOP_FRAC*100)}% per biotype: {len(top10):,} ORFs across "
          f"{len(top_counts)} biotypes", flush=True)

    stage2.to_csv(f"{a.out_dir}/{a.prefix}_stage2_passed.tsv", sep="\t", index=False)
    len50.to_csv(f"{a.out_dir}/{a.prefix}_filtered_len50aa.tsv", sep="\t", index=False)
    top10.to_csv(f"{a.out_dir}/{a.prefix}_top10_per_biotype.tsv", sep="\t", index=False)

    stats = {
        "stage1_passed": len(st1),
        "any_sample_psite_gt9": psite9,
        "any_sample_psite_gt9_and_pct_ge_0.5": pctpass,
        "stage2_passed": int(len(stage2)),
        "stage2_pct_of_stage1": round(100.0 * len(stage2) / len(keep), 2),
        "length_le50aa": int(len(len50)),
        "top10_per_biotype": int(len(top10)),
        "n_biotypes": len(top_counts),
        "pct_definition": "p_site_GSE / {sample}_reads (raw counts, same scale)",
        "pos_criterion": "global_p_site_pos > 2 (ORF-level)",
        "stage2_biotype_counts": stage2["biotype_final"]
                                        .value_counts(dropna=False).to_dict(),
        "len50_biotype_counts": len50["biotype_final"]
                                       .value_counts(dropna=False).to_dict(),
        "top10_biotype_counts": top_counts,
    }
    with open(f"{a.out_dir}/{a.prefix}_filter_stats.json", "w") as f:
        json.dump(stats, f, indent=2)
    print(json.dumps(stats, indent=2), flush=True)
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
