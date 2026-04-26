#!/usr/bin/env python3
"""
collect_qc_stats.py - Aggregate per-sample QC metrics for nf-core/riboseq.

Parses:
  - STAR *Log.final.out         → alignment statistics
  - *.sorf.filter_stats.tsv     → sORF filter statistics
  - *_P_sites_calcs             → RiboseQC P-site periodicity
  - *_all.txt                   → Ribo-TISH ORF counts
  - *_translating_ORFs.tsv      → Ribotricer ORF counts
  - *_final_ORFquant_results    → ORFquant ORF counts

Outputs (in --output_dir):
  alignment_stats.csv, sorf_stats.csv, psite_periodicity_stats.csv,
  orf_counts.csv, qc_summary.csv
"""

import argparse
import csv
import glob
import os
import re
import sys


# ─── Sample ID extraction ────────────────────────────────────────────────────

def _sample_id(path: str, suffix: str) -> str:
    """Extract sample ID by stripping suffix from filename."""
    name = os.path.basename(path)
    if name.endswith(suffix):
        return name[: -len(suffix)]
    return name


# ─── Parsers ─────────────────────────────────────────────────────────────────

def parse_star_log(path: str) -> dict:
    """Parse STAR *Log.final.out into a flat dict."""
    d = {}
    key_map = {
        "Number of input reads": "total_reads",
        "Uniquely mapped reads number": "unique_mapped",
        "Uniquely mapped reads %": "unique_pct",
        "% of reads mapped to multiple loci": "multi_pct",
        "% of reads mapped to too many loci": "multi_too_many_pct",
        "% of reads unmapped: too many mismatches": "unmapped_mismatch_pct",
        "% of reads unmapped: too short": "unmapped_short_pct",
        "% of reads unmapped: other": "unmapped_other_pct",
        "Number of splices: Total": "splices_total",
        "Average mapped length": "avg_mapped_length",
        "Average input read length": "avg_input_length",
    }
    try:
        with open(path) as fh:
            for line in fh:
                line = line.strip()
                if "|" not in line:
                    continue
                k, v = [x.strip() for x in line.split("|", 1)]
                if k in key_map:
                    col = key_map[k]
                    v = v.rstrip("%").strip()
                    try:
                        d[col] = int(v) if "." not in v else float(v)
                    except ValueError:
                        d[col] = v
    except Exception as e:
        print(f"WARNING: Could not parse STAR log {path}: {e}", file=sys.stderr)
    return d


def parse_sorf_stats(path: str) -> dict:
    """Parse *.sorf.filter_stats.tsv."""
    try:
        with open(path) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                return {
                    "total_primary_mapped": int(row["total_primary_mapped"]),
                    "kept_primary_mapped": int(row["kept_primary_mapped"]),
                    "pct_kept": float(row["pct_kept"]),
                }
    except Exception as e:
        print(f"WARNING: Could not parse sORF stats {path}: {e}", file=sys.stderr)
    return {}


def parse_psites_calcs(path: str) -> list:
    """Parse *_P_sites_calcs; returns list of row dicts."""
    rows = []
    try:
        with open(path) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                rows.append({
                    "read_length": int(row["read_length"]),
                    "offset": int(row["cutoff"]),
                    "frame_preference": round(float(row["frame_preference"]), 4),
                    "pct_map": float(row["pct_map"]),
                    "max_coverage": row["max_coverage"],
                    "max_inframe": row["max_inframe"],
                    "all": row["all"],
                })
    except Exception as e:
        print(f"WARNING: Could not parse P_sites_calcs {path}: {e}", file=sys.stderr)
    return rows


def count_lines_minus_header(path: str) -> int:
    """Count data lines (excluding header) in a text file."""
    try:
        with open(path) as fh:
            return sum(1 for ln in fh if ln.strip() and not ln.startswith("#")) - 1
    except Exception:
        return 0


def psite_summary(rows: list) -> dict:
    """Derive per-sample P-site summary from parsed psites rows."""
    max_cov = [r for r in rows if r["max_coverage"] == "TRUE"]
    if not max_cov:
        return {}
    dominant = max(max_cov, key=lambda r: r["pct_map"])
    return {
        "dominant_rl": dominant["read_length"],
        "dominant_offset": dominant["offset"],
        "dominant_frame_pref": dominant["frame_preference"],
        "best_frame_pref": round(max(r["frame_preference"] for r in max_cov), 4),
        "n_valid_rl": len(max_cov),
    }


# ─── File discovery ───────────────────────────────────────────────────────────

def discover_files(paths: list, suffix: str) -> dict:
    """
    Given a list of file paths (possibly with globs), return a
    {sample_id: path} dict for files matching *<suffix>.
    """
    result = {}
    for p in (paths or []):
        for match in glob.glob(p) if "*" in p else [p]:
            if not os.path.isfile(match):
                continue
            sample = _sample_id(match, suffix)
            result[sample] = match
    return result


# ─── Main ─────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--star_logs",       nargs="*", default=[], help="STAR *Log.final.out files")
    ap.add_argument("--sorf_stats",      nargs="*", default=[], help="*.sorf.filter_stats.tsv files")
    ap.add_argument("--psites_calcs",    nargs="*", default=[], help="*_P_sites_calcs files")
    ap.add_argument("--ribotish_all",    nargs="*", default=[], help="Ribo-TISH *_all.txt files")
    ap.add_argument("--ribotricer_orfs", nargs="*", default=[], help="Ribotricer *_translating_ORFs.tsv files")
    ap.add_argument("--ribocode",        nargs="*", default=[], help="RiboCode *_collapsed.txt files")
    ap.add_argument("--orfquant",        nargs="*", default=[], help="ORFquant *_final_ORFquant_results files")
    ap.add_argument("--output_dir",      default=".",          help="Output directory (default: .)")
    args = ap.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    # ── Discover files by suffix ──────────────────────────────────────────────
    star_map    = discover_files(args.star_logs,       ".Log.final.out")
    sorf_map    = discover_files(args.sorf_stats,      ".sorf.filter_stats.tsv")
    psite_map   = discover_files(args.psites_calcs,    "_P_sites_calcs")
    rtish_map   = discover_files(args.ribotish_all,    "_all.txt")
    rtricer_map = discover_files(args.ribotricer_orfs, "_translating_ORFs.tsv")
    ribocode_map = discover_files(args.ribocode,       "_collapsed.txt")
    orfq_map    = discover_files(args.orfquant,        "_final_ORFquant_results")

    all_samples = sorted(set(
        list(star_map) + list(sorf_map) + list(psite_map) +
        list(rtish_map) + list(rtricer_map) + list(ribocode_map) + list(orfq_map)
    ))

    if not all_samples:
        print("WARNING: No input files found. Producing empty output files.", file=sys.stderr)

    print(f"Collecting QC stats for {len(all_samples)} samples", file=sys.stderr)

    # ── Build per-sample data ─────────────────────────────────────────────────
    aln_rows    = []
    sorf_rows   = []
    psite_rows  = []
    orf_rows    = []

    for s in all_samples:
        # --- Alignment ---
        if s in star_map:
            d = parse_star_log(star_map[s])
            if d:
                aln_rows.append({"sample": s, **d})

        # --- sORF filter ---
        if s in sorf_map:
            d = parse_sorf_stats(sorf_map[s])
            if d:
                sorf_rows.append({"sample": s, **d})

        # --- P-site periodicity ---
        if s in psite_map:
            for r in parse_psites_calcs(psite_map[s]):
                psite_rows.append({"sample": s, **r})

        # --- ORF counts ---
        orf_row = {"sample": s}
        orf_row["ribotish_all"]   = count_lines_minus_header(rtish_map[s])   if s in rtish_map   else ""
        orf_row["ribotricer"]     = count_lines_minus_header(rtricer_map[s]) if s in rtricer_map else ""
        orf_row["ribocode"]       = count_lines_minus_header(ribocode_map[s]) if s in ribocode_map else ""
        orf_row["orfquant"]       = count_lines_minus_header(orfq_map[s])    if s in orfq_map    else ""
        ps_sum = psite_summary(parse_psites_calcs(psite_map[s])) if s in psite_map else {}
        orf_row.update(ps_sum)
        orf_rows.append(orf_row)

    # ── Write CSV files ───────────────────────────────────────────────────────
    def write_csv(name, rows, fields):
        out = os.path.join(args.output_dir, name)
        with open(out, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=fields)
            w.writeheader()
            for r in rows:
                w.writerow({k: r.get(k, "") for k in fields})
        print(f"  Wrote {out} ({len(rows)} rows)", file=sys.stderr)

    write_csv("alignment_stats.csv", aln_rows, [
        "sample", "total_reads", "unique_mapped", "unique_pct", "multi_pct",
        "multi_too_many_pct", "unmapped_mismatch_pct", "unmapped_short_pct",
        "unmapped_other_pct", "avg_mapped_length", "avg_input_length", "splices_total",
    ])

    write_csv("sorf_stats.csv", sorf_rows, [
        "sample", "total_primary_mapped", "kept_primary_mapped", "pct_kept",
    ])

    write_csv("psite_periodicity_stats.csv", psite_rows, [
        "sample", "read_length", "offset", "frame_preference", "pct_map",
        "max_coverage", "max_inframe", "all",
    ])

    write_csv("orf_counts.csv", orf_rows, [
        "sample", "ribotish_all", "ribotricer", "ribocode", "orfquant",
        "dominant_rl", "dominant_offset", "dominant_frame_pref", "best_frame_pref", "n_valid_rl",
    ])

    # ── Write merged qc_summary.csv ───────────────────────────────────────────
    aln_d   = {r["sample"]: r for r in aln_rows}
    sorf_d  = {r["sample"]: r for r in sorf_rows}
    orf_d   = {r["sample"]: r for r in orf_rows}
    psite_best = {}
    for r in psite_rows:
        s = r["sample"]
        if r["max_coverage"] == "TRUE":
            if s not in psite_best or r["pct_map"] > psite_best[s]["pct_map"]:
                psite_best[s] = r

    summary_rows = []
    for s in all_samples:
        row = {"sample": s}
        a = aln_d.get(s, {})
        sf = sorf_d.get(s, {})
        o = orf_d.get(s, {})
        pb = psite_best.get(s, {})
        row.update({
            "total_reads":        a.get("total_reads", ""),
            "unique_pct":         a.get("unique_pct", ""),
            "multi_pct":          a.get("multi_pct", ""),
            "unmapped_short_pct": a.get("unmapped_short_pct", ""),
            "total_primary_mapped":  sf.get("total_primary_mapped", ""),
            "kept_primary_mapped":   sf.get("kept_primary_mapped", ""),
            "pct_kept":              sf.get("pct_kept", ""),
            "dominant_rl":           o.get("dominant_rl", ""),
            "dominant_offset":       o.get("dominant_offset", ""),
            "best_frame_pref":       o.get("best_frame_pref", ""),
            "n_valid_rl":            o.get("n_valid_rl", ""),
            "ribotish_all":          o.get("ribotish_all", ""),
            "ribotricer":            o.get("ribotricer", ""),
            "ribocode":              o.get("ribocode", ""),
            "orfquant":              o.get("orfquant", ""),
        })
        summary_rows.append(row)

    write_csv("qc_summary.csv", summary_rows, [
        "sample",
        "total_reads", "unique_pct", "multi_pct", "unmapped_short_pct",
        "total_primary_mapped", "kept_primary_mapped", "pct_kept",
        "dominant_rl", "dominant_offset", "best_frame_pref", "n_valid_rl",
        "ribotish_all", "ribotricer", "ribocode", "orfquant",
    ])

    print(f"\nDone. Output written to: {args.output_dir}", file=sys.stderr)


if __name__ == "__main__":
    main()
