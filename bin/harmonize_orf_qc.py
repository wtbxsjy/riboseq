#!/usr/bin/env python3
"""
harmonize_orf_qc.py — Phase 2 of ORF QC Module

Harmonizes QC metrics across all tools into a unified schema:
  1. P-site offset harmonization (cross-tool consensus)
  2. Periodicity score normalization (0-1 scale)
  3. Read-level QC summary

Input:  tool_data.json (from extract_orf_qc_metrics.py)
Output: psite_harmonized.tsv, periodicity.tsv, read_qc_summary.json

Usage:
  harmonize_orf_qc.py --input tool_data.json --output-prefix SAMPLE
"""

import argparse
import json
import sys
from collections import defaultdict
from typing import Any, Dict, List, Optional, Tuple


# ---------------------------------------------------------------------------
# P-site offset harmonization
# ---------------------------------------------------------------------------

# Tool priority order for P-site authority (higher = more trusted)
PSITE_AUTHORITY_ORDER = [
    "ribowaltz",    # Two-step coherence correction — best algorithm
    "riboseqc",     # Peak-based with frame distribution validation
    "ribocode",     # Wilcoxon-based periodicity + peak
    "ribotricer",   # Cross-correlation method
    "rpbp",         # Bayes factor peak selection
    "ribotish",     # Simple aggregate profile peak
    "price",        # GEDI internal method
]


def extract_psite_offsets(tool_data: Dict) -> Dict[str, Dict[int, Dict[str, Any]]]:
    """Extract per-tool, per-read-length P-site offsets.
    Returns: {tool_name: {read_length: {offset: int, ...}}}
    """
    psites: Dict[str, Dict[int, Dict]] = defaultdict(dict)

    for tool_name, samples in tool_data.get("samples", {}).items():
        for sample_id, sample_data in samples.items():
            if sample_data.get("status") != "OK":
                continue

            # Different tools store P-site info differently
            if tool_name == "ribowaltz":
                for p in sample_data.get("psites", []):
                    rl = p.get("read_length")
                    if rl:
                        psites[tool_name][rl] = {
                            "offset": p.get("corrected_offset_from_5"),
                            "total_pct": p.get("total_percentage"),
                            "start_pct": p.get("start_percentage"),
                        }

            elif tool_name == "riboseqc":
                for p in sample_data.get("psites", []):
                    rl = p.get("read_length")
                    if rl:
                        psites[tool_name][rl] = {
                            "offset": p.get("p_site_offset"),
                            "f0_pct": p.get("f0_percent"),
                        }

            elif tool_name == "ribocode":
                # RiboCode stores P-site in _pre_config.txt which isn't parsed
                # by the extractor yet; ORF-level frame data is available
                # From ORFs, we can infer typical offsets
                pass

            elif tool_name == "ribotish":
                offsets = sample_data.get("psite_offsets", {})
                for rl, offset in offsets.items():
                    psites[tool_name][int(rl)] = {"offset": offset}

            elif tool_name == "rpbp":
                for p in sample_data.get("psite_offsets", []):
                    rl = p.get("read_length")
                    if rl:
                        psites[tool_name][rl] = {
                            "offset": p.get("p_site_offset"),
                            "bf_mean": p.get("bf_mean"),
                        }

    return dict(psites)


def harmonize_psites(
    psites: Dict[str, Dict[int, Dict]],
    max_delta: int = 1
) -> Tuple[List[Dict], Dict[str, Any]]:
    """Harmonize P-site offsets across tools.

    Returns:
        table: List of per-read-length rows for TSV output.
        summary: Dict with consensus statistics.
    """
    # Collect all read lengths across all tools
    all_lengths: set = set()
    for tool_offsets in psites.values():
        all_lengths.update(tool_offsets.keys())
    all_lengths = sorted(all_lengths)

    if not all_lengths:
        return [], {"status": "NO_DATA"}

    table = []
    concordant = 0
    total_with_data = 0

    for rl in all_lengths:
        row: Dict[str, Any] = {"read_length": rl}
        offsets = {}

        for tool_name in PSITE_AUTHORITY_ORDER:
            tool_data = psites.get(tool_name, {}).get(rl, {})
            offset = tool_data.get("offset")
            row[f"{tool_name}_offset"] = offset
            if offset is not None:
                offsets[tool_name] = offset

        if len(offsets) >= 2:
            total_with_data += 1
            values = list(offsets.values())
            max_val = max(values)
            min_val = min(values)
            row["n_tools"] = len(offsets)
            row["max_delta"] = max_val - min_val

            # Consensus: median of available values
            sorted_vals = sorted(values)
            n = len(sorted_vals)
            if n % 2 == 1:
                row["consensus_offset"] = sorted_vals[n // 2]
            else:
                row["consensus_offset"] = (sorted_vals[n // 2 - 1] + sorted_vals[n // 2]) // 2

            # Flag
            if row["max_delta"] <= max_delta:
                row["flag"] = "OK"
                concordant += 1
            elif row["max_delta"] <= 2:
                row["flag"] = "WARN"
            else:
                row["flag"] = "FAIL"

        elif len(offsets) == 1:
            total_with_data += 1
            row["n_tools"] = 1
            row["max_delta"] = 0
            row["consensus_offset"] = list(offsets.values())[0]
            row["flag"] = "OK_SINGLE"
            concordant += 1
        else:
            row["n_tools"] = 0
            row["max_delta"] = None
            row["consensus_offset"] = None
            row["flag"] = "NO_DATA"

        # Primary authority offset
        for tool_name in PSITE_AUTHORITY_ORDER:
            if tool_name in offsets:
                row["authority_offset"] = offsets[tool_name]
                row["authority_tool"] = tool_name
                break

        table.append(row)

    consensus_rate = concordant / total_with_data if total_with_data > 0 else 0

    summary = {
        "total_read_lengths": len(all_lengths),
        "lengths_with_data": total_with_data,
        "lengths_concordant": concordant,
        "consensus_rate": round(consensus_rate, 3),
        "status": "OK" if consensus_rate >= 0.7 else "WARN" if consensus_rate >= 0.4 else "FAIL",
    }

    return table, summary


# ---------------------------------------------------------------------------
# Periodicity harmonization
# ---------------------------------------------------------------------------

def normalize_periodicity(raw_value: float, tool_name: str, metric: str = "periodicity") -> float:
    """Normalize a tool-specific periodicity metric to [0, 1] scale.

    Higher = more periodic.
    """
    if tool_name in ("ribocode", "riboseqc", "ribowaltz"):
        # f0_percent or frame-0% — already 0-1
        return min(1.0, max(0.0, raw_value))

    elif tool_name == "ribotricer":
        # phase_score — already 0-1 (Fourier coherence)
        return min(1.0, max(0.0, raw_value))

    elif tool_name == "rpbp":
        # Bayes factor — log10(BF) / 3, clamp to [0,1]
        if raw_value <= 0:
            return 0.0
        import math
        return min(1.0, math.log10(raw_value + 1) / 3.0)

    elif tool_name == "ribotish":
        # FrameQvalue — convert to 1 - q (q越小越periodic)
        return min(1.0, max(0.0, 1.0 - raw_value))

    elif tool_name == "orfquant":
        # p-value — convert via -log10(p)/10
        if raw_value <= 0:
            return 1.0
        import math
        return min(1.0, -math.log10(raw_value) / 10.0)

    else:
        return 0.5  # Unknown — neutral


def assess_read_periodicity(tool_data: Dict) -> Dict[str, Any]:
    """Assess sample-wide read-level periodicity from all available tools.

    Returns:
        Dict with per-tool periodicity scores and aggregate summary.
    """
    periodicity: Dict[str, Any] = {
        "per_tool": {},
        "aggregate_score": None,
        "classification": "UNKNOWN",
    }

    scores = []

    for tool_name, samples in tool_data.get("samples", {}).items():
        if tool_name == "price":
            continue  # PRICE doesn't report read-level periodicity

        tool_scores = []
        for sample_id, sample_data in samples.items():
            if sample_data.get("status") != "OK":
                continue

            if tool_name == "ribowaltz":
                # Frame distribution gives per-region frame %
                # For now, use P-site total_percentage as proxy
                for p in sample_data.get("psites", []):
                    total_pct = p.get("total_percentage", 0)
                    if total_pct > 0:
                        tool_scores.append(min(1.0, total_pct / 10.0))

            elif tool_name == "riboseqc":
                for p in sample_data.get("psites", []):
                    f0 = p.get("f0_percent", 0)
                    if f0 > 0:
                        tool_scores.append(f0 / 100.0 if f0 > 1 else f0)

            elif tool_name == "ribocode":
                for orf in sample_data.get("orfs", [])[:100]:  # Top 100 ORFs
                    pval = orf.get("pval_combined", 1)
                    if pval < 1:
                        import math
                        tool_scores.append(min(1.0, -math.log10(max(pval, 1e-10)) / 10.0))

            elif tool_name == "ribotricer":
                for orf in sample_data.get("orfs", [])[:100]:
                    ps = orf.get("phase_score", 0)
                    tool_scores.append(ps)

            elif tool_name == "rpbp":
                for p in sample_data.get("psite_offsets", []):
                    bf = p.get("bf_mean", 0)
                    if bf > 0:
                        import math
                        tool_scores.append(min(1.0, math.log10(bf + 1) / 3.0))

        if tool_scores:
            mean_score = sum(tool_scores) / len(tool_scores)
            periodicity["per_tool"][tool_name] = round(mean_score, 3)
            scores.append(mean_score)

    if scores:
        periodicity["aggregate_score"] = round(sum(scores) / len(scores), 3)
        agg = periodicity["aggregate_score"]
        if agg >= 0.65:
            periodicity["classification"] = "GOOD"
        elif agg >= 0.5:
            periodicity["classification"] = "MARGINAL"
        else:
            periodicity["classification"] = "POOR"

    return periodicity


# ---------------------------------------------------------------------------
# Sample QC flags
# ---------------------------------------------------------------------------

def compute_sample_flags(
    psite_summary: Dict,
    periodicity: Dict,
    tool_data: Dict,
) -> Dict[str, Any]:
    """Compute sample-level quality flags."""
    flags = []

    # Periodicity flags
    agg_score = periodicity.get("aggregate_score")
    if agg_score is not None:
        if agg_score < 0.3:
            flags.append({
                "flag": "NO_PERIODICITY",
                "severity": "CRITICAL",
                "detail": f"Aggregate periodicity score {agg_score:.2f} — no read length shows clear 3-nt periodicity",
            })
        elif agg_score < 0.5:
            flags.append({
                "flag": "LOW_PERIODICITY",
                "severity": "WARNING",
                "detail": f"Aggregate periodicity score {agg_score:.2f} — <50% read lengths periodic",
            })

    # P-site discordance
    if psite_summary.get("status") == "FAIL":
        flags.append({
            "flag": "P_SITE_DISCORDANCE",
            "severity": "WARNING",
            "detail": f"Only {psite_summary.get('consensus_rate', 0):.0%} read lengths have consistent P-site offsets",
        })

    # Tool failure
    failed_tools = [
        t for t, s in tool_data.get("tool_status", {}).items()
        if s in ("FAILED", "NOT_FOUND")
    ]
    if failed_tools:
        flags.append({
            "flag": "TOOL_FAILURE",
            "severity": "INFO",
            "detail": f"Tools failed or not found: {', '.join(failed_tools)}",
        })

    # ORF yield
    total_orfs = 0
    for tool_name in tool_data.get("tool_status", {}):
        if tool_data["tool_status"][tool_name] == "OK":
            for sample_data in tool_data["samples"].get(tool_name, {}).values():
                total_orfs += sample_data.get("orf_count", 0)
    if total_orfs < 100:
        flags.append({
            "flag": "LOW_YIELD",
            "severity": "INFO",
            "detail": f"Only {total_orfs} ORFs detected across all tools",
        })

    return {
        "flags": flags,
        "summary": {
            "aggregate_periodicity": agg_score,
            "psite_consensus_rate": psite_summary.get("consensus_rate"),
            "total_orfs_raw": total_orfs,
            "tools_ran": sum(
                1 for s in tool_data.get("tool_status", {}).values()
                if s == "OK"
            ),
            "tools_failed": len(failed_tools),
        },
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Harmonize ORF QC metrics across tools"
    )
    parser.add_argument("--input", "-i", required=True, help="tool_data.json from extract_orf_qc_metrics.py")
    parser.add_argument("--output-prefix", "-o", default="qc", help="Output file prefix")
    parser.add_argument("--max-offset-delta", type=int, default=1, help="Max allowed P-site offset difference")
    args = parser.parse_args()

    with open(args.input) as fh:
        tool_data = json.load(fh)

    # 1. P-site harmonization
    print("Harmonizing P-site offsets...")
    psites = extract_psite_offsets(tool_data)
    psite_table, psite_summary = harmonize_psites(psites, args.max_offset_delta)

    # 2. Periodicity assessment
    print("Assessing read-level periodicity...")
    periodicity = assess_read_periodicity(tool_data)

    # 3. Sample flags
    print("Computing sample flags...")
    flags = compute_sample_flags(psite_summary, periodicity, tool_data)

    # Write outputs
    # P-site table
    if psite_table:
        with open(f"{args.output_prefix}_psite_harmonized.tsv", "w") as fh:
            headers = ["read_length"] + [f"{t}_offset" for t in PSITE_AUTHORITY_ORDER] + \
                      ["n_tools", "consensus_offset", "authority_offset", "authority_tool", "max_delta", "flag"]
            fh.write("\t".join(headers) + "\n")
            for row in psite_table:
                fh.write("\t".join(str(row.get(h, "")) for h in headers) + "\n")
        print(f"  Wrote {args.output_prefix}_psite_harmonized.tsv ({len(psite_table)} read lengths)")
    else:
        print("  WARNING: No P-site data found")

    # Periodicity summary
    with open(f"{args.output_prefix}_periodicity.json", "w") as fh:
        json.dump(periodicity, fh, indent=2)
    print(f"  Periodicity aggregate score: {periodicity.get('aggregate_score', 'N/A')} ({periodicity.get('classification', 'N/A')})")

    # Sample flags
    with open(f"{args.output_prefix}_sample_flags.json", "w") as fh:
        json.dump(flags, fh, indent=2)
    print(f"  Flags: {len(flags['flags'])} issue(s)")

    for f in flags["flags"]:
        print(f"    [{f['severity']}] {f['flag']}: {f['detail']}")


if __name__ == "__main__":
    main()
