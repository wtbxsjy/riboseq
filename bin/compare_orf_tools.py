#!/usr/bin/env python3
"""
compare_orf_tools.py — Phase 3+4 of ORF QC Module

Phase 3: Cross-tool ORF comparison
  - Pairwise overlap and Jaccard index (using bedtools intersect)
  - Classification agreement matrix
  - Consensus ORF identification

Phase 4: ORF Confidence Score (OCS) computation
  - S_translation: normalized significance from detecting tool(s)
  - S_agreement: cross-tool support
  - S_coverage: coverage completeness
  - S_periodicity: ORF-level periodicity
  - S_readlevel: sample-wide read QC modifier

Input:  tool_data.json, unified_orfs.bed, unified_orfs.metadata.tsv, periodicity.json
Output: orf_confidence.tsv, tool_agreement.tsv

Usage:
  compare_orf_tools.py \
    --tool-data tool_data.json \
    --unified-bed unified_orfs.bed \
    --unified-meta unified_orfs.metadata.tsv \
    --periodicity periodicity.json \
    --output-prefix qc
"""

import argparse
import json
import math
import os
import re
import subprocess
import sys
import tempfile
from collections import defaultdict
from typing import Any, Dict, List, Optional, Set, Tuple


# ---------------------------------------------------------------------------
# Phase 3: Cross-Tool Comparison
# ---------------------------------------------------------------------------

def _parse_genome_pos(genome_pos: str) -> tuple:
    """Parse 'chr:start-end:strand' or 'chr+strand:start-end' format.
    Returns (chrom, start, end, strand) or (None, 0, 0, '.').
    """
    try:
        # Format 1: "chr+strand:start-end"  (PRICE)
        if "+:" in genome_pos or "-:" in genome_pos:
            m = re.match(r'([^+:]+)([+-]):(\d+)-(\d+)', genome_pos)
            if m:
                return (m.group(1), int(m.group(3)) - 1, int(m.group(4)), m.group(2))
        # Format 2: "chr:start-end:strand"  (Ribo-TISH)
        m = re.match(r'([^:]+):(\d+)-(\d+):([+-])', genome_pos)
        if m:
            return (m.group(1), int(m.group(2)) - 1, int(m.group(3)), m.group(4))
    except Exception:
        pass
    return (None, 0, 0, ".")


def _orfs_to_bed6(orfs: List[Dict], tool_name: str) -> str:
    """Convert ORF list to temporary BED6 file. Returns file path."""
    tmp = tempfile.NamedTemporaryFile(mode="w", suffix=".bed", delete=False)
    written = 0
    for i, orf in enumerate(orfs):
        chrom = orf.get("chrom", "unknown")
        strand = orf.get("strand", ".")

        # Try multiple coordinate field names based on tool
        start = orf.get("orf_gstart", orf.get("start", None))
        end = orf.get("orf_gstop", orf.get("end", None))

        # If start/end are missing, try parsing genome_pos or location
        if start is None or end is None:
            pos_str = orf.get("genome_pos", "") or orf.get("location", "")
            if pos_str and ":" in str(pos_str):
                c, s, e, st = _parse_genome_pos(str(pos_str))
                if c:
                    chrom = c
                    strand = st
                    start = s + 1  # Convert back to 1-based
                    end = e

        orf_id = orf.get("orf_id", f"{tool_name}_{i}")

        try:
            s = int(start) - 1 if start is not None else 0  # BED is 0-based
            e = int(end) if end is not None else s + 1
            if s < 0:
                s = 0
            if e <= s:
                continue
        except (ValueError, TypeError):
            continue

        tmp.write(f"{chrom}\t{s}\t{e}\t{orf_id}\t.\t{strand}\n".replace("\t\t", "\t.\t"))
        written += 1
    tmp.close()
    if written == 0:
        os.unlink(tmp.name)
        return ""
    return tmp.name


def _bedtools_intersect(bed_a: str, bed_b: str, reciprocal: bool = True) -> List[Tuple[str, str]]:
    """Run bedtools intersect between two BED files.
    Returns list of (orf_id_a, orf_id_b) overlapping pairs.
    """
    cmd = ["bedtools", "intersect", "-wa", "-wb"]
    if reciprocal:
        cmd += ["-f", "0.5", "-r"]
    cmd += ["-a", bed_a, "-b", bed_b]

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
        if result.returncode != 0:
            print(f"  WARNING: bedtools intersect failed: {result.stderr}", file=sys.stderr)
            return []
        pairs = []
        for line in result.stdout.strip().split("\n"):
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) >= 10:
                pairs.append((fields[3], fields[9]))  # name columns (4th and 10th)
        return pairs
    except FileNotFoundError:
        print("  WARNING: bedtools not found — skipping pairwise comparison", file=sys.stderr)
        return []
    except Exception as e:
        print(f"  WARNING: bedtools error: {e}", file=sys.stderr)
        return []


def compute_tool_agreement(tool_data: Dict) -> Dict[str, Any]:
    """Compute pairwise ORF overlap between tools.

    Returns:
        Dict with pairwise_jaccard, pairwise_overlap_count, etc.
    """
    # Collect ORFs per tool (prediction tools only)
    tool_orfs: Dict[str, List[Dict]] = {}
    for tool_name, samples in tool_data.get("samples", {}).items():
        if tool_data.get("tool_status", {}).get(tool_name) != "OK":
            continue
        # Skip pure QC tools
        tool_info = {}
        for t in ["ribocode", "ribotricer", "ribotish", "price", "rpbp", "orfquant"]:
            if t == tool_name:
                tool_info["type"] = "orf_predictor"
        if tool_name in ("riboseqc", "ribowaltz"):
            continue

        all_orfs = []
        for sample_data in samples.values():
            all_orfs.extend(sample_data.get("orfs", []))
        tool_orfs[tool_name] = all_orfs

    if len(tool_orfs) < 2:
        return {"status": "INSUFFICIENT_DATA", "pairwise": []}

    # Convert each tool's ORFs to BED6
    bed_files: Dict[str, str] = {}
    for tool_name, orfs in tool_orfs.items():
        bf = _orfs_to_bed6(orfs, tool_name)
        if bf:
            bed_files[tool_name] = bf

    # Compute pairwise comparisons
    tool_names = sorted(tool_orfs.keys())
    pairwise = []

    for i, tool_a in enumerate(tool_names):
        for j, tool_b in enumerate(tool_names):
            if i >= j:
                continue
            if tool_a not in bed_files or tool_b not in bed_files:
                continue
            pairs = _bedtools_intersect(bed_files[tool_a], bed_files[tool_b])
            overlap = len(pairs)
            total_a = len(tool_orfs[tool_a])
            total_b = len(tool_orfs[tool_b])
            union = total_a + total_b - overlap
            jaccard = overlap / union if union > 0 else 0.0

            pairwise.append({
                "tool_a": tool_a,
                "tool_b": tool_b,
                "jaccard": round(jaccard, 4),
                "overlap_count": overlap,
                "a_total": total_a,
                "b_total": total_b,
                "a_only": total_a - overlap,
                "b_only": total_b - overlap,
            })

    # Clean up temp files
    for f in bed_files.values():
        try:
            os.unlink(f)
        except OSError:
            pass

    # Compute classification agreement where possible
    classification_agreement = _compute_classification_agreement(tool_orfs)

    return {
        "status": "OK",
        "tools_compared": len(tool_names),
        "tool_names": tool_names,
        "pairwise": pairwise,
        "classification_agreement": classification_agreement,
    }


def _compute_classification_agreement(tool_orfs: Dict[str, List[Dict]]) -> Dict[str, Any]:
    """Compare ORF classifications across tools."""
    # Collect classification counts per tool
    tool_types: Dict[str, Dict[str, int]] = {}
    for tool_name, orfs in tool_orfs.items():
        counts = defaultdict(int)
        for orf in orfs:
            orf_type = orf.get("orf_type", "unknown")
            # Map to harmonized categories
            harmonized = _harmonize_orf_type(orf_type)
            counts[harmonized] += 1
        tool_types[tool_name] = dict(counts)

    return {"per_tool_types": tool_types}


# ---------------------------------------------------------------------------
# Harmonized ORF type mapping
# ---------------------------------------------------------------------------

ORF_TYPE_MAP = {
    "annotated": "CDS",
    "CDS": "CDS",
    "cds": "CDS",
    "uORF": "uORF",
    "dORF": "dORF",
    "Overlap_uORF": "Overlap_uORF",
    "overlap_uORF": "Overlap_uORF",
    "Overlap_dORF": "Overlap_dORF",
    "overlap_dORF": "Overlap_dORF",
    "internal": "Internal",
    "novel": "Novel",
    "ncRNA": "Novel",
    "Novel": "Novel",
    "Intergenic": "Novel",
}


def _harmonize_orf_type(raw: str) -> str:
    """Map tool-specific ORF type to harmonized category."""
    return ORF_TYPE_MAP.get(raw, raw)


# ---------------------------------------------------------------------------
# Phase 4: ORF Confidence Score (OCS)
# ---------------------------------------------------------------------------

def compute_ocs(
    unified_meta_path: str,
    tool_data: Dict,
    periodicity_data: Dict,
    weights: Tuple[float, ...] = (0.30, 0.30, 0.20, 0.15, 0.05),
) -> List[Dict]:
    """Compute ORF Confidence Score for each unified ORF.

    Args:
        unified_meta_path: Path to unified_orfs.metadata.tsv
        tool_data: Parsed tool data from extract_orf_qc_metrics.py
        periodicity_data: Read-level periodicity assessment
        weights: (w_translation, w_agreement, w_coverage, w_periodicity, w_readlevel)

    Returns:
        List of per-ORF confidence records.
    """
    # Read unified metadata
    unified_orfs = _read_unified_metadata(unified_meta_path)
    if not unified_orfs:
        return []

    # Count available ORF prediction tools from tool_data
    available_tools = [
        t for t, s in tool_data.get("tool_status", {}).items()
        if s == "OK" and t not in ("riboseqc", "ribowaltz")
    ]
    n_available = len(available_tools)

    # Read-level modifier
    s_readlevel = _compute_readlevel_modifier(periodicity_data)

    results = []
    for orf in unified_orfs:
        orf_id = orf.get("orf_id", "")

        # ---- S_agreement: from unified metadata 'tools' column ----
        detecting_tools = _parse_detecting_tools(
            orf.get("tools", ""), orf.get("sources", ""))
        n_detecting = len(detecting_tools)
        if n_available > 0:
            agreement_raw = n_detecting / n_available
            if agreement_raw >= 0.8:    s_agreement = 1.0
            elif agreement_raw >= 0.6:  s_agreement = 0.75
            elif agreement_raw >= 0.4:  s_agreement = 0.5
            elif agreement_raw >= 0.2:  s_agreement = 0.25
            else:                       s_agreement = 0.1
        else:
            s_agreement = 0.0

        # ---- S_translation: from unified metadata 'tool_scores' column ----
        # Format: "Ribotricer:0.678" or "Ribo-TISH:0.5,Ribotricer:0.8"
        best_translation = 0.0
        tool_scores_str = orf.get("tool_scores", "")
        if tool_scores_str and tool_scores_str != "NA":
            for part in tool_scores_str.split(","):
                part = part.strip()
                if ":" in part:
                    try:
                        score = float(part.split(":", 1)[1])
                        if score > best_translation:
                            best_translation = score
                    except (ValueError, IndexError):
                        pass

        # ---- S_periodicity: from unified metadata 'pN' column ----
        # pN is already a 0-1 P-site periodicity score
        best_periodicity = 0.0
        try:
            pn_val = float(orf.get("pN", 0))
            if pn_val > 0:
                best_periodicity = min(1.0, pn_val)
        except (ValueError, TypeError):
            pass

        # ---- S_coverage: from unified metadata 'unique_psites' / length_aa ----
        best_coverage = 0.0
        try:
            aa_len = int(orf.get("length_aa", 0))
            uq_psites = int(float(orf.get("unique_psites", 0)))
            if aa_len > 0 and uq_psites > 0:
                # Density: P-sites per codon. Cap at 10 (very high coverage)
                best_coverage = min(1.0, (uq_psites / aa_len) / 10.0)
        except (ValueError, TypeError):
            pass

        # Compose OCS
        w_t, w_a, w_c, w_p, w_r = weights
        ocs = (
            w_t * best_translation +
            w_a * s_agreement +
            w_c * best_coverage +
            w_p * best_periodicity +
            w_r * s_readlevel
        )
        ocs = round(min(1.0, max(0.0, ocs)), 4)
        tier = _assign_tier(ocs)

        results.append({
            "orf_id": orf_id,
            "gene_id": orf.get("gene_id", ""),
            "gene_name": orf.get("gene_name", ""),
            "chrom": orf.get("chrom", ""),
            "orf_type": _harmonize_orf_type(orf.get("orf_type", "")),
            "ocs": ocs,
            "tier": tier,
            "s_translation": round(best_translation, 3),
            "s_agreement": round(s_agreement, 3),
            "s_coverage": round(best_coverage, 3),
            "s_periodicity": round(best_periodicity, 3),
            "s_readlevel": round(s_readlevel, 3),
            "detecting_tools": ",".join(sorted(detecting_tools)),
            "n_detecting": n_detecting,
        })

    # Sort by OCS descending
    results.sort(key=lambda x: x["ocs"], reverse=True)
    return results


def _read_unified_metadata(path: str) -> List[Dict]:
    """Read unified ORF metadata TSV."""
    rows = []
    if not os.path.exists(path):
        print(f"  WARNING: unified metadata not found at {path}", file=sys.stderr)
        return rows
    with open(path) as fh:
        header = None
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if header is None:
                header = fields
            else:
                rows.append({header[i]: fields[i] if i < len(fields) else "" for i in range(len(header))})
    return rows


def _parse_detecting_tools(tools_str: str, sources_str: str) -> Set[str]:
    """Parse tool names from metadata columns. Handle comma/semicolon-separated formats."""
    tools: Set[str] = set()
    for s in [tools_str, sources_str]:
        if not s or s in (".", "-", "NA", "None"):
            continue
        for part in s.replace(";", ",").split(","):
            part = part.strip().strip("'\"")
            if part:
                # Normalize tool names
                part_lower = part.lower()
                if "ribocode" in part_lower:
                    tools.add("ribocode")
                elif "ribotricer" in part_lower:
                    tools.add("ribotricer")
                elif "ribotish" in part_lower or "ribo-tish" in part_lower:
                    tools.add("ribotish")
                elif "orfquant" in part_lower:
                    tools.add("orfquant")
                elif "price" in part_lower and "rpbp" not in part_lower:
                    tools.add("price")
                elif "rpbp" in part_lower:
                    tools.add("rpbp")
                else:
                    tools.add(part)  # Keep as-is if unrecognized
    return tools


def _extract_coords(orf: Dict) -> Optional[Tuple[str, int, int, str]]:
    """Extract (chrom, start, end, strand) from an ORF dict.
    Returns None if coordinates can't be determined.
    """
    chrom = orf.get("chrom", "")
    strand = orf.get("strand", ".")

    # Try multiple coordinate field names
    start = orf.get("orf_gstart") or orf.get("start")
    end = orf.get("orf_gstop") or orf.get("end")

    # If start/end are missing, try parsing genome_pos or location
    if start is None or end is None:
        pos_str = str(orf.get("genome_pos", "") or orf.get("location", ""))
        if ":" in pos_str:
            c, s, e, st = _parse_genome_pos(pos_str)
            if c:
                return (c, s + 1, e, st)  # Convert back to 1-based

    try:
        s = int(start)
        e = int(end)
        if s > 0 and e > s:
            return (str(chrom), s, e, str(strand))
    except (ValueError, TypeError, AttributeError):
        pass
    return None


def _build_tool_spatial_index(tool_data: Dict) -> Dict[str, List[Tuple[str, int, int, str, Dict]]]:
    """Build spatial index for each tool: {tool_name: [(chrom, start, end, strand, metrics), ...]}."""
    index: Dict[str, List[Tuple[str, int, int, str, Dict]]] = defaultdict(list)
    for tool_name, samples in tool_data.get("samples", {}).items():
        for sample_data in samples.values():
            for orf in sample_data.get("orfs", []):
                coords = _extract_coords(orf)
                if coords:
                    index[tool_name].append((*coords, orf))
    return dict(index)


def _overlap_fraction(a_start: int, a_end: int, b_start: int, b_end: int) -> float:
    """Compute overlap fraction of the shorter interval."""
    overlap = max(0, min(a_end, b_end) - max(a_start, b_start))
    shorter = min(a_end - a_start, b_end - b_start)
    return overlap / shorter if shorter > 0 else 0.0


def _find_matching_orfs(
    unified_orf: Dict,
    tool_spatial_index: Dict[str, List[Tuple]],
    tool_name: str,
    min_overlap: float = 0.5,
) -> List[Dict]:
    """Find tool ORFs that overlap a unified ORF on the same chrom/strand."""
    coords = _extract_coords(unified_orf)
    if not coords:
        return []
    u_chrom, u_start, u_end, u_strand = coords

    matches = []
    for entry in tool_spatial_index.get(tool_name, []):
        t_chrom, t_start, t_end, t_strand, metrics = entry
        if t_chrom != u_chrom or t_strand != u_strand:
            continue
        if _overlap_fraction(u_start, u_end, t_start, t_end) >= min_overlap:
            matches.append(metrics)
    return matches


def _normalize_translation(metrics: Dict, tool_name: str) -> float:
    """Normalize translation significance to [0,1]."""
    if tool_name == "ribocode":
        pval = metrics.get("adjusted_pval") or metrics.get("pval_combined", 1)
        if pval <= 0:
            return 1.0
        return min(1.0, -math.log10(max(pval, 1e-300)) / 10.0)

    elif tool_name == "ribotricer":
        # phase_score directly
        return min(1.0, max(0.0, metrics.get("phase_score", 0)))

    elif tool_name == "ribotish":
        qval = metrics.get("fisher_qvalue", 1)
        if qval <= 0:
            return 1.0
        return min(1.0, -math.log10(max(qval, 1e-300)) / 10.0)

    elif tool_name == "price":
        pval = metrics.get("p_value", 1)
        if pval <= 0:
            return 1.0
        return min(1.0, -math.log10(max(pval, 1e-300)) / 10.0)

    elif tool_name == "rpbp":
        bf = metrics.get("bayes_factor", 0) or 0
        return min(1.0, math.log10(max(bf + 1, 1)) / 5.0)

    elif tool_name == "orfquant":
        # Use orf_type as proxy — annotated CDS have higher confidence
        orf_type = metrics.get("orf_type", "")
        if orf_type == "CDS":
            return 0.9
        elif orf_type in ("uORF", "dORF"):
            return 0.6
        else:
            return 0.4

    return 0.5


def _normalize_periodicity(metrics: Dict, tool_name: str) -> float:
    """Normalize ORF-level periodicity to [0,1]."""
    if tool_name == "ribocode":
        pval = metrics.get("pval_combined", 1)
        if pval <= 0:
            return 1.0
        return min(1.0, -math.log10(max(pval, 1e-300)) / 10.0)

    elif tool_name == "ribotricer":
        return min(1.0, max(0.0, metrics.get("phase_score", 0)))

    elif tool_name == "ribotish":
        qval = metrics.get("frame_qvalue", 1)
        if qval <= 0:
            return 1.0
        return min(1.0, -math.log10(max(qval, 1e-300)) / 10.0)

    elif tool_name == "rpbp":
        bf = metrics.get("bayes_factor", 0) or 0
        return min(1.0, math.log10(max(bf + 1, 1)) / 5.0)

    elif tool_name == "orfquant":
        return 0.5  # Not directly available

    elif tool_name == "price":
        return 0.5  # Not directly available

    return 0.3


def _normalize_coverage(metrics: Dict, tool_name: str) -> float:
    """Normalize coverage completeness to [0,1]."""
    if tool_name == "ribocode":
        return min(1.0, max(0.0, metrics.get("coverage_frame0", 0)))

    elif tool_name == "ribotricer":
        return min(1.0, max(0.0, metrics.get("valid_codons_ratio", 0)))

    # For other tools, impute from read count / length ratio
    read_count = metrics.get("read_count", 0)
    length = metrics.get("length_nt", metrics.get("orf_length_nt", 1))
    if length > 0 and read_count > 0:
        density = read_count / length
        return min(1.0, density / 10.0)  # >10 reads/nt = full coverage

    return 0.0


def _compute_readlevel_modifier(periodicity_data: Dict) -> float:
    """Compute global read-level modifier from periodicity assessment."""
    score = periodicity_data.get("aggregate_score", 0.5)
    if score is None:
        return 0.5
    if score >= 0.7:
        return 1.0
    elif score >= 0.5:
        return 0.8
    elif score >= 0.3:
        return 0.5
    else:
        return 0.2


def _assign_tier(ocs: float) -> str:
    """Assign confidence tier based on OCS."""
    if ocs >= 0.7:
        return "High"
    elif ocs >= 0.4:
        return "Medium"
    elif ocs >= 0.2:
        return "Low"
    else:
        return "Uncertain"


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Cross-tool ORF comparison and confidence scoring"
    )
    parser.add_argument("--tool-data", required=True, help="tool_data.json")
    parser.add_argument("--unified-bed", help="unified_orfs.bed (for overlap computation)")
    parser.add_argument("--unified-meta", required=True, help="unified_orfs.metadata.tsv")
    parser.add_argument("--periodicity", help="periodicity.json from harmonize_orf_qc.py")
    parser.add_argument("--output-prefix", "-o", default="qc", help="Output file prefix")
    parser.add_argument("--weights", type=float, nargs=5, default=[0.30, 0.30, 0.20, 0.15, 0.05],
                        help="OCS weights: translation agreement coverage periodicity readlevel")
    args = parser.parse_args()

    with open(args.tool_data) as fh:
        tool_data = json.load(fh)

    periodicity_data = {}
    if args.periodicity and os.path.exists(args.periodicity):
        with open(args.periodicity) as fh:
            periodicity_data = json.load(fh)

    # Phase 3: Cross-tool agreement
    print("Computing cross-tool agreement...")
    agreement = compute_tool_agreement(tool_data)
    with open(f"{args.output_prefix}_tool_agreement.json", "w") as fh:
        json.dump(agreement, fh, indent=2)
    if agreement.get("pairwise"):
        with open(f"{args.output_prefix}_tool_agreement.tsv", "w") as fh:
            fh.write("tool_a\ttool_b\tjaccard\toverlap_count\ta_total\tb_total\n")
            for p in agreement["pairwise"]:
                fh.write(f"{p['tool_a']}\t{p['tool_b']}\t{p['jaccard']}\t{p['overlap_count']}\t{p['a_total']}\t{p['b_total']}\n")
        print(f"  Pairwise comparisons: {len(agreement['pairwise'])}")
    else:
        print(f"  {agreement.get('status', 'UNKNOWN')}")

    # Phase 4: OCS scoring
    print("Computing ORF confidence scores...")
    weights = tuple(args.weights)
    orf_scores = compute_ocs(args.unified_meta, tool_data, periodicity_data, weights)

    with open(f"{args.output_prefix}_orf_confidence.tsv", "w") as fh:
        headers = [
            "orf_id", "gene_id", "gene_name", "chrom", "orf_type",
            "ocs", "tier", "s_translation", "s_agreement", "s_coverage",
            "s_periodicity", "s_readlevel", "detecting_tools", "n_detecting"
        ]
        fh.write("\t".join(headers) + "\n")
        for r in orf_scores:
            fh.write("\t".join(str(r.get(h, "")) for h in headers) + "\n")

    # Summary
    tiers = defaultdict(int)
    for r in orf_scores:
        tiers[r["tier"]] += 1
    print(f"  Total ORFs scored: {len(orf_scores)}")
    for tier in ["High", "Medium", "Low", "Uncertain"]:
        print(f"    {tier}: {tiers.get(tier, 0)}")

    mean_ocs = sum(r["ocs"] for r in orf_scores) / len(orf_scores) if orf_scores else 0
    print(f"  Mean OCS: {mean_ocs:.3f}")


if __name__ == "__main__":
    main()
