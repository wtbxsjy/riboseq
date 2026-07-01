#!/usr/bin/env python3
"""
extract_orf_qc_metrics.py — Phase 1 of ORF QC Module

Auto-detects which tools produced output in a results directory, parses each
tool's output format, and extracts QC metrics into a standardized JSON schema.

Output: tool_data.json

Supported tools (8 total):
  Pure QC (2):    RiboseQC, riboWaltz
  ORF Predictors (6): RiboCode, ORFquant, Ribotricer, Ribo-TISH, PRICE, rp-bp

Usage:
  extract_orf_qc_metrics.py --results-dir results/orf_predictions/ --output tool_data.json
  extract_orf_qc_metrics.py --file-list samples.txt --output tool_data.json
"""

import argparse
import json
import os
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple


# ---------------------------------------------------------------------------
# Tool registry — each tool has: name, category, glob patterns, parser function
# ---------------------------------------------------------------------------

TOOL_REGISTRY: Dict[str, Dict[str, Any]] = {}


def register_tool(name: str, category: str, patterns: List[str]):
    """Decorator to register a tool parser."""
    def decorator(func):
        TOOL_REGISTRY[name] = {
            "name": name,
            "category": category,
            "patterns": patterns,
            "parser": func,
        }
        return func
    return decorator


# ---------------------------------------------------------------------------
# Utility functions
# ---------------------------------------------------------------------------

def find_files(directory: str, patterns: List[str]) -> Dict[str, List[str]]:
    """Find files matching glob patterns in a directory tree.
    Returns {pattern: [matching_paths]}.
    """
    results: Dict[str, List[str]] = defaultdict(list)
    base = Path(directory)
    if not base.exists():
        return results
    for pattern in patterns:
        for match in base.rglob(pattern):
            results[pattern].append(str(match))
    return results


def safe_open(path: str, mode: str = "r"):
    """Open a file that may be gzipped."""
    if path.endswith(".gz"):
        import gzip
        return gzip.open(path, mode + "t")
    return open(path, mode)


def read_tsv(path: str, comment_char: str = "#") -> List[Dict[str, str]]:
    """Read a TSV file into a list of dicts, skipping comment lines."""
    rows = []
    with safe_open(path) as fh:
        header = None
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith(comment_char):
                continue
            fields = line.split("\t")
            if header is None:
                header = [h.strip() for h in fields]
            else:
                row = {header[i]: fields[i] if i < len(fields) else ""
                       for i in range(len(header))}
                rows.append(row)
    return rows


# Known tool output suffixes — stripped from filename to get clean sample_id
_TOOL_SUFFIXES = [
    "_collapsed.txt", "_collapsed.bed.gz", "_collapsed.gtf.gz",
    ".txt", ".gtf.gz", ".bed.gz",
    "_P_sites_calcs", "_P_sites_plus.bedgraph", "_P_sites_minus.bedgraph",
    "_for_ORFquant", "_results_RiboseQC", "_results_RiboseQC_all",
    "_psite_offset.tsv", "_psite_offset.txt",
    "_cds_coverage.tsv", "_codon_usage.tsv",
    "_frame_distribution.tsv", "_region_distribution.tsv",
    "_translating_ORFs.tsv", "_psite_offsets.txt",
    "_metagene_profiles_5p.tsv", "_metagene_profiles_3p.tsv",
    "_pred.txt", "_all.txt", "_transprofile.py",
    ".para.py", ".orfs.tsv",
    "_Detected_ORFs.gtf.gz", "_Detected_ORFs.gtf",
    "_Protein_sequences.fasta", "_final_ORFquant_results",
    "_tmp_ORFquant_results",
    ".bayes-factors.bed.gz", ".predicted-orfs.bed.gz",
    ".predicted-orfs.dna.fa", ".predicted-orfs.protein.fa",
    "-bayes-factors.bed.gz", "-predicted-orfs.bed.gz",
    ".periodic-offsets.csv.gz",
    "_qual.txt", "_qual.pdf",
    "_ribowaltz_plots",
]


def extract_sample_id(filepath: str) -> str:
    """Extract a clean sample ID from a tool output filename by stripping
    all known tool-specific suffixes. Returns the original filename stem
    if no known suffix matches.
    """
    basename = os.path.basename(filepath)
    # Sort by length descending so longer suffixes match first
    for suffix in sorted(_TOOL_SUFFIXES, key=len, reverse=True):
        if basename.endswith(suffix):
            return basename[:-len(suffix)]
    # Fallback: remove extension
    return os.path.splitext(basename)[0]


def parse_python_dict_file(path: str, key: str = "offdict") -> Optional[Dict[int, int]]:
    """Parse a Python file containing a dict assignment like 'offdict = {28: 12, ...}'.
    Handles nested dicts (e.g., 'm0': {...}) by extracting only integer-keyed entries.
    """
    try:
        with open(path) as fh:
            content = fh.read()
        # Extract the outer dict
        match = re.search(rf'{key}\s*=\s*(\{{.+?\}})', content, re.DOTALL)
        if match:
            raw = eval(match.group(1))
            # Filter to integer keys only (skip 'm0', etc.)
            return {int(k): int(v) for k, v in raw.items()
                    if isinstance(k, int) or (isinstance(k, str) and k.isdigit())}
    except Exception:
        pass
    return None


# ---------------------------------------------------------------------------
# Tool parsers
# ---------------------------------------------------------------------------

@register_tool("ribocode", "orf_predictor",
               ["*_collapsed.txt", "ribocode/*_collapsed.txt"])
def parse_ribocode(paths: Dict[str, List[str]]) -> Dict[str, Any]:
    """Parse RiboCode _collapsed.txt output."""
    results = {}
    all_files = paths.get("*_collapsed.txt", []) + paths.get("ribocode/*_collapsed.txt", [])
    for fpath in all_files:
        sample_id = extract_sample_id(fpath)
        rows = read_tsv(fpath)
        if not rows:
            results[sample_id] = {"status": "EMPTY", "orf_count": 0, "orfs": []}
            continue

        orfs = []
        for row in rows:
            orf = {
                "orf_id": row.get("ORF_ID", ""),
                "orf_type": row.get("ORF_type", ""),
                "transcript_id": row.get("transcript_id", ""),
                "gene_id": row.get("gene_id", ""),
                "gene_name": row.get("gene_name", ""),
                "chrom": row.get("chrom", ""),
                "strand": row.get("strand", ""),
                "orf_length_nt": int(row.get("ORF_length", 0)),
                # Frame-0/1/2 P-site sums
                "psites_frame0": int(float(row.get("Psites_sum_frame0", 0))),
                "psites_frame1": int(float(row.get("Psites_sum_frame1", 0))),
                "psites_frame2": int(float(row.get("Psites_sum_frame2", 0))),
                # Coverage
                "coverage_frame0": float(row.get("Psites_coverage_frame0", "0").rstrip("%")) / 100.0 if row.get("Psites_coverage_frame0") else 0.0,
                "coverage_frame1": float(row.get("Psites_coverage_frame1", "0").rstrip("%")) / 100.0 if row.get("Psites_coverage_frame1") else 0.0,
                "coverage_frame2": float(row.get("Psites_coverage_frame2", "0").rstrip("%")) / 100.0 if row.get("Psites_coverage_frame2") else 0.0,
                # Significance
                "pval_frame0_vs_frame1": float(row.get("pval_frame0_vs_frame1", 1)),
                "pval_frame0_vs_frame2": float(row.get("pval_frame0_vs_frame2", 1)),
                "pval_combined": float(row.get("pval_combined", 1)),
                "adjusted_pval": float(row.get("adjusted_pval", row.get("pval_combined", 1))),
                # Note: older RiboCode versions don't output adjusted_pval separately
                # Use pval_combined as fallback
                # Abundance
                "rpkm_frame0": float(row.get("Psites_frame0_RPKM", 0)),
                # Coordinates
                "orf_gstart": int(row.get("ORF_gstart", 0)) if row.get("ORF_gstart") else None,
                "orf_gstop": int(row.get("ORF_gstop", 0)) if row.get("ORF_gstop") else None,
                "orf_tstart": int(row.get("ORF_tstart", 0)) if row.get("ORF_tstart") else None,
                "orf_tstop": int(row.get("ORF_tstop", 0)) if row.get("ORF_tstop") else None,
            }
            orfs.append(orf)

        results[sample_id] = {
            "status": "OK",
            "orf_count": len(orfs),
            "orf_type_counts": _count_by(orfs, "orf_type"),
            "orfs": orfs,
        }
    return results


@register_tool("riboseqc", "pure_qc",
               ["*_P_sites_calcs", "riboseqc/*_P_sites_calcs"])
def parse_riboseqc(paths: Dict[str, List[str]]) -> Dict[str, Any]:
    """Parse RiboseQC _P_sites_calcs for P-site offset data."""
    results = {}
    all_files = paths.get("*_P_sites_calcs", []) + paths.get("riboseqc/*_P_sites_calcs", [])
    for fpath in all_files:
        sample_id = extract_sample_id(fpath)
        rows = read_tsv(fpath)

        # RiboseQC format: cutoff, frame_preference, read_length, ...
        # cutoff = P-site offset, frame_preference = f0% (0-100 scale)
        psites = []
        for row in rows:
            try:
                rl = int(float(row.get("read_length", 0)))
                if rl <= 0:
                    continue
                offset = int(float(row.get("cutoff", 0)))
                f0_pct = float(row.get("frame_preference", 0))
                # frame_preference is 0-100; normalize
                if f0_pct > 1:
                    f0_pct = f0_pct / 100.0
                psites.append({
                    "read_length": rl,
                    "p_site_offset": offset,
                    "f0_percent": f0_pct,
                })
            except (ValueError, KeyError):
                continue

        results[sample_id] = {
            "status": "OK" if psites else "EMPTY",
            "psite_count": len(psites),
            "psites": psites,
        }
    return results


@register_tool("ribowaltz", "pure_qc",
               ["*_psite_offset.tsv", "*_frame_distribution.tsv",
                "ribowaltz/*_psite_offset.tsv", "ribowaltz/*_frame_distribution.tsv"])
def parse_ribowaltz(paths: Dict[str, List[str]]) -> Dict[str, Any]:
    """Parse riboWaltz P-site offset and frame distribution outputs."""
    results: Dict[str, Dict] = defaultdict(lambda: {"status": "EMPTY", "psites": [], "frames": {}})

    # Parse P-site offsets
    offset_files = paths.get("*_psite_offset.tsv", []) + paths.get("ribowaltz/*_psite_offset.tsv", [])
    for fpath in offset_files:
        sample_id = extract_sample_id(fpath)
        rows = read_tsv(fpath)
        for row in rows:
            try:
                rl = int(float(row.get("length", 0)))
                results[sample_id]["psites"].append({
                    "read_length": rl,
                    "offset_from_5": int(float(row.get("offset_from_5", 0))),
                    "corrected_offset_from_5": int(float(row.get("corrected_offset_from_5", 0))),
                    "total_percentage": float(row.get("total_percentage", 0)),
                    "start_percentage": float(row.get("start_percentage", 0)),
                })
            except (ValueError, KeyError):
                continue
        if results[sample_id]["psites"]:
            results[sample_id]["status"] = "OK"

    # Parse frame distributions
    frame_files = paths.get("*_frame_distribution.tsv", []) + paths.get("ribowaltz/*_frame_distribution.tsv", [])
    for fpath in frame_files:
        sample_id = extract_sample_id(fpath)
        if sample_id not in results:
            results[sample_id] = {"status": "EMPTY", "psites": [], "frames": {}}
        # frame_psite() output format: complex; try to parse
        try:
            content = Path(fpath).read_text()
            results[sample_id]["frames_raw"] = content[:5000]  # First 5K chars
        except Exception:
            pass

    return dict(results)


@register_tool("ribotricer", "orf_predictor",
               ["*_translating_ORFs.tsv", "*_psite_offsets.txt",
                "ribotricer/*_translating_ORFs.tsv", "ribotricer/*_psite_offsets.txt"])
def parse_ribotricer(paths: Dict[str, List[str]]) -> Dict[str, Any]:
    """Parse Ribotricer _translating_ORFs.tsv and _psite_offsets.txt."""
    results: Dict[str, Dict] = {}

    # Parse ORF data (primary)
    orf_files = paths.get("*_translating_ORFs.tsv", []) + paths.get("ribotricer/*_translating_ORFs.tsv", [])
    for fpath in orf_files:
        sample_id = extract_sample_id(fpath)
        rows = read_tsv(fpath)
        if not rows:
            if sample_id not in results:
                results[sample_id] = {"status": "EMPTY", "orf_count": 0, "orfs": []}
            continue

        orfs = []
        for row in rows:
            # Parse genomic coordinates from ORF_ID format:
            #   {transcript_id}_{gstart}_{gstop}_{cds_length}
            # or from dedicated chrom/start/end columns if present
            orf_id = row.get("ORF_ID", "")
            gstart = row.get("start") or row.get("orf_gstart")
            gstop = row.get("end") or row.get("orf_gstop")
            if not gstart or not gstop:
                parts = orf_id.rsplit("_", 3)
                if len(parts) == 4:
                    try:
                        gstart = int(parts[1])
                        gstop = int(parts[2])
                    except (ValueError, TypeError):
                        pass
            orf = {
                "orf_id": orf_id,
                "orf_type": row.get("ORF_type", ""),
                "status": row.get("status", "nontranslating"),
                "phase_score": float(row.get("phase_score", 0)),
                "read_count": int(float(row.get("read_count", 0))),
                "length_nt": int(float(row.get("length", 0))),
                "valid_codons": int(float(row.get("valid_codons", 0))),
                "valid_codons_ratio": float(row.get("valid_codons_ratio", 0)),
                "read_density": float(row.get("read_density", 0)),
                "transcript_id": row.get("transcript_id", ""),
                "gene_id": row.get("gene_id", ""),
                "gene_name": row.get("gene_name", ""),
                "chrom": row.get("chrom", ""),
                "strand": row.get("strand", ""),
                "start_codon": row.get("start_codon", ""),
                # Genomic coordinates (for cross-tool overlap)
                "orf_gstart": gstart,
                "orf_gstop": gstop,
                # P-value from phase score
                "phase_score_pvalue": _phase_score_to_pvalue(
                    float(row.get("phase_score", 0)),
                    int(float(row.get("valid_codons", 1))),
                ) if row.get("phase_score") else None,
            }
            orfs.append(orf)

        if sample_id not in results:
            results[sample_id] = {"status": "EMPTY", "orf_count": 0, "orfs": []}
        # Only update if we have actual ORF data (don't overwrite with less data)
        if len(orfs) > 0:
            results[sample_id].update({
                "status": "OK",
                "orf_count": len(orfs),
                "translating_count": sum(1 for o in orfs if o["status"] == "translating"),
                "orf_type_counts": _count_by(orfs, "orf_type"),
                "orfs": orfs,
            })

    # Parse P-site offset data (secondary) and merge
    offset_files = paths.get("*_psite_offsets.txt", []) + paths.get("ribotricer/*_psite_offsets.txt", [])
    for fpath in offset_files:
        sid = extract_sample_id(fpath)
        if sid not in results:
            results[sid] = {"status": "EMPTY", "orf_count": 0, "orfs": []}
        try:
            with open(fpath) as fh:
                content = fh.read()
            # Parse format: "lag of LENGTH: OFFSET" lines
            offsets = {}
            for line in content.split("\n"):
                m = re.search(r'lag of (\d+):\s*(-?\d+)', line)
                if m:
                    rl = int(m.group(1))
                    lag = int(m.group(2))
                    offsets[rl] = lag
            if offsets:
                results[sid]["psite_offsets"] = offsets
        except Exception:
            pass

    return results


@register_tool("ribotish", "orf_predictor",
               ["*_pred.txt", "*.para.py", "ribotish/*_pred.txt", "ribotish/*.para.py"])
def parse_ribotish(paths: Dict[str, List[str]]) -> Dict[str, Any]:
    """Parse Ribo-TISH _pred.txt and .para.py outputs."""
    results: Dict[str, Dict] = {}

    # Parse predictions
    pred_files = paths.get("*_pred.txt", []) + paths.get("ribotish/*_pred.txt", [])
    for fpath in pred_files:
        sample_id = extract_sample_id(fpath)
        rows = read_tsv(fpath)
        if not rows:
            results[sample_id] = {"status": "EMPTY", "orf_count": 0, "orfs": []}
            continue

        orfs = []
        for row in rows:
            orf = {
                "gene_id": row.get("Gid", ""),
                "transcript_id": row.get("Tid", ""),
                "gene_name": row.get("Symbol", ""),
                "gene_type": row.get("GeneType", ""),
                "genome_pos": row.get("GenomePos", ""),
                "strand": _extract_strand(row.get("GenomePos", "")),
                "start_codon": row.get("StartCodon", ""),
                "aa_length": int(float(row.get("AALen", 0))),
                "orf_type": row.get("TisType", ""),
                "tis_pvalue": _safe_float(row.get("TISPvalue")),
                "ribo_pvalue": _safe_float(row.get("RiboPvalue")),
                "fisher_pvalue": _safe_float(row.get("FisherPvalue")),
                "frame_qvalue": _safe_float(row.get("FrameQvalue")),
                "fisher_qvalue": _safe_float(row.get("FisherQvalue")),
            }
            orfs.append(orf)

        if sample_id not in results:
            results[sample_id] = {"status": "EMPTY", "orf_count": 0, "orfs": []}
        # Only update if we have actual ORF data
        if len(orfs) > 0:
            results[sample_id].update({
                "status": "OK",
                "orf_count": len(orfs),
                "orf_type_counts": _count_by(orfs, "orf_type"),
                "orfs": orfs,
            })

    # Parse P-site offsets from .para.py (merge into existing entries)
    para_files = paths.get("*.para.py", []) + paths.get("ribotish/*.para.py", [])
    for fpath in para_files:
        sample_id = extract_sample_id(fpath)
        if sample_id not in results:
            results[sample_id] = {"status": "EMPTY", "orf_count": 0, "orfs": []}
        offdict = parse_python_dict_file(fpath, "offdict")
        if offdict:
            results[sample_id]["psite_offsets"] = offdict

    return results


@register_tool("price", "orf_predictor",
               ["*.orfs.tsv", "price/*.orfs.tsv"])
def parse_price(paths: Dict[str, List[str]]) -> Dict[str, Any]:
    """Parse PRICE (GEDI) *.orfs.tsv output."""
    results = {}
    all_files = paths.get("*.orfs.tsv", []) + paths.get("price/*.orfs.tsv", [])
    for fpath in all_files:
        sample_id = extract_sample_id(fpath)
        rows = read_tsv(fpath)
        if not rows:
            results[sample_id] = {"status": "EMPTY", "orf_count": 0, "orfs": []}
            continue

        orfs = []
        for row in rows:
            # PRICE columns: Gene, Id, Location, Candidate Location, Codon,
            #                Type, Start, Range, p value, Total
            orf = {
                "gene_id": row.get("Gene", ""),
                "orf_id": row.get("Id", ""),
                "location": row.get("Location", ""),
                "candidate_location": row.get("Candidate Location", ""),
                "start_codon": row.get("Codon", ""),
                "orf_type": row.get("Type", ""),
                "start_pos": row.get("Start", ""),
                "range": row.get("Range", ""),
                "p_value": float(row.get("p value", 1)),
                "total": float(row.get("Total", 0)) if row.get("Total") else None,
            }
            orfs.append(orf)

        results[sample_id] = {
            "status": "OK",
            "orf_count": len(orfs),
            "orf_type_counts": _count_by(orfs, "orf_type"),
            "orfs": orfs,
        }
    return results


@register_tool("rpbp", "orf_predictor",
               ["*bayes-factors.bed.gz", "*predicted-orfs.bed.gz",
                "rpbp/*bayes-factors.bed.gz", "rpbp/*predicted-orfs.bed.gz"])
def parse_rpbp(paths: Dict[str, List[str]]) -> Dict[str, Any]:
    """Parse rp-bp bayes-factors.bed.gz and predicted-orfs.bed.gz."""
    results: Dict[str, Dict] = defaultdict(lambda: {"status": "EMPTY", "orf_count": 0, "orfs": []})

    # Parse Bayes factors BED
    bf_files = paths.get("*bayes-factors.bed.gz", []) + paths.get("rpbp/*bayes-factors.bed.gz", [])
    for fpath in bf_files:
        sample_id = extract_sample_id(fpath)
        orfs = []
        with safe_open(fpath) as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) >= 6:
                    try:
                        score = float(fields[4]) if fields[4] != "." else None
                    except ValueError:
                        score = None
                    orfs.append({
                        "chrom": fields[0],
                        "start": int(fields[1]),
                        "end": int(fields[2]),
                        "name": fields[3],
                        "bayes_factor": score,
                        "strand": fields[5],
                    })
        results[sample_id]["orfs"] = orfs
        results[sample_id]["orf_count"] = len(orfs)
        if orfs:
            results[sample_id]["status"] = "OK"

    # Parse periodicity offsets
    offset_files = paths.get("*periodic-offsets.csv.gz", []) or \
                   [f for p in ["rpbp/*periodic-offsets.csv.gz", "*periodic-offsets.csv.gz"]
                    for f in paths.get(p, [])]
    # Try to find periodic offsets near BFs
    for fpath in bf_files:
        offset_f = fpath.replace("bayes-factors.bed.gz", "periodic-offsets.csv.gz")
        if os.path.exists(offset_f):
            offset_files.append(offset_f)

    for fpath in offset_files:
        sample_id = extract_sample_id(fpath)
        if sample_id not in results:
            results[sample_id] = {"status": "EMPTY", "orf_count": 0, "orfs": []}
        try:
            import pandas as pd
            df = pd.read_csv(fpath)
            results[sample_id]["psite_offsets"] = []
            for _, row in df.iterrows():
                results[sample_id]["psite_offsets"].append({
                    "read_length": int(row.get("length", 0)),
                    "p_site_offset": int(row.get("highest_peak_offset", 0)),
                    "bf_mean": float(row.get("highest_peak_bf_mean", 0)),
                    "bf_var": float(row.get("highest_peak_bf_var", 0)),
                })
        except Exception:
            pass

    return dict(results)


@register_tool("orfquant", "orf_predictor",
               ["*_Detected_ORFs.gtf.gz", "orfquant/*_Detected_ORFs.gtf.gz"])
def parse_orfquant(paths: Dict[str, List[str]]) -> Dict[str, Any]:
    """Parse ORFquant _Detected_ORFs.gtf.gz for ORF attributes.
    ORFquant stores ORFs as CDS features with ORF_id attribute.
    """
    results = {}
    all_files = paths.get("*_Detected_ORFs.gtf.gz", []) + paths.get("orfquant/*_Detected_ORFs.gtf.gz", [])
    for fpath in all_files:
        sample_id = extract_sample_id(fpath)
        orfs = {}
        with safe_open(fpath) as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 9:
                    continue
                # ORFquant stores ORFs as CDS features with ORF_id attribute
                if fields[2] not in ("CDS", "ORF"):
                    continue

                # Parse GTF attributes (space-separated key-value pairs with semicolons)
                attrs = {}
                attr_str = fields[8]
                for part in attr_str.split(";"):
                    part = part.strip()
                    if not part:
                        continue
                    # Handle format: key "value" or key value
                    m = re.match(r'(\S+)\s+"([^"]*)"', part)
                    if m:
                        attrs[m.group(1)] = m.group(2)
                    else:
                        m2 = re.match(r'(\S+)\s+(\S+)', part)
                        if m2:
                            attrs[m2.group(1)] = m2.group(2)

                orf_id = attrs.get("ORF_id", "")
                if not orf_id:
                    continue

                if orf_id not in orfs:
                    # Collect P-site metrics from ORFquant attributes
                    orfs[orf_id] = {
                        "orf_id": orf_id,
                        "gene_id": attrs.get("gene_id", ""),
                        "gene_name": attrs.get("gene_name", ""),
                        "transcript_id": attrs.get("transcript_id", ""),
                        "orf_type": attrs.get("ORF_category_Tx_compatible",
                                   attrs.get("ORF_category_Gen", "")),
                        "chrom": fields[0],
                        "start": int(fields[3]),   # GTF col 4 (1-based)
                        "end": int(fields[4]),     # GTF col 5 (1-based)
                        "strand": fields[6],
                        "p_sites": float(attrs.get("P_sites", 0)),
                        "pct_p_sites": float(attrs.get("ORF_pct_P_sites", 0)),
                        "orfs_pm": float(attrs.get("ORFs_pM", 0)),
                    }

        orf_list = list(orfs.values())
        results[sample_id] = {
            "status": "OK" if orf_list else "EMPTY",
            "orf_count": len(orf_list),
            "orf_type_counts": _count_by(orf_list, "orf_type"),
            "orfs": orf_list,
        }
    return results


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _count_by(items: List[Dict], key: str) -> Dict[str, int]:
    """Count items by a key value."""
    counts: Dict[str, int] = defaultdict(int)
    for item in items:
        counts[item.get(key, "unknown")] += 1
    return dict(counts)


def _safe_float(val: Any) -> float:
    """Convert to float, returning 1.0 for None/'None'/empty."""
    if val is None or val == "None" or val == "":
        return 1.0
    try:
        return float(val)
    except (ValueError, TypeError):
        return 1.0


def _extract_strand(genome_pos: str) -> str:
    """Extract strand from 'chr:start-end:strand' format."""
    if ":" in genome_pos:
        parts = genome_pos.split(":")
        if len(parts) >= 3 and parts[-1] in ("+", "-"):
            return parts[-1]
    return "."


def _phase_score_to_pvalue(phase_score: float, N: int) -> float:
    """Approximate p-value from phase score using non-central chi-squared.
    From ribotricer statistics.py: x = 2 * N^2 * phase_score / (N - 1).
    """
    try:
        from scipy import stats
        if phase_score <= 0 or N <= 1:
            return 1.0
        x = 2.0 * N * N * phase_score / (N - 1)
        return float(stats.ncx2.sf(x, df=2, nc=2.0 / (N - 1)))
    except ImportError:
        # Fallback: rough approximation
        return max(0.0, min(1.0, 1.0 - phase_score))


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def extract_metrics(results_dir: str, file_list: Optional[str] = None) -> Dict[str, Any]:
    """Main extraction pipeline.

    Args:
        results_dir: Root directory to scan for tool outputs.
        file_list: Optional file listing specific files to process.

    Returns:
        Nested dict: {tool_name: {sample_id: {status, metrics, ...}}}
    """
    tool_data: Dict[str, Any] = {
        "tool_status": {},   # {tool_name: "OK"|"EMPTY"|"NOT_FOUND"}
        "samples": {},       # {tool_name: {sample_id: {...}}}
    }

    # If file_list provided, build paths from it
    file_map: Dict[str, List[str]] = defaultdict(list)
    if file_list and os.path.exists(file_list):
        with open(file_list) as fh:
            for line in fh:
                path = line.strip()
                if path:
                    file_map["__all__"].append(path)
        # Classify by tool
        for path in file_map["__all__"]:
            for tool_name, tool_info in TOOL_REGISTRY.items():
                for pattern in tool_info["patterns"]:
                    # Simple match: if path contains tool name
                    if tool_name in path.lower() or pattern.replace("*", "") in path:
                        file_map[tool_name].append(path)
                        break
    else:
        # Scan directory
        for tool_name, tool_info in TOOL_REGISTRY.items():
            found = find_files(results_dir, tool_info["patterns"])
            for matches in found.values():
                file_map[tool_name].extend(matches)

    # Parse each tool
    for tool_name, tool_info in TOOL_REGISTRY.items():
        # Get files for this tool: from file_map if available, else scan directory
        tool_files = file_map.get(tool_name, [])
        if not tool_files:
            paths = find_files(results_dir, tool_info["patterns"])
        else:
            # Distribute files among patterns based on filename matching
            import fnmatch
            paths = defaultdict(list)
            for f in tool_files:
                matched = False
                basename = os.path.basename(f)
                for pattern in tool_info["patterns"]:
                    # Handle both glob patterns (e.g. "*_collapsed.txt") and
                    # path-prefixed patterns (e.g. "ribocode/*_collapsed.txt")
                    pat = pattern.split("/")[-1]  # strip directory prefix
                    if fnmatch.fnmatch(basename, pat):
                        paths[pattern].append(f)
                        matched = True
                        break
                if not matched:
                    paths[tool_info["patterns"][0]].append(f)

        if not any(paths.values()):
            tool_data["tool_status"][tool_name] = "NOT_FOUND"
            tool_data["samples"][tool_name] = {}
            continue

        try:
            parsed = tool_info["parser"](paths)
            if not parsed:
                tool_data["tool_status"][tool_name] = "EMPTY"
                tool_data["samples"][tool_name] = {}
            else:
                all_empty = all(
                    v.get("status") == "EMPTY"
                    for v in parsed.values()
                )
                tool_data["tool_status"][tool_name] = "EMPTY" if all_empty else "OK"
                tool_data["samples"][tool_name] = parsed
        except Exception as e:
            print(f"WARNING: Failed to parse {tool_name}: {e}", file=sys.stderr)
            tool_data["tool_status"][tool_name] = "FAILED"
            tool_data["samples"][tool_name] = {"error": str(e)}

    # Collect sample IDs across all tools
    all_samples: set = set()
    for samples in tool_data["samples"].values():
        all_samples.update(samples.keys())
    tool_data["sample_ids"] = sorted(all_samples)

    return tool_data


def main():
    parser = argparse.ArgumentParser(
        description="Extract ORF QC metrics from all Ribo-seq ORF prediction tools"
    )
    parser.add_argument(
        "--results-dir", "-d",
        default=".",
        help="Root results directory containing tool outputs",
    )
    parser.add_argument(
        "--file-list", "-f",
        help="File containing list of specific files to process (one per line)",
    )
    parser.add_argument(
        "--output", "-o",
        default="tool_data.json",
        help="Output JSON file path",
    )
    parser.add_argument(
        "--pretty", "-p",
        action="store_true",
        help="Pretty-print JSON output",
    )
    args = parser.parse_args()

    print(f"Scanning {args.results_dir} for ORF prediction tool outputs...")
    tool_data = extract_metrics(args.results_dir, args.file_list)

    # Summary
    print("\nTool status summary:")
    for tool, status in tool_data["tool_status"].items():
        n_samples = len(tool_data["samples"].get(tool, {}))
        print(f"  {tool:15s}  {status:10s}  ({n_samples} samples)")

    print(f"\nTotal samples detected: {len(tool_data['sample_ids'])}")

    # Write output
    indent = 2 if args.pretty else None
    with open(args.output, "w") as fh:
        json.dump(tool_data, fh, indent=indent, default=str)
    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    print(f"\nOutput written to: {args.output}")


if __name__ == "__main__":
    main()
