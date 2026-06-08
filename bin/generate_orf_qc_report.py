#!/usr/bin/env python3
"""
generate_orf_qc_report.py — Phase 5 of ORF QC Module

Generates an interactive HTML QC report using Plotly with comprehensive
methodology documentation.

Inputs:
  - psite_harmonized.tsv
  - orf_confidence.tsv
  - periodicity.json
  - tool_agreement.json (or .tsv)
  - sample_flags.json
  - tool_data.json (optional, for raw data tables)

Output: qc_report.html  +  MultiQC-compatible YAML/TSV

Usage:
  generate_orf_qc_report.py \
    --psite psite_harmonized.tsv \
    --confidence orf_confidence.tsv \
    --periodicity periodicity.json \
    --agreement tool_agreement.json \
    --flags sample_flags.json \
    --output qc_report.html
"""

import argparse
import json
import os
import sys
from typing import Any, Dict, List, Optional


def _load_tsv(path: str) -> List[Dict]:
    """Load a TSV file into list of dicts."""
    rows = []
    if not os.path.exists(path):
        return rows
    with open(path) as fh:
        header = None
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t")
            if header is None:
                header = fields
            else:
                row = {header[i]: fields[i] if i < len(fields) else "" for i in range(len(header))}
                rows.append(row)
    return rows


# ---------------------------------------------------------------------------
# Methodology markdown — embedded in the report as a self-contained reference
# ---------------------------------------------------------------------------

_METHODOLOGY_HTML = r"""
<div class="methodology">

<h2>1. Analysis Pipeline Overview</h2>
<p>
The ORF QC module runs in <strong>five sequential phases</strong>, each building
on the output of the previous one. Together they transform raw per-tool
prediction files into a single unified confidence score for every ORF.
</p>

<table class="method-table">
  <tr><th>Phase</th><th>Script</th><th>What It Does</th></tr>
  <tr>
    <td>1 &nbsp;Extract</td>
    <td><code>extract_orf_qc_metrics.py</code></td>
    <td>Scans the results directory for output files from up to 8 integrated
        ORF-prediction / QC tools (RiboCode, RiboseQC, riboWaltz, Ribotricer,
        Ribo-TISH, PRICE, rp-bp, ORFquant). For each tool it parses
        sample-level metrics: P-site offsets, periodicity (frame-0 %),
        read counts, ORF counts, and translation scores. Results are written
        to a single JSON file (<code>tool_data.json</code>).</td>
  </tr>
  <tr>
    <td>2 &nbsp;Harmonize</td>
    <td><code>harmonize_orf_qc.py</code></td>
    <td>Harmonises P-site offsets across tools using a fixed authority order
        (riboWaltz → RiboseQC → RiboCode → Ribotricer → rp-bp → Ribo-TISH →
        PRICE). The consensus offset per read length is the median of available
        values; individual offsets deviating by &gt;3 nt are flagged.
        Periodicity is normalised to a 0–1 scale and an aggregate read-level
        score is computed. Sample-level quality flags are raised.</td>
  </tr>
  <tr>
    <td>3 &nbsp;Compare</td>
    <td><code>compare_orf_tools.py</code></td>
    <td>Computes pairwise Jaccard similarity between tools by comparing the
        sets of unified ORFs each tool contributed to. A heatmap of tool
        agreement is produced.</td>
  </tr>
  <tr>
    <td>4 &nbsp;Score (OCS)</td>
    <td><code>compare_orf_tools.py</code></td>
    <td>Assigns an <strong>ORF Confidence Score (OCS)</strong> to every
        unified ORF using five orthogonal sub-scores drawn from the unified
        metadata table (see Section 2 below).</td>
  </tr>
  <tr>
    <td>5 &nbsp;Report</td>
    <td><code>generate_orf_qc_report.py</code></td>
    <td>Generates this interactive HTML report plus MultiQC-compatible
        YAML/TSV files for integration into the final pipeline report.</td>
  </tr>
</table>


<h2>2. ORF Confidence Score (OCS)</h2>
<p>
The OCS is a weighted linear combination of five sub-scores, each in the range
[0,&thinsp;1]. The weights reflect the relative importance of each dimension
for distinguishing <em>bona fide</em> translated ORFs from noise.
</p>

<div class="formula-box">
  <strong>OCS</strong> = 0.30·<em>S</em><sub>translation</sub>
  + 0.30·<em>S</em><sub>agreement</sub>
  + 0.20·<em>S</em><sub>coverage</sub>
  + 0.15·<em>S</em><sub>periodicity</sub>
  + 0.05·<em>S</em><sub>readlevel</sub>
</div>

<table class="method-table">
  <tr><th>Component</th><th>Weight</th><th>Source</th><th>Description</th></tr>
  <tr>
    <td><em>S</em><sub>translation</sub></td>
    <td>0.30</td>
    <td><code>tool_scores</code> column in unified metadata</td>
    <td>Best per-tool translation score across all tools that detected the ORF.
        For Ribotricer this is <code>phase_score</code>; for ORFquant it is
        the <code>ORF_score</code>; for Ribo-TISH it is derived from the
        <code>pvalue</code> (−log<sub>10</sub>-normalised).  A score &ge; 0.7 is
        considered strong evidence of translation.</td>
  </tr>
  <tr>
    <td><em>S</em><sub>agreement</sub></td>
    <td>0.30</td>
    <td><code>tools</code> column in unified metadata</td>
    <td><strong>Cross-tool reproducibility.</strong> Defined as the fraction of
        available ORF-prediction tools that independently called the same ORF
        locus (after coordinate-based merging). The raw fraction is discretised:
        &ge;80% → 1.0, &ge;60% → 0.75, &ge;40% → 0.50, &ge;20% → 0.25,
        &lt;20% → 0.10.  An ORF detected by 3+ tools is the strongest signal.</td>
  </tr>
  <tr>
    <td><em>S</em><sub>coverage</sub></td>
    <td>0.20</td>
    <td><code>unique_psites</code> and <code>length_aa</code> in unified metadata</td>
    <td><strong>P-site density</strong> (unique P-site positions per codon),
        normalised by a ceiling of 10 P-sites/codon.  An ORF with fewer than
        one P-site per codon on average is unlikely to be actively translated.</td>
  </tr>
  <tr>
    <td><em>S</em><sub>periodicity</sub></td>
    <td>0.15</td>
    <td><code>pN</code> column in unified metadata</td>
    <td><strong>3-nt periodicity</strong> — the fraction of P-sites that fall
        in the dominant reading frame (frame&nbsp;0).  Values near 1.0 indicate
        strong triplet periodicity, which is the hallmark of ribosome-protected
        fragments.  Values near 0.33 (random) suggest noise.</td>
  </tr>
  <tr>
    <td><em>S</em><sub>readlevel</sub></td>
    <td>0.05</td>
    <td>Aggregate periodicity score from Phase 2</td>
    <td><strong>Global read-level quality modifier.</strong> Derived from the
        aggregate periodicity classification across all read lengths:
        &ge;0.7 → 1.0, &ge;0.5 → 0.8, &ge;0.3 → 0.5, &lt;0.3 → 0.20.
        This is intentionally low weight — it modulates confidence slightly
        when overall data quality is poor but does not dominate.</td>
  </tr>
</table>

<h3>Confidence Tiers</h3>
<table class="method-table">
  <tr><th>Tier</th><th>OCS Range</th><th>Interpretation</th></tr>
  <tr><td class="tier-High">High</td><td>&ge; 0.70</td><td>Strong evidence — multiple tools agree, good translation
      score, solid periodicity and coverage. Suitable for downstream validation
      or functional studies.</td></tr>
  <tr><td class="tier-Medium">Medium</td><td>0.40 – 0.69</td><td>Moderate evidence — plausible ORF, but may lack
      cross-tool confirmation or have marginal translation metrics.</td></tr>
  <tr><td class="tier-Low">Low</td><td>0.20 – 0.39</td><td>Weak evidence — detected by a single tool, low coverage,
      or poor periodicity. Use with caution; may include false positives.</td></tr>
  <tr><td class="tier-Uncertain">Uncertain</td><td>&lt; 0.20</td><td>Insufficient evidence — typically single-tool detections
      with minimal supporting metrics. Likely noise or very low-expression
      ORFs.</td></tr>
</table>


<h2>3. P-site Offset Harmonization</h2>
<p>
Each ORF prediction tool estimates the distance from the 5′ end of the read to
the ribosomal P-site independently.  These estimates can differ by several
nucleotides between tools because of different algorithms, training data, or
read-length binning strategies.
</p>
<p>
The harmonization step (Phase 2) resolves these discrepancies:
</p>
<ol>
  <li><strong>Authority order:</strong> riboWaltz &gt; RiboseQC &gt; RiboCode
      &gt; Ribotricer &gt; rp-bp &gt; Ribo-TISH &gt; PRICE.  Earlier tools
      in this list provide more reliable or experimentally validated offsets.</li>
  <li><strong>Consensus:</strong> For each read length, the consensus offset
      is the <em>median</em> of all available tool estimates.</li>
  <li><strong>Flagging:</strong> An individual tool offset is flagged
      (<code>WARN</code>) if it deviates from the consensus by more than 3 nt.
      Read lengths where no tool provides an offset are flagged
      <code>NO_DATA</code>.</li>
</ol>


<h2>4. Cross-Tool Agreement</h2>
<p>
Pairwise Jaccard similarity between tools is computed on the final unified ORF
set.  The Jaccard index for tools A and B is:
</p>
<div class="formula-box">
  <strong>J(A,&thinsp;B)</strong> = |ORF<sub>A</sub> &cap; ORF<sub>B</sub>|
  / |ORF<sub>A</sub> &cup; ORF<sub>B</sub>|
</div>
<p>
Two ORFs are considered "the same" if they survived the frame-aware merging
step in unification — i.e., they share the same reading frame and have
substantial coordinate overlap (&ge;90% by default).  A high Jaccard between
two tools indicates they tend to call the same ORFs; a low Jaccard suggests
complementary detection (different algorithmic biases) or poor data quality.
</p>
<p>
<strong>Note:</strong> Because the unification pipeline merges overlapping ORFs
from different tools into a single representative, the final unified set
<em>under</em>-counts tool-specific ORFs that were subsumed.  The Jaccard
values should therefore be interpreted as a lower bound on true agreement.
</p>


<h2>5. Integrated Tools</h2>
<table class="method-table">
  <tr><th>Tool</th><th>Category</th><th>Key QC Metrics Provided</th></tr>
  <tr>
    <td>Ribo-TISH</td><td>ORF predictor</td>
    <td>Predicted ORF coordinates, p-value, read count, TisType classification.
        Uses a negative-binomial model on triplet-periodicity and read coverage.</td>
  </tr>
  <tr>
    <td>Ribotricer</td><td>ORF predictor</td>
    <td>Phase score (0–1, measures 3-nt periodicity), ORF coordinates, read
        count, codon validity, P-site offsets per read length.</td>
  </tr>
  <tr>
    <td>RiboCode</td><td>ORF predictor</td>
    <td>P-site counts in frames 0/1/2, combined p-value (modified
        Wilcoxon rank-sum test), ORF coordinates in transcript and genome
        space, amino-acid sequence.</td>
  </tr>
  <tr>
    <td>ORFquant</td><td>ORF predictor</td>
    <td>Splice-aware ORF detection; ORF score, read count, P-site count,
        entropy-based periodicity score.  Requires RiboseQC output.</td>
  </tr>
  <tr>
    <td>PRICE</td><td>ORF predictor</td>
    <td>GEDI/PRICE platform (EM-based); ORF type (CDS, ncRNA, Trunc, iORF,
        dORF, uORF), p-value, start-codon usage score (Start), read count.</td>
  </tr>
  <tr>
    <td>RiboseQC</td><td>Pure QC</td>
    <td>P-site offsets per read length (cutoff column), frame-0 preference (%),
        metagene profiles, P-site bedgraph tracks.</td>
  </tr>
  <tr>
    <td>riboWaltz</td><td>Pure QC</td>
    <td>P-site offset from 5′ and 3′ ends (corrected_offset_from_5),
        frame distribution per read length, CDS coverage, codon usage,
        region distribution (5′ UTR / CDS / 3′ UTR).</td>
  </tr>
  <tr>
    <td>rp-bp</td><td>ORF predictor</td>
    <td>Bayes-factor based ORF detection; Bayes factor, posterior probability,
        ORF coordinates.  (Optional; disabled by default.)</td>
  </tr>
</table>


<h2>6. Interpreting the Report</h2>
<ul>
  <li><strong>Dashboard:</strong> High-level summary — total ORFs, mean OCS,
      confidence-tier breakdown, aggregate periodicity, and any quality flags.
      Start here for a quick overview.</li>
  <li><strong>Read-Level QC:</strong> Per-read-length P-site offsets and their
      cross-tool agreement.  A wide spread in offsets (max_delta &gt; 3 nt)
      for a given read length suggests one or more tools have problematic
      offset estimates.  The periodicity bar chart shows per-tool normalised
      periodicity scores.</li>
  <li><strong>ORF Confidence:</strong> OCS histogram and ORF-type breakdown.
      A healthy experiment should show a bimodal or right-skewed OCS
      distribution (many low-score noise ORFs, a distinct population of
      high-score ORFs).  A flat or uniformly low distribution suggests
      poor data quality or insufficient tool diversity.</li>
  <li><strong>Cross-Tool:</strong> Pairwise Jaccard heatmap.  Warm colours
      indicate high overlap between tools.  Cold colours (low Jaccard) may
      indicate that tools are detecting genuinely different ORF populations
      (e.g., one tool is biased toward short ORFs, another toward long ones).</li>
  <li><strong>Detail Table:</strong> Per-ORF scores.  Use the search box to
      filter by gene name, ORF ID, or tier.  Sort by clicking column headers
      (if enabled).</li>
</ul>

</div><!-- methodology -->
"""


def _generate_html(
    sample_id: str,
    psite_data: List[Dict],
    confidence_data: List[Dict],
    periodicity: Dict,
    agreement: Dict,
    flags: Dict,
    tool_data_meta: Optional[Dict] = None,
) -> str:
    """Generate complete HTML report with embedded Plotly charts."""

    # ---- Compute statistics ----
    n_orfs = len(confidence_data)
    tier_counts = {"High": 0, "Medium": 0, "Low": 0, "Uncertain": 0}
    for r in confidence_data:
        tier_counts[r.get("tier", "Uncertain")] = tier_counts.get(r.get("tier", "Uncertain"), 0) + 1
    mean_ocs = sum(float(r.get("ocs", 0)) for r in confidence_data) / max(n_orfs, 1)
    n_periodic = sum(1 for r in psite_data if r.get("flag", "").startswith("OK"))

    # Tool agreement summary
    pairwise = agreement.get("pairwise", [])
    mean_jaccard = sum(p.get("jaccard", 0) for p in pairwise) / max(len(pairwise), 1)

    # Available tools
    tool_names = sorted(set(
        p["tool_a"] for p in pairwise
    ) | set(
        p["tool_b"] for p in pairwise
    )) if pairwise else []

    # Per-tool contributions from tool_data_meta if available
    tool_contrib = {}
    if tool_data_meta:
        for tname, tdata in tool_data_meta.get("tools", {}).items():
            if tname in ("riboseqc", "ribowaltz"):
                continue
            n = tdata.get("n_samples", 0)
            if n > 0:
                tool_contrib[tname] = n

    # ---- Quality flags HTML ----
    flags_html = ""
    flag_severities = {"CRITICAL": "#dc3545", "WARNING": "#ffc107", "INFO": "#17a2b8"}
    for f in flags.get("flags", []):
        sev = f.get("severity", "INFO")
        sev_color = flag_severities.get(sev, "#6c757d")
        flags_html += f"""
        <div style="margin:4px 0; padding:6px 10px; border-left:4px solid {sev_color}; background:#f8f9fa;">
          <strong>[{sev}]</strong> {f.get('flag', '')}: {f.get('detail', '')}
        </div>"""

    # ---- P-site table rows ----
    flag_colors = {
        "OK": "#d4edda", "OK_SINGLE": "#cce5ff",
        "WARN": "#fff3cd", "FAIL": "#f8d7da", "NO_DATA": "#e2e3e5"
    }
    psite_rows = ""
    for r in psite_data:
        fc = flag_colors.get(r.get("flag", ""), "")
        psite_rows += f"""<tr style="background:{fc}">
            <td>{r.get('read_length', '')}</td>
            <td>{r.get('n_tools', '')}</td>
            <td>{r.get('consensus_offset', '')}</td>
            <td>{r.get('max_delta', '')}</td>
            <td>{r.get('flag', '')}</td>
            <td style="font-size:12px">{r.get('tool_values', '')}</td>
        </tr>"""

    # ---- OCS histogram bins ----
    ocs_values = [float(r.get("ocs", 0)) for r in confidence_data]
    ocs_bins = [(0, 0.1), (0.1, 0.2), (0.2, 0.3), (0.3, 0.4),
                (0.4, 0.5), (0.5, 0.6), (0.6, 0.7), (0.7, 0.8),
                (0.8, 0.9), (0.9, 1.0)]
    ocs_hist_labels = []
    ocs_hist_counts = []
    for lo, hi in ocs_bins:
        label = f"{lo:.1f}-{hi:.1f}"
        cnt = sum(1 for v in ocs_values if lo <= v < hi)
        ocs_hist_labels.append(label)
        ocs_hist_counts.append(cnt)
    # Bin colours: transitional from red (0-0.2 = Uncertain) to green (0.7+ = High)
    bin_colors = ["#dc3545", "#dc3545", "#fd7e14", "#fd7e14",
                  "#ffc107", "#ffc107", "#ffc107", "#28a745",
                  "#28a745", "#28a745"]

    # ---- ORF type distribution (top 15) ----
    type_counts: Dict[str, int] = {}
    for r in confidence_data:
        t = r.get("orf_type", "Unknown")
        type_counts[t] = type_counts.get(t, 0) + 1
    top_types = sorted(type_counts.items(), key=lambda x: -x[1])[:15]

    # ---- Per-tool OCS stats ----
    tool_ocs: Dict[str, List[float]] = {}
    for r in confidence_data:
        for t in r.get("detecting_tools", "").split(","):
            t = t.strip()
            if t:
                tool_ocs.setdefault(t, []).append(float(r.get("ocs", 0)))
    tool_ocs_stats = {}
    for t, vals in tool_ocs.items():
        if vals:
            tool_ocs_stats[t] = {
                "mean": round(sum(vals) / len(vals), 3),
                "n": len(vals),
            }

    # =========================================================================
    # Build HTML
    # =========================================================================
    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>ORF QC Report — {sample_id}</title>
<script src="https://cdn.plot.ly/plotly-3.0.0.min.js"></script>
<style>
  body {{ font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
         margin: 0; padding: 20px; background: #f5f5f5; color: #333; }}
  .container {{ max-width: 1500px; margin: 0 auto; }}
  h1 {{ color: #1a1a2e; border-bottom: 3px solid #0f3460; padding-bottom: 10px; }}
  h2 {{ color: #16213e; margin-top: 30px; }}
  h3 {{ color: #0f3460; }}
  .card {{ background: white; border-radius: 8px; padding: 20px; margin: 15px 0;
          box-shadow: 0 2px 4px rgba(0,0,0,0.1); }}
  .metric-grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(180px, 1fr)); gap: 12px; }}
  .metric {{ text-align: center; padding: 15px; background: #f8f9fa; border-radius: 6px; }}
  .metric .value {{ font-size: 28px; font-weight: bold; color: #0f3460; }}
  .metric .label {{ font-size: 11px; color: #666; text-transform: uppercase; letter-spacing: 0.5px; }}
  .metric-sm {{ text-align: center; padding: 10px; background: #f0f4ff; border-radius: 6px; }}
  .metric-sm .value {{ font-size: 20px; font-weight: bold; color: #16213e; }}
  .metric-sm .label {{ font-size: 10px; color: #666; }}

  table {{ width: 100%; border-collapse: collapse; margin: 10px 0; font-size: 13px; }}
  th {{ background: #0f3460; color: white; padding: 8px 12px; text-align: left; }}
  td {{ padding: 6px 12px; border-bottom: 1px solid #eee; }}
  tr:hover {{ background: #f0f4ff; }}

  .method-table th {{ background: #16213e; font-size: 12px; }}
  .method-table td {{ font-size: 12px; vertical-align: top; }}
  .method-table code {{ background: #e8e8e8; padding: 1px 4px; border-radius: 3px; font-size: 11px; }}

  .tier-High {{ color: #28a745; font-weight: bold; }}
  .tier-Medium {{ color: #e6a817; font-weight: bold; }}
  .tier-Low {{ color: #fd7e14; font-weight: bold; }}
  .tier-Uncertain {{ color: #dc3545; font-weight: bold; }}

  .formula-box {{ background: #f0f4ff; border-left: 4px solid #0f3460; padding: 12px 16px;
                  margin: 12px 0; font-family: 'Courier New', monospace; font-size: 14px; }}

  .tab {{ overflow: hidden; border-bottom: 2px solid #0f3460; margin-bottom: 20px; }}
  .tab button {{ background: none; border: none; padding: 10px 20px; cursor: pointer;
                 font-size: 13px; font-weight: 500; color: #666; }}
  .tab button.active {{ color: #0f3460; border-bottom: 3px solid #0f3460; font-weight: 600; }}
  .tab-content {{ display: none; }}
  .tab-content.active {{ display: block; }}

  .methodology {{ font-size: 13px; line-height: 1.7; }}
  .methodology ol, .methodology ul {{ margin: 6px 0; padding-left: 20px; }}
  .methodology li {{ margin: 3px 0; }}

  #search {{ width: 100%; padding: 8px 12px; margin-bottom: 10px;
             border: 1px solid #ddd; border-radius: 4px; font-size: 13px; }}
  .note {{ font-size: 12px; color: #888; margin-top: 4px; }}

  .footer {{ text-align: center; color: #999; font-size: 11px; margin-top: 40px;
             padding-top: 20px; border-top: 1px solid #eee; }}
</style>
</head>
<body>
<div class="container">

<h1>&#x1F52C; ORF Prediction QC Report</h1>
<p style="color:#666">Sample: <strong>{sample_id}</strong> &nbsp;|&nbsp;
Generated by the ORF QC Module (riboseq pipeline) &nbsp;|&nbsp;
{n_orfs:,} unified ORFs scored &nbsp;|&nbsp;
{len(tool_names)} tools compared</p>

<!-- ====== Tabs ====== -->
<div class="tab">
  <button class="active" onclick="showTab(event, 'dashboard')">&#x1F4CA; Dashboard</button>
  <button onclick="showTab(event, 'readqc')">&#x1F52C; Read-Level QC</button>
  <button onclick="showTab(event, 'orfqc')">&#x1F9EC; ORF Confidence</button>
  <button onclick="showTab(event, 'crosstool')">&#x1F517; Cross-Tool</button>
  <button onclick="showTab(event, 'detail')">&#x1F4CB; Detail Table</button>
  <button onclick="showTab(event, 'methods')">&#x1F4D6; Methodology</button>
</div>

<!-- ================================================================= -->
<!-- DASHBOARD                                                          -->
<!-- ================================================================= -->
<div id="dashboard" class="tab-content active">

  <div class="card">
    <h2>Sample Summary</h2>
    <div class="metric-grid">
      <div class="metric"><div class="value">{n_orfs:,}</div><div class="label">Unified ORFs</div></div>
      <div class="metric"><div class="value">{mean_ocs:.3f}</div><div class="label">Mean OCS</div></div>
      <div class="metric"><div class="value">{tier_counts['High']:,}</div><div class="label">High Confidence</div></div>
      <div class="metric"><div class="value">{tier_counts['Medium']:,}</div><div class="label">Medium Confidence</div></div>
      <div class="metric"><div class="value">{periodicity.get('aggregate_score', 'N/A')}</div><div class="label">Periodicity Score</div></div>
      <div class="metric"><div class="value">{n_periodic}/{len(psite_data)}</div><div class="label">Periodic Read Lengths</div></div>
      <div class="metric"><div class="value">{mean_jaccard:.3f}</div><div class="label">Mean Tool Jaccard</div></div>
      <div class="metric"><div class="value">{len(tool_names)}</div><div class="label">Tools Compared</div></div>
    </div>
    <p class="note">
      <strong>Periodicity Score:</strong> {periodicity.get('classification', 'N/A')} &mdash;
      aggregate 3-nt periodicity across all read lengths.
      &ge;0.7 = GOOD, 0.5–0.7 = MARGINAL, &lt;0.5 = POOR.
    </p>
  </div>

  <div class="card">
    <h2>Quality Flags</h2>
    {flags_html if flags_html else '<p style="color:#28a745">&#x2705; No issues detected</p>'}
  </div>

  <div style="display:grid; grid-template-columns: 1fr 1fr; gap:15px">
    <div class="card">
      <h2>Confidence Tier Distribution</h2>
      <div id="tier_chart" style="height:350px"></div>
      <p class="note">Donut chart showing the proportion of ORFs in each confidence tier.
      A healthy experiment typically has &lt;5% High, 10–30% Medium+Low,
      and the remainder Uncertain (mostly noise / single-tool detections).</p>
    </div>
    <div class="card">
      <h2>OCS Distribution</h2>
      <div id="ocs_hist_dash" style="height:350px"></div>
      <p class="note">Right-skewed or bimodal distributions indicate a distinct
      population of high-confidence ORFs above the noise floor.</p>
    </div>
  </div>

  <div class="card">
    <h2>Per-Tool Mean OCS</h2>
    <div id="tool_ocs_chart" style="height:350px"></div>
    <p class="note">Mean OCS of ORFs detected by each tool. Tools with higher means
    tend to produce more reproducible, higher-quality predictions.
    Only ORF-prediction tools shown; pure-QC tools (RiboseQC, riboWaltz) excluded.</p>
  </div>
</div>

<!-- ================================================================= -->
<!-- READ-LEVEL QC                                                      -->
<!-- ================================================================= -->
<div id="readqc" class="tab-content">
  <div class="card">
    <h2>P-site Offset Harmonization</h2>
    <p style="font-size:13px;color:#666;">
    The <strong>consensus offset</strong> per read length is the median of all
    available tool estimates. A <strong>max delta &gt; 3 nt</strong> indicates
    disagreement between tools. The authority order for offset selection is:
    riboWaltz &gt; RiboseQC &gt; RiboCode &gt; Ribotricer &gt; rp-bp &gt;
    Ribo-TISH &gt; PRICE.
    </p>
    <div style="overflow-x:auto">
    <table>
      <tr><th>Read Length (nt)</th><th># Tools</th><th>Consensus Offset</th><th>Max &Delta;</th><th>Flag</th><th>Per-Tool Values</th></tr>
      {psite_rows if psite_rows else '<tr><td colspan="6">No P-site data available</td></tr>'}
    </table>
    </div>
  </div>

  <div class="card">
    <h2>Periodicity Scores by Tool</h2>
    <div id="periodicity_chart" style="height:350px"></div>
    <p class="note">Normalised periodicity (frame-0 preference) per tool, 0–1 scale.
    Values &lt;0.5 suggest the tool detected weak or no 3-nt periodicity in the data.
    Pure-QC tools (RiboseQC, riboWaltz) assess the raw reads; ORF-predictor tools
    assess the reads within their called ORFs only.</p>
  </div>

  <div class="card">
    <h2>Read-Length Coverage</h2>
    <div id="readlen_chart" style="height:300px"></div>
    <p class="note">Number of tools providing P-site offset data per read length.
    A gap at certain lengths may indicate missing data from key tools.</p>
  </div>
</div>

<!-- ================================================================= -->
<!-- ORF CONFIDENCE                                                     -->
<!-- ================================================================= -->
<div id="orfqc" class="tab-content">
  <div class="card">
    <h2>OCS Full Distribution</h2>
    <div id="ocs_hist" style="height:400px"></div>
    <p class="note">Density of ORFs across OCS bins (0–1). Bin colours match
    confidence tiers: <span style="color:#dc3545">Uncertain (0–0.2)</span>,
    <span style="color:#fd7e14">Low (0.2–0.4)</span>,
    <span style="color:#ffc107">Medium (0.4–0.7)</span>,
    <span style="color:#28a745">High (0.7–1.0)</span>.</p>
  </div>

  <div style="display:grid; grid-template-columns: 1fr 1fr; gap:15px">
    <div class="card">
      <h2>ORF Type Distribution (Top 15)</h2>
      <div id="type_chart" style="height:450px"></div>
    </div>
    <div class="card">
      <h2>Tier by Number of Detecting Tools</h2>
      <div id="tier_ntools_chart" style="height:450px"></div>
      <p class="note">ORFs detected by more tools tend to have higher confidence.
      Single-tool ORFs dominate the Uncertain tier.</p>
    </div>
  </div>
</div>

<!-- ================================================================= -->
<!-- CROSS-TOOL                                                         -->
<!-- ================================================================= -->
<div id="crosstool" class="tab-content">
  <div class="card">
    <h2>Pairwise Tool Agreement (Jaccard Similarity)</h2>
    <div id="agreement_chart" style="height:500px"></div>
    <p class="note">Jaccard index between tools computed on the final unified ORF
    set. Values range from 0 (no shared ORFs) to 1 (identical ORF sets).
    Cold cells (low Jaccard) may reflect genuinely complementary detection or
    differing algorithmic biases rather than poor quality.</p>
  </div>

  <div class="card">
    <h2>Tool Contributions</h2>
    <div class="metric-grid">
""" + "".join(
    f'<div class="metric-sm"><div class="value">{tool_contrib.get(t, "N/A")}</div><div class="label">{t} samples</div></div>'
    for t in sorted(tool_contrib.keys())
) + f"""
    </div>
  </div>
</div>

<!-- ================================================================= -->
<!-- DETAIL TABLE                                                       -->
<!-- ================================================================= -->
<div id="detail" class="tab-content">
  <div class="card">
    <h2>All Unified ORFs <span style="font-weight:normal;font-size:12px;color:#888">(showing first 1,000 of {n_orfs:,} total; full data in <code>orf_confidence.tsv</code>)</span></h2>
    <input type="text" id="search" placeholder="&#x1F50D; Search by ORF ID, gene name, tier, or detecting tools..."
           onkeyup="filterTable()">
    <div style="overflow-x:auto; max-height:600px; overflow-y:auto">
    <table id="orfTable">
      <thead>
        <tr>
          <th>ORF ID</th><th>Gene</th><th>Type</th><th>OCS</th><th>Tier</th>
          <th>S<sub>tr</sub></th><th>S<sub>ag</sub></th><th>S<sub>co</sub></th><th>S<sub>pe</sub></th><th>S<sub>rl</sub></th>
          <th>Detecting Tools</th>
        </tr>
      </thead>
      <tbody>
""" + "".join(
    f"""<tr>
        <td style="font-family:monospace;font-size:11px">{r.get('orf_id', '')}</td>
        <td>{r.get('gene_name', r.get('gene_id', ''))}</td>
        <td>{r.get('orf_type', '')}</td>
        <td style="font-weight:bold">{r.get('ocs', '')}</td>
        <td class="tier-{r.get('tier', '')}">{r.get('tier', '')}</td>
        <td>{r.get('s_translation', '')}</td>
        <td>{r.get('s_agreement', '')}</td>
        <td>{r.get('s_coverage', '')}</td>
        <td>{r.get('s_periodicity', '')}</td>
        <td>{r.get('s_readlevel', '')}</td>
        <td style="font-size:11px">{r.get('detecting_tools', '')}</td>
    </tr>""" for r in confidence_data[:1000]
) + """
      </tbody>
    </table>
    </div>
  </div>
</div>

<!-- ================================================================= -->
<!-- METHODOLOGY                                                        -->
<!-- ================================================================= -->
<div id="methods" class="tab-content">
  <div class="card">
""" + _METHODOLOGY_HTML + """
  </div>
</div>

<div class="footer">
  ORF QC Report &mdash; riboseq pipeline &mdash; nf-core/riboseq<br>
  Generated by <code>generate_orf_qc_report.py</code> (Phase 5)
</div>

</div><!-- container -->

<!-- ================================================================= -->
<!-- JAVASCRIPT                                                         -->
<!-- ================================================================= -->
<script>
function showTab(evt, tabName) {
  document.querySelectorAll('.tab-content').forEach(el => el.classList.remove('active'));
  document.querySelectorAll('.tab button').forEach(el => el.classList.remove('active'));
  document.getElementById(tabName).classList.add('active');
  evt.currentTarget.classList.add('active');
  // Trigger relayout for charts in newly shown tab
  setTimeout(function() { Plotly.relayout(tabName + ' .js-plotly-plot', {}); }, 100);
}

function filterTable() {
  const input = document.getElementById('search');
  const filter = input.value.toUpperCase();
  const table = document.getElementById('orfTable');
  const rows = table.getElementsByTagName('tr');
  for (let i = 1; i < rows.length; i++) {
    const text = rows[i].textContent || rows[i].innerText;
    rows[i].style.display = text.toUpperCase().includes(filter) ? '' : 'none';
  }
}

// ---- Tier distribution donut chart ----
Plotly.newPlot('tier_chart', [{
  type: 'pie',
  values: """ + json.dumps([tier_counts[t] for t in ["High", "Medium", "Low", "Uncertain"]]) + """,
  labels: ['High (≥0.7)', 'Medium (0.4–0.7)', 'Low (0.2–0.4)', 'Uncertain (<0.2)'],
  marker: {colors: ['#28a745', '#ffc107', '#fd7e14', '#dc3545']},
  textinfo: 'label+percent',
  hole: 0.45,
  sort: false
}], {margin: {t: 10}, showlegend: false});

// ---- Dashboard OCS histogram (compact) ----
Plotly.newPlot('ocs_hist_dash', [{
  type: 'bar',
  x: """ + json.dumps(ocs_hist_labels) + """,
  y: """ + json.dumps(ocs_hist_counts) + """,
  marker: {color: """ + json.dumps(bin_colors) + """}
}], {margin: {t: 10, b: 50}, xaxis: {title: 'OCS Range'}, yaxis: {title: 'Count'}});

// ---- Per-tool mean OCS ----
Plotly.newPlot('tool_ocs_chart', [{
  type: 'bar',
  x: """ + json.dumps(list(tool_ocs_stats.keys())) + """,
  y: """ + json.dumps([s["mean"] for s in tool_ocs_stats.values()]) + """,
  text: """ + json.dumps([f"n={s['n']:,}" for s in tool_ocs_stats.values()]) + """,
  textposition: 'auto',
  marker: {color: '#0f3460'},
  textfont: {size: 11}
}], {margin: {t: 10}, yaxis: {title: 'Mean OCS', range: [0, 1]}});

// ---- Periodicity chart ----
Plotly.newPlot('periodicity_chart', [{
  type: 'bar',
  x: """ + json.dumps(list(periodicity.get("per_tool", {}).keys())) + """,
  y: """ + json.dumps(list(periodicity.get("per_tool", {}).values())) + """,
  marker: {color: '#0f3460'},
  text: """ + json.dumps([f"{v:.3f}" for v in periodicity.get("per_tool", {}).values()]) + """,
  textposition: 'auto',
  textfont: {size: 11}
}], {margin: {t: 10}, yaxis: {title: 'Normalised Periodicity', range: [0, 1.1]}});

// ---- Read length coverage ----
const psiteRows = """ + json.dumps([{"rl": r.get("read_length", ""), "n": int(r.get("n_tools", 0))} for r in psite_data]) + """;
if (psiteRows.length > 0) {
  Plotly.newPlot('readlen_chart', [{
    type: 'bar',
    x: psiteRows.map(r => r.rl),
    y: psiteRows.map(r => r.n),
    marker: {color: '#16213e'}
  }], {margin: {t: 10}, xaxis: {title: 'Read Length (nt)'}, yaxis: {title: '# Tools with Data', dtick: 1}});
} else {
  document.getElementById('readlen_chart').innerHTML =
    '<p style="text-align:center;color:#999;padding:30px">No read-length data available</p>';
}

// ---- OCS full histogram ----
Plotly.newPlot('ocs_hist', [{
  type: 'bar',
  x: """ + json.dumps(ocs_hist_labels) + """,
  y: """ + json.dumps(ocs_hist_counts) + """,
  marker: {color: """ + json.dumps(bin_colors) + """},
  hovertemplate: '%{x}: %{y:,} ORFs<extra></extra>'
}], {margin: {t: 10, b: 60}, xaxis: {title: 'OCS Range'}, yaxis: {title: 'Number of ORFs'}});

// ---- ORF type distribution ----
Plotly.newPlot('type_chart', [{
  type: 'bar',
  orientation: 'h',
  y: """ + json.dumps([t[0] for t in reversed(top_types)]) + """,
  x: """ + json.dumps([t[1] for t in reversed(top_types)]) + """,
  marker: {color: '#16213e'}
}], {margin: {t: 10, l: 180}, xaxis: {title: 'Count'}});

// ---- Tier by number of detecting tools ----
const tierNtools = {};
""" + json.dumps([{
    "tier": r.get("tier", "Uncertain"),
    "n": int(r.get("n_detecting", 0))
} for r in confidence_data]) + """.forEach(function(r) {
  const key = r.n + ' tool' + (r.n !== 1 ? 's' : '');
  if (!tierNtools[key]) tierNtools[key] = {High: 0, Medium: 0, Low: 0, Uncertain: 0};
  tierNtools[key][r.tier] = (tierNtools[key][r.tier] || 0) + 1;
});
const ntLabels = Object.keys(tierNtools).sort((a,b) => parseInt(a) - parseInt(b));
Plotly.newPlot('tier_ntools_chart', [
  {type: 'bar', name: 'High',       x: ntLabels, y: ntLabels.map(k => tierNtools[k].High || 0),       marker: {color: '#28a745'}},
  {type: 'bar', name: 'Medium',     x: ntLabels, y: ntLabels.map(k => tierNtools[k].Medium || 0),     marker: {color: '#ffc107'}},
  {type: 'bar', name: 'Low',        x: ntLabels, y: ntLabels.map(k => tierNtools[k].Low || 0),        marker: {color: '#fd7e14'}},
  {type: 'bar', name: 'Uncertain',  x: ntLabels, y: ntLabels.map(k => tierNtools[k].Uncertain || 0),  marker: {color: '#dc3545'}}
], {margin: {t: 10}, barmode: 'stack', xaxis: {title: 'Number of Detecting Tools'}, yaxis: {title: 'ORF Count'}});

// ---- Tool agreement heatmap ----
const agreementData = """ + json.dumps(pairwise) + """;
if (agreementData.length > 0) {
  const tools = [...new Set(agreementData.flatMap(d => [d.tool_a, d.tool_b]))].sort();
  const z = tools.map(ta => tools.map(tb => {
    if (ta === tb) return 1.0;
    const match = agreementData.find(d =>
      (d.tool_a === ta && d.tool_b === tb) || (d.tool_a === tb && d.tool_b === ta));
    return match ? match.jaccard : 0;
  }));
  Plotly.newPlot('agreement_chart', [{
    type: 'heatmap',
    z: z, x: tools, y: tools,
    colorscale: [
      [0.0, '#f0f0f0'], [0.2, '#c6dbef'], [0.5, '#6baed6'],
      [0.7, '#2171b5'], [1.0, '#08306b']
    ],
    zmin: 0, zmax: 1,
    text: z.map(row => row.map(v => v.toFixed(3))),
    texttemplate: '%{text}',
    hovertemplate: '%{y} × %{x}: J=%{z:.3f}<extra></extra>'
  }], {margin: {t: 10, b: 80}, xaxis: {tickangle: 45}, yaxis: {}});
} else {
  document.getElementById('agreement_chart').innerHTML =
    '<p style="text-align:center;color:#999;padding:50px">Insufficient tool data for cross-tool comparison</p>';
}
</script>

</body>
</html>"""

    return html


def _generate_multiqc_data(sample_id: str, confidence_data: List[Dict],
                           periodicity: Dict, flags: Dict) -> tuple:
    """Generate MultiQC-compatible YAML and TSV data."""
    n_orfs = len(confidence_data)
    tier_counts = {"High": 0, "Medium": 0, "Low": 0, "Uncertain": 0}
    for r in confidence_data:
        tier_counts[r.get("tier", "Uncertain")] = tier_counts.get(r.get("tier", "Uncertain"), 0) + 1
    mean_ocs = sum(float(r.get("ocs", 0)) for r in confidence_data) / max(n_orfs, 1)

    agg_score = periodicity.get("aggregate_score", "N/A")

    yaml = f"""id: 'orf_qc'
section_name: 'ORF Prediction QC'
description: 'Unified quality control across all ORF prediction tools — OCS scoring, P-site harmonization, cross-tool agreement'
plot_type: 'table'
pconfig:
    id: 'orf_qc_table'
    title: 'ORF Prediction QC Summary'
data:
    sample_name: '{sample_id}'
    unified_orfs: {n_orfs}
    mean_ocs: {mean_ocs:.3f}
    high_confidence: {tier_counts['High']}
    medium_confidence: {tier_counts['Medium']}
    low_confidence: {tier_counts['Low']}
    uncertain_confidence: {tier_counts['Uncertain']}
    periodicity_score: {agg_score}
    periodicity_class: '{periodicity.get("classification", "N/A")}'
    flags_count: {len(flags.get('flags', []))}
"""

    tsv_headers = ["sample_name", "unified_orfs", "mean_ocs",
                   "high_confidence", "medium_confidence", "low_confidence",
                   "uncertain_confidence", "periodicity_score",
                   "periodicity_class", "flags_count"]
    tsv = "\t".join(tsv_headers) + "\n"
    tsv += "\t".join(str(x) for x in [
        sample_id, n_orfs, f"{mean_ocs:.3f}",
        tier_counts['High'], tier_counts['Medium'], tier_counts['Low'],
        tier_counts['Uncertain'], agg_score,
        periodicity.get('classification', 'N/A'),
        len(flags.get('flags', []))
    ]) + "\n"

    return yaml, tsv


def main():
    parser = argparse.ArgumentParser(description="Generate ORF QC HTML report with methodology documentation")
    parser.add_argument("--psite", help="psite_harmonized.tsv")
    parser.add_argument("--confidence", required=True, help="orf_confidence.tsv")
    parser.add_argument("--periodicity", help="periodicity.json")
    parser.add_argument("--agreement", help="tool_agreement.json")
    parser.add_argument("--flags", help="sample_flags.json")
    parser.add_argument("--tool-data", help="tool_data.json (optional, for richer tool stats)")
    parser.add_argument("--sample-id", default="sample", help="Sample identifier for report title")
    parser.add_argument("--output", "-o", default="qc_report.html", help="Output HTML file")
    parser.add_argument("--mqc-prefix", default="joint_riboseq_qc", help="MultiQC output prefix")
    args = parser.parse_args()

    # Load data
    psite_data = _load_tsv(args.psite) if args.psite else []
    confidence_data = _load_tsv(args.confidence)

    periodicity = {}
    if args.periodicity and os.path.exists(args.periodicity):
        with open(args.periodicity) as fh:
            periodicity = json.load(fh)

    agreement = {}
    if args.agreement and os.path.exists(args.agreement):
        with open(args.agreement) as fh:
            agreement = json.load(fh)

    flags = {}
    if args.flags and os.path.exists(args.flags):
        with open(args.flags) as fh:
            flags = json.load(fh)

    tool_data_meta = None
    if args.tool_data and os.path.exists(args.tool_data):
        with open(args.tool_data) as fh:
            tool_data_meta = json.load(fh)

    if not confidence_data:
        print("ERROR: No ORF confidence data found — cannot generate report", file=sys.stderr)
        sys.exit(1)

    html = _generate_html(args.sample_id, psite_data, confidence_data,
                          periodicity, agreement, flags, tool_data_meta)

    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    with open(args.output, "w") as fh:
        fh.write(html)

    # Generate MultiQC-compatible output (replaces JOINT_QC_REPORT)
    mqc_yaml, mqc_tsv = _generate_multiqc_data(args.sample_id, confidence_data, periodicity, flags)
    with open(f"{args.mqc_prefix}_mqc.yaml", "w") as fh:
        fh.write(mqc_yaml)
    with open(f"{args.mqc_prefix}_mqc.txt", "w") as fh:
        fh.write(mqc_tsv)

    print(f"Report written to: {args.output}")
    print(f"  {len(confidence_data):,} ORFs scored")
    print(f"  {len(psite_data)} read lengths with P-site data")
    print(f"  {len(agreement.get('pairwise', []))} pairwise tool comparisons")
    print(f"  {len(flags.get('flags', []))} quality flags")
    print(f"  MultiQC output: {args.mqc_prefix}_mqc.yaml, {args.mqc_prefix}_mqc.txt")


if __name__ == "__main__":
    main()
