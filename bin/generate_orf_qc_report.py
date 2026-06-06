#!/usr/bin/env python3
"""
generate_orf_qc_report.py — Phase 5 of ORF QC Module

Generates an interactive HTML QC report using Plotly.

Inputs:
  - psite_harmonized.tsv
  - orf_confidence.tsv
  - periodicity.json
  - tool_agreement.json (or .tsv)
  - sample_flags.json
  - tool_data.json (optional, for raw data tables)

Output: qc_report.html

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


def _generate_html(
    sample_id: str,
    psite_data: List[Dict],
    confidence_data: List[Dict],
    periodicity: Dict,
    agreement: Dict,
    flags: Dict,
) -> str:
    """Generate complete HTML report with embedded Plotly charts."""

    # Count statistics
    n_orfs = len(confidence_data)
    tier_counts = {"High": 0, "Medium": 0, "Low": 0, "Uncertain": 0}
    for r in confidence_data:
        tier_counts[r.get("tier", "Uncertain")] = tier_counts.get(r.get("tier", "Uncertain"), 0) + 1
    mean_ocs = sum(float(r.get("ocs", 0)) for r in confidence_data) / max(n_orfs, 1)
    n_periodic = sum(1 for r in psite_data if r.get("flag", "").startswith("OK"))

    # Tool agreement summary
    pairwise = agreement.get("pairwise", [])
    mean_jaccard = sum(p.get("jaccard", 0) for p in pairwise) / max(len(pairwise), 1)

    flags_html = ""
    for f in flags.get("flags", []):
        sev_color = {"CRITICAL": "#dc3545", "WARNING": "#ffc107", "INFO": "#17a2b8"}.get(
            f.get("severity", "INFO"), "#6c757d"
        )
        flags_html += f"""
        <div style="margin:4px 0; padding:6px 10px; border-left:4px solid {sev_color}; background:#f8f9fa;">
          <strong>[{f.get('severity', 'INFO')}]</strong> {f.get('flag', '')}: {f.get('detail', '')}
        </div>"""

    # P-site table
    psite_rows = ""
    for r in psite_data:
        flag_color = {"OK": "#d4edda", "WARN": "#fff3cd", "FAIL": "#f8d7da", "OK_SINGLE": "#cce5ff", "NO_DATA": "#e2e3e5"}.get(r.get("flag", ""), "")
        psite_rows += f"""<tr style="background:{flag_color}">
            <td>{r.get('read_length', '')}</td>
            <td>{r.get('n_tools', '')}</td>
            <td>{r.get('consensus_offset', '')}</td>
            <td>{r.get('max_delta', '')}</td>
            <td>{r.get('flag', '')}</td>
        </tr>"""

    # Confidence tier distribution data
    tier_json = json.dumps([
        {"tier": t, "count": tier_counts[t]} for t in ["High", "Medium", "Low", "Uncertain"]
    ])

    # ORF type distribution (if available)
    type_counts = {}
    for r in confidence_data:
        t = r.get("orf_type", "Unknown")
        type_counts[t] = type_counts.get(t, 0) + 1
    type_json = json.dumps([{"type": k, "count": v} for k, v in sorted(type_counts.items(), key=lambda x: -x[1])[:10]])

    # OCS histogram data
    ocs_values = [float(r.get("ocs", 0)) for r in confidence_data]
    ocs_bins = [0, 0.2, 0.4, 0.7, 1.0]
    ocs_hist = {}
    for v in ocs_values:
        for i in range(len(ocs_bins) - 1):
            if ocs_bins[i] <= v < ocs_bins[i + 1]:
                label = f"{ocs_bins[i]:.1f}-{ocs_bins[i+1]:.1f}"
                ocs_hist[label] = ocs_hist.get(label, 0) + 1
                break
    ocs_hist_json = json.dumps([{"range": k, "count": v} for k, v in ocs_hist.items()])

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
  .container {{ max-width: 1400px; margin: 0 auto; }}
  h1 {{ color: #1a1a2e; border-bottom: 3px solid #0f3460; padding-bottom: 10px; }}
  h2 {{ color: #16213e; margin-top: 30px; }}
  .card {{ background: white; border-radius: 8px; padding: 20px; margin: 15px 0;
          box-shadow: 0 2px 4px rgba(0,0,0,0.1); }}
  .metric-grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 15px; }}
  .metric {{ text-align: center; padding: 15px; background: #f8f9fa; border-radius: 6px; }}
  .metric .value {{ font-size: 32px; font-weight: bold; color: #0f3460; }}
  .metric .label {{ font-size: 12px; color: #666; text-transform: uppercase; }}
  table {{ width: 100%; border-collapse: collapse; margin: 10px 0; }}
  th {{ background: #0f3460; color: white; padding: 8px 12px; text-align: left; }}
  td {{ padding: 6px 12px; border-bottom: 1px solid #eee; }}
  tr:hover {{ background: #f8f9fa; }}
  .tier-High {{ color: #28a745; font-weight: bold; }}
  .tier-Medium {{ color: #ffc107; font-weight: bold; }}
  .tier-Low {{ color: #fd7e14; font-weight: bold; }}
  .tier-Uncertain {{ color: #dc3545; font-weight: bold; }}
  .tab {{ overflow: hidden; border-bottom: 2px solid #0f3460; margin-bottom: 20px; }}
  .tab button {{ background: none; border: none; padding: 10px 20px; cursor: pointer;
                 font-size: 14px; font-weight: 500; color: #666; }}
  .tab button.active {{ color: #0f3460; border-bottom: 3px solid #0f3460; }}
  .tab-content {{ display: none; }}
  .tab-content.active {{ display: block; }}
</style>
</head>
<body>
<div class="container">

<h1>🔬 ORF Prediction QC Report</h1>
<p style="color:#666">Sample: <strong>{sample_id}</strong> | Generated by riboseq pipeline ORF QC Module</p>

<!-- Tabs -->
<div class="tab">
  <button class="active" onclick="showTab(event, 'dashboard')">📊 Dashboard</button>
  <button onclick="showTab(event, 'readqc')">🔬 Read-Level QC</button>
  <button onclick="showTab(event, 'orfqc')">🧬 ORF Confidence</button>
  <button onclick="showTab(event, 'crosstool')">🔗 Cross-Tool</button>
  <button onclick="showTab(event, 'detail')">📋 Detail Table</button>
</div>

<!-- Dashboard -->
<div id="dashboard" class="tab-content active">
  <div class="card">
    <h2>Sample Summary</h2>
    <div class="metric-grid">
      <div class="metric"><div class="value">{n_orfs}</div><div class="label">Unified ORFs</div></div>
      <div class="metric"><div class="value">{mean_ocs:.2f}</div><div class="label">Mean OCS</div></div>
      <div class="metric"><div class="value">{tier_counts['High']}</div><div class="label">High Confidence</div></div>
      <div class="metric"><div class="value">{periodicity.get('aggregate_score', 'N/A')}</div><div class="label">Periodicity Score</div></div>
      <div class="metric"><div class="value">{n_periodic}/{len(psite_data)}</div><div class="label">Periodic Read Lengths</div></div>
      <div class="metric"><div class="value">{mean_jaccard:.2f}</div><div class="label">Mean Jaccard</div></div>
    </div>
  </div>

  <div class="card">
    <h2>Quality Flags</h2>
    {flags_html if flags_html else '<p style="color:#28a745">✅ No issues detected</p>'}
  </div>

  <div class="card">
    <h2>Confidence Tier Distribution</h2>
    <div id="tier_chart" style="height:350px"></div>
  </div>
</div>

<!-- Read-Level QC -->
<div id="readqc" class="tab-content">
  <div class="card">
    <h2>P-site Offset Harmonization</h2>
    <div style="overflow-x:auto">
    <table>
      <tr><th>Read Length</th><th>N Tools</th><th>Consensus Offset</th><th>Max Delta</th><th>Flag</th></tr>
      {psite_rows if psite_rows else '<tr><td colspan="5">No P-site data available</td></tr>'}
    </table>
    </div>
  </div>

  <div class="card">
    <h2>Periodicity Scores by Tool</h2>
    <div id="periodicity_chart" style="height:300px"></div>
  </div>
</div>

<!-- ORF Confidence -->
<div id="orfqc" class="tab-content">
  <div class="card">
    <h2>OCS Distribution</h2>
    <div id="ocs_hist" style="height:350px"></div>
  </div>
  <div class="card">
    <h2>ORF Type Distribution</h2>
    <div id="type_chart" style="height:400px"></div>
  </div>
</div>

<!-- Cross-Tool -->
<div id="crosstool" class="tab-content">
  <div class="card">
    <h2>Tool Agreement</h2>
    <div id="agreement_chart" style="height:500px"></div>
  </div>
</div>

<!-- Detail Table -->
<div id="detail" class="tab-content">
  <div class="card">
    <h2>All Unified ORFs</h2>
    <input type="text" id="search" placeholder="🔍 Search ORFs..." style="width:100%;padding:8px;margin-bottom:10px;border:1px solid #ddd;border-radius:4px"
           onkeyup="filterTable()">
    <div style="overflow-x:auto; max-height:600px; overflow-y:auto">
    <table id="orfTable">
      <thead>
        <tr><th>ORF ID</th><th>Gene</th><th>Type</th><th>OCS</th><th>Tier</th><th>Translation</th><th>Agreement</th><th>Coverage</th><th>Periodicity</th><th>Detecting Tools</th></tr>
      </thead>
      <tbody>
""" + "".join(
    f"""<tr>
        <td>{r.get('orf_id', '')}</td>
        <td>{r.get('gene_name', r.get('gene_id', ''))}</td>
        <td>{r.get('orf_type', '')}</td>
        <td>{r.get('ocs', '')}</td>
        <td class="tier-{r.get('tier', '')}">{r.get('tier', '')}</td>
        <td>{r.get('s_translation', '')}</td>
        <td>{r.get('s_agreement', '')}</td>
        <td>{r.get('s_coverage', '')}</td>
        <td>{r.get('s_periodicity', '')}</td>
        <td>{r.get('detecting_tools', '')}</td>
    </tr>""" for r in confidence_data[:500]  # Limit to first 500 for performance
) + """
      </tbody>
    </table>
    </div>
  </div>
</div>

</div><!-- container -->

<script>
function showTab(evt, tabName) {
  document.querySelectorAll('.tab-content').forEach(el => el.classList.remove('active'));
  document.querySelectorAll('.tab button').forEach(el => el.classList.remove('active'));
  document.getElementById(tabName).classList.add('active');
  evt.currentTarget.classList.add('active');
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

// Tier distribution chart
Plotly.newPlot('tier_chart', [{
  type: 'pie',
  values: """ + json.dumps([tier_counts[t] for t in ["High", "Medium", "Low", "Uncertain"]]) + """,
  labels: ['High', 'Medium', 'Low', 'Uncertain'],
  marker: {colors: ['#28a745', '#ffc107', '#fd7e14', '#dc3545']},
  textinfo: 'label+percent',
  hole: 0.4
}], {margin: {t: 0}});

// Periodicity chart
Plotly.newPlot('periodicity_chart', [{
  type: 'bar',
  x: """ + json.dumps(list(periodicity.get("per_tool", {}).keys())) + """,
  y: """ + json.dumps(list(periodicity.get("per_tool", {}).values())) + """,
  marker: {color: '#0f3460'}
}], {margin: {t: 0}, yaxis: {title: 'Normalized Score', range: [0, 1]}});

// OCS histogram
Plotly.newPlot('ocs_hist', [{
  type: 'bar',
  x: """ + json.dumps(list(ocs_hist.keys())) + """,
  y: """ + json.dumps(list(ocs_hist.values())) + """,
  marker: {color: ['#dc3545', '#fd7e14', '#ffc107', '#28a745']}
}], {margin: {t: 0}, xaxis: {title: 'OCS Range'}, yaxis: {title: 'Count'}});

// ORF type distribution
Plotly.newPlot('type_chart', [{
  type: 'bar',
  x: """ + json.dumps([t["type"] for t in json.loads(type_json)]) + """,
  y: """ + json.dumps([t["count"] for t in json.loads(type_json)]) + """,
  marker: {color: '#16213e'}
}], {margin: {t: 0, b: 100}, yaxis: {title: 'Count'}});

// Tool agreement heatmap
const agreementData = """ + json.dumps(pairwise) + """;
if (agreementData.length > 0) {
  const tools = [...new Set(agreementData.flatMap(d => [d.tool_a, d.tool_b]))];
  const z = tools.map(ta => tools.map(tb => {
    if (ta === tb) return 1.0;
    const match = agreementData.find(d =>
      (d.tool_a === ta && d.tool_b === tb) || (d.tool_a === tb && d.tool_b === ta));
    return match ? match.jaccard : 0;
  }));
  Plotly.newPlot('agreement_chart', [{
    type: 'heatmap',
    z: z, x: tools, y: tools,
    colorscale: 'Blues',
    zmin: 0, zmax: 1,
    text: z.map(row => row.map(v => v.toFixed(3))),
    texttemplate: '%{text}'
  }], {margin: {t: 0}, xaxis: {title: ''}, yaxis: {title: ''}});
} else {
  document.getElementById('agreement_chart').innerHTML =
    '<p style="text-align:center;color:#999;padding:50px">Insufficient tool data for comparison</p>';
}
</script>

</body>
</html>"""

    return html


def _generate_multiqc_data(sample_id: str, confidence_data: List[Dict],
                           periodicity: Dict, flags: Dict) -> tuple:
    """Generate MultiQC-compatible YAML and TSV data.

    Returns (yaml_str, tsv_str).
    """
    n_orfs = len(confidence_data)
    tier_counts = {"High": 0, "Medium": 0, "Low": 0, "Uncertain": 0}
    for r in confidence_data:
        tier_counts[r.get("tier", "Uncertain")] = tier_counts.get(r.get("tier", "Uncertain"), 0) + 1
    mean_ocs = sum(float(r.get("ocs", 0)) for r in confidence_data) / max(n_orfs, 1)

    agg_score = periodicity.get("aggregate_score", "N/A")
    psite_ok = sum(1 for f in flags.get("summary", {}).get("flags", []) if f == "OK")

    yaml = f"""id: 'orf_qc'
section_name: 'ORF Prediction QC'
description: 'Unified quality control across all ORF prediction tools'
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
    periodicity_score: {agg_score}
    periodicity_class: '{periodicity.get("classification", "N/A")}'
    flags_count: {len(flags.get('flags', []))}
"""

    tsv = f"sample_name\tunified_orfs\tmean_ocs\thigh_confidence\tmedium_confidence\tlow_confidence\tperiodicity_score\tperiodicity_class\tflags_count\n"
    tsv += f"{sample_id}\t{n_orfs}\t{mean_ocs:.3f}\t{tier_counts['High']}\t{tier_counts['Medium']}\t{tier_counts['Low']}\t{agg_score}\t{periodicity.get('classification', 'N/A')}\t{len(flags.get('flags', []))}\n"

    return yaml, tsv


def main():
    parser = argparse.ArgumentParser(description="Generate ORF QC HTML report")
    parser.add_argument("--psite", help="psite_harmonized.tsv")
    parser.add_argument("--confidence", required=True, help="orf_confidence.tsv")
    parser.add_argument("--periodicity", help="periodicity.json")
    parser.add_argument("--agreement", help="tool_agreement.json")
    parser.add_argument("--flags", help="sample_flags.json")
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

    if not confidence_data:
        print("ERROR: No ORF confidence data found — cannot generate report", file=sys.stderr)
        sys.exit(1)

    html = _generate_html(args.sample_id, psite_data, confidence_data, periodicity, agreement, flags)

    with open(args.output, "w") as fh:
        fh.write(html)

    # Generate MultiQC-compatible output (replaces JOINT_QC_REPORT)
    mqc_yaml, mqc_tsv = _generate_multiqc_data(args.sample_id, confidence_data, periodicity, flags)
    with open(f"{args.mqc_prefix}_mqc.yaml", "w") as fh:
        fh.write(mqc_yaml)
    with open(f"{args.mqc_prefix}_mqc.txt", "w") as fh:
        fh.write(mqc_tsv)

    print(f"Report written to: {args.output}")
    print(f"  {len(confidence_data)} ORFs scored")
    print(f"  {len(psite_data)} read lengths with P-site data")
    print(f"  MultiQC output: {args.mqc_prefix}_mqc.yaml, {args.mqc_prefix}_mqc.txt")


if __name__ == "__main__":
    main()
