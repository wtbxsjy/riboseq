from __future__ import annotations

import html
import re
import shlex
import subprocess
import time
from collections import Counter
from pathlib import Path
from typing import Callable

try:
    from IPython.display import HTML, display, clear_output
    MISSING_NOTEBOOK_DEPS: list[str] = []
except ModuleNotFoundError as exc:  # pragma: no cover - import-time dependency guard
    HTML = None  # type: ignore[assignment]
    display = None  # type: ignore[assignment]
    clear_output = None  # type: ignore[assignment]
    MISSING_NOTEBOOK_DEPS = [exc.name or "IPython"]

try:
    import ipywidgets as widgets
    import matplotlib.pyplot as plt
    import pandas as pd
except ModuleNotFoundError as exc:  # pragma: no cover - import-time dependency guard
    widgets = None  # type: ignore[assignment]
    plt = None  # type: ignore[assignment]
    pd = None  # type: ignore[assignment]
    MISSING_NOTEBOOK_DEPS.append(exc.name or "unknown")

TYPE_OPTIONS = ["riboseq", "rnaseq", "tiseq", "unknown"]
DOWNLOAD_METHODS = ["ena-ftp-wget", "ena-ftp-curl", "ena-ftp-ascp", "ncbi-prefetch"]


def find_repo_root(start: Path | None = None) -> Path:
    candidate = (start or Path.cwd()).resolve()
    for path in [candidate, *candidate.parents]:
        if (path / "main.nf").exists() and (path / "scripts").exists():
            return path
    raise RuntimeError("Could not locate repo root from current working directory.")


REPO_ROOT = find_repo_root()
DEMO_DIR = REPO_ROOT / "test_data" / "tutorial_demo_public_data"
SCRIPT_PATH = REPO_ROOT / "scripts" / "fetch_public_metadata.py"
OUTPUT_DIR = REPO_ROOT / "tutorial_outputs"
OUTPUT_DIR.mkdir(exist_ok=True)

STATE: dict[str, object] = {}


def read_table(path: Path, sep: str = "\t") -> pd.DataFrame:
    return pd.read_csv(path, sep=sep, dtype=str).fillna("")


def save_text_file(path: Path, content: str) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")
    return path


def clear_state(*keys: str) -> None:
    for key in keys:
        STATE.pop(key, None)


def normalize_accessions(raw_text: str) -> list[str]:
    accessions: list[str] = []
    seen: set[str] = set()
    for raw_line in raw_text.splitlines():
        line = raw_line.split("#", 1)[0]
        for token in re.split(r"[\s,;]+", line):
            cleaned = token.strip()
            if cleaned and cleaned not in seen:
                accessions.append(cleaned)
                seen.add(cleaned)
    return accessions


def render_message(text: str, kind: str = "info") -> None:
    color = {
        "info": "#1f5aa6",
        "warn": "#8a5a00",
        "error": "#a61f2d",
        "ok": "#1d6d2e",
    }.get(kind, "#1f5aa6")
    display(
        HTML(
            f"<div style='padding:8px 10px;border-left:4px solid {color};"
            f"background:#f7f9fc;margin:8px 0'>{text}</div>"
        )
    )


def select_columns(df: pd.DataFrame, columns: list[str]) -> pd.DataFrame:
    return df[[column for column in columns if column in df.columns]].copy()


def value_counts_table(series: pd.Series, label: str) -> pd.DataFrame:
    normalized = series.fillna("").astype(str).str.strip().replace("", "(blank)")
    counts = normalized.value_counts(dropna=False).rename_axis(label).reset_index(name="count")
    return counts


def render_cards(cards: list[tuple[str, str, str]]) -> widgets.HTML:
    blocks = []
    for title, value, note in cards:
        blocks.append(
            "<div style='flex:1 1 180px;min-width:180px;padding:12px 14px;"
            "border:1px solid #dbe4f0;border-radius:10px;background:#ffffff'>"
            f"<div style='font-size:12px;color:#5b6778;text-transform:uppercase;letter-spacing:.04em'>{html.escape(title)}</div>"
            f"<div style='font-size:24px;font-weight:700;color:#1c2f52;margin-top:4px'>{html.escape(value)}</div>"
            f"<div style='font-size:12px;color:#667085;margin-top:6px'>{html.escape(note)}</div>"
            "</div>"
        )
    return widgets.HTML(
        value=(
            "<div style='display:flex;flex-wrap:wrap;gap:12px;margin:8px 0 14px 0'>"
            + "".join(blocks)
            + "</div>"
        )
    )


def dataframe_html(df: pd.DataFrame, max_rows: int | None = None) -> str:
    preview = df.head(max_rows) if max_rows is not None else df
    return preview.to_html(index=False, escape=True)


def dataframe_accordion(title: str, df: pd.DataFrame) -> widgets.Accordion:
    output = widgets.Output()
    with output:
        display(df)
    accordion = widgets.Accordion(children=[output])
    accordion.set_title(0, title)
    return accordion


def current_state_df() -> pd.DataFrame | None:
    if "reviewed_df" in STATE:
        return STATE["reviewed_df"]  # type: ignore[return-value]
    if "metadata_df" in STATE:
        return STATE["metadata_df"]  # type: ignore[return-value]
    return None


def build_state_summary_html() -> str:
    metadata_df = current_state_df()
    metadata_note = "No metadata loaded yet."
    metadata_status = "Pending"
    if metadata_df is not None:
        review_count = int(metadata_df["needs_manual_review"].astype(str).str.lower().eq("true").sum())
        metadata_status = f"{len(metadata_df)} runs"
        metadata_note = (
            f"{review_count} rows still need manual review."
            if review_count
            else "Type and group look ready for downstream steps."
        )

    commands = STATE.get("download_commands")
    download_status = "Pending"
    download_note = "Generate commands after Step 0."
    if isinstance(commands, list) and commands:
        download_status = f"{len(commands)} commands"
        download_note = "Review or export the shell script before running."

    samplesheet_df = STATE.get("samplesheet_df")
    samplesheet_status = "Pending"
    samplesheet_note = "Build after all unknown sample types are resolved."
    if isinstance(samplesheet_df, pd.DataFrame):
        samplesheet_status = f"{len(samplesheet_df)} rows"
        samplesheet_note = f"Saved to {STATE.get('samplesheet_path', 'configured path')}."

    nextflow_command = STATE.get("nextflow_command")
    command_status = "Pending"
    command_note = "Command builder comes after samplesheet generation."
    if isinstance(nextflow_command, str) and nextflow_command:
        command_status = "Ready"
        command_note = "Suggested run command has been generated."

    cards = [
        ("Step 0 Metadata", metadata_status, metadata_note),
        ("Step 1 Downloads", download_status, download_note),
        ("Step 2 Samplesheet", samplesheet_status, samplesheet_note),
        ("Step 3 Command", command_status, command_note),
    ]
    card_html = render_cards(cards).value
    return (
        "<div style='padding:14px 16px;border:1px solid #ccd7e6;border-radius:12px;"
        "background:linear-gradient(135deg,#f7fbff 0%,#eef5ff 100%);margin:10px 0 18px 0'>"
        "<div style='font-size:18px;font-weight:700;color:#15345b'>Notebook state / 当前进度</div>"
        "<div style='font-size:13px;color:#526074;margin-top:4px'>"
        f"Guided outputs default to <code>{html.escape(str(OUTPUT_DIR))}</code>. "
        "Re-running an upstream step will clear downstream cached outputs so the notebook stays consistent."
        "</div>"
        f"{card_html}</div>"
    )


def build_metadata_summary_widget(df: pd.DataFrame) -> widgets.Widget:
    review_mask = df["needs_manual_review"].astype(str).str.lower().eq("true")
    unknown_mask = df["inferred_type"].astype(str).eq("unknown")
    blank_group_mask = df["suggested_group"].astype(str).str.strip().eq("")

    cards = render_cards(
        [
            ("Runs", str(len(df)), "Run-level records ready for review"),
            ("Need review", str(int(review_mask.sum())), "Rows flagged by heuristic checks"),
            ("Unknown type", str(int(unknown_mask.sum())), "Must be fixed before samplesheet build"),
            ("Groups", str(df["suggested_group"].astype(str).replace("", pd.NA).dropna().nunique()), "Distinct suggested groups"),
        ]
    )

    type_counts = value_counts_table(df["inferred_type"], "inferred_type")
    layout_counts = value_counts_table(df["library_layout"], "library_layout")
    group_counts = value_counts_table(df["suggested_group"], "suggested_group").head(8)
    preview = select_columns(
        df,
        [
            "run_accession",
            "sample_title",
            "library_strategy",
            "library_layout",
            "inferred_type",
            "suggested_group",
            "needs_manual_review",
        ],
    )

    warnings: list[str] = []
    if int(unknown_mask.sum()):
        warnings.append("At least one run is still `unknown`; Step 2 will intentionally block until you review it.")
    if int(blank_group_mask.sum()):
        warnings.append("Some runs have empty `group`; this is allowed, but replicate merge will not help those rows.")
    if not warnings:
        warnings.append("No blocking issues detected from the current heuristic summary.")

    summary_html = (
        "<div style='display:flex;flex-wrap:wrap;gap:18px'>"
        "<div style='flex:1 1 280px'>"
        "<h4 style='margin:4px 0'>Type summary</h4>"
        f"{dataframe_html(type_counts, max_rows=10)}"
        "</div>"
        "<div style='flex:1 1 220px'>"
        "<h4 style='margin:4px 0'>Layout summary</h4>"
        f"{dataframe_html(layout_counts, max_rows=10)}"
        "</div>"
        "<div style='flex:1 1 280px'>"
        "<h4 style='margin:4px 0'>Top groups</h4>"
        f"{dataframe_html(group_counts, max_rows=8)}"
        "</div>"
        "</div>"
        "<div style='margin-top:10px'>"
        "<h4 style='margin:4px 0'>What to check next</h4>"
        "<ul style='margin:6px 0 0 18px'>"
        + "".join(f"<li>{html.escape(item)}</li>" for item in warnings)
        + "</ul></div>"
        "<div style='margin-top:12px'>"
        "<h4 style='margin:4px 0'>Key columns for manual review</h4>"
        f"{dataframe_html(preview, max_rows=12)}"
        "</div>"
    )
    return widgets.VBox([cards, widgets.HTML(value=summary_html), dataframe_accordion("Full metadata table", df)])


def build_download_summary_widget(
    manifest: pd.DataFrame,
    commands: list[str],
    method: str,
    outdir: str,
) -> widgets.Widget:
    ftp_rows = int(manifest["ftp_url"].astype(str).ne("").sum())
    sra_rows = int(manifest["ftp_url"].astype(str).eq("").sum())
    cards = render_cards(
        [
            ("Files", str(len(manifest)), "Download manifest rows"),
            ("Commands", str(len(commands)), "Includes directory setup"),
            ("FTP-backed", str(ftp_rows), "Can pull directly from ENA FTP"),
            ("SRA fallback", str(sra_rows), "Will use prefetch/fasterq-dump"),
        ]
    )
    preview = select_columns(
        manifest,
        ["run_accession", "file_role", "file_name", "recommended_method", "fallback_sra_accession"],
    )
    notes = [
        f"Method selected: {method}",
        f"Output directory: {outdir}",
        "Recommended workflow: generate, inspect, export, then execute only if paths and quotas look right.",
    ]
    html_block = (
        "<div>"
        "<h4 style='margin:4px 0'>Manifest preview</h4>"
        f"{dataframe_html(preview, max_rows=12)}"
        "<h4 style='margin:12px 0 4px 0'>Review checklist</h4>"
        "<ul style='margin:6px 0 0 18px'>"
        + "".join(f"<li>{html.escape(note)}</li>" for note in notes)
        + "</ul></div>"
    )
    return widgets.VBox([cards, widgets.HTML(value=html_block), dataframe_accordion("Full download manifest", manifest)])


def build_samplesheet_summary_widget(df: pd.DataFrame, strandedness: str) -> widgets.Widget:
    group_counts = value_counts_table(df["group"], "group").head(10)
    type_counts = value_counts_table(df["type"], "type")
    paired_rows = int(df["fastq_2"].astype(str).ne("").sum())
    cards = render_cards(
        [
            ("Rows", str(len(df)), "Samplesheet rows written to disk"),
            ("Paired", str(paired_rows), "Rows with both FASTQ mates"),
            ("Single-end", str(len(df) - paired_rows), "Rows with only FASTQ_1"),
            ("Strandedness", strandedness, "Current notebook setting"),
        ]
    )
    html_block = (
        "<div style='display:flex;flex-wrap:wrap;gap:18px'>"
        "<div style='flex:1 1 220px'>"
        "<h4 style='margin:4px 0'>Type summary</h4>"
        f"{dataframe_html(type_counts, max_rows=10)}"
        "</div>"
        "<div style='flex:1 1 280px'>"
        "<h4 style='margin:4px 0'>Group summary</h4>"
        f"{dataframe_html(group_counts, max_rows=10)}"
        "</div>"
        "</div>"
        "<div style='margin-top:12px'>"
        "<h4 style='margin:4px 0'>Key samplesheet columns</h4>"
        f"{dataframe_html(select_columns(df, ['sample', 'fastq_1', 'fastq_2', 'strandedness', 'type', 'group']), max_rows=12)}"
        "</div>"
    )
    return widgets.VBox([cards, widgets.HTML(value=html_block), dataframe_accordion("Full samplesheet", df)])


def build_demo_orf_summary_widget(df: pd.DataFrame) -> widgets.Widget:
    multi_tool = int(df["tools"].astype(str).str.contains(",").sum())
    multi_sample = int(df["samples"].astype(str).str.contains(",").sum())
    top_tool = value_counts_table(df["tools"].astype(str).str.split(",").explode().str.strip(), "tool").head(1)
    top_tool_name = top_tool.iloc[0]["tool"] if not top_tool.empty else "n/a"
    cards = render_cards(
        [
            ("Demo ORFs", str(len(df)), "Rows in the bundled demo table"),
            ("Multi-tool", str(multi_tool), "Supported by more than one caller"),
            ("Multi-sample", str(multi_sample), "Detected in multiple samples"),
            ("Top tool", str(top_tool_name), "Most frequent caller label"),
        ]
    )
    ranked = df.copy()
    ranked["unique_psites_num"] = pd.to_numeric(ranked["unique_psites"], errors="coerce")
    ranked["length_aa_num"] = pd.to_numeric(ranked["length_aa"], errors="coerce")
    preview = ranked.sort_values(["unique_psites_num", "length_aa_num"], ascending=[False, False])
    preview = select_columns(
        preview,
        ["orf_id", "tools", "samples", "length_aa", "start_codon", "unique_psites", "pN"],
    ).head(10)
    html_block = (
        "<div>"
        "<h4 style='margin:4px 0'>High-support ORFs in the demo set</h4>"
        f"{dataframe_html(preview, max_rows=10)}"
        "</div>"
    )
    return widgets.VBox([cards, widgets.HTML(value=html_block), dataframe_accordion("Full demo ORF table", df)])


def load_demo_metadata() -> pd.DataFrame:
    return read_table(DEMO_DIR / "metadata_curated.tsv")


def load_demo_orfs() -> pd.DataFrame:
    return read_table(DEMO_DIR / "unified_orfs_demo.tsv")


def run_metadata_query(
    accessions: list[str],
    source_strategy: str,
    output_prefix: Path,
    progress_callback: Callable[[str], None] | None = None,
) -> pd.DataFrame:
    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        "python",
        str(SCRIPT_PATH),
        "--output-prefix",
        str(output_prefix),
        "--source-strategy",
        source_strategy,
        "--emit-download-manifest",
        "--emit-samplesheet-template",
    ]
    for accession in accessions:
        cmd.extend(["--accession", accession])
    process = subprocess.Popen(
        cmd,
        cwd=REPO_ROOT,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )

    start_time = time.time()
    while process.poll() is None:
        if progress_callback:
            elapsed = int(time.time() - start_time)
            progress_callback(f"Querying metadata from remote services... {elapsed}s elapsed")
        time.sleep(0.5)

    stdout, stderr = process.communicate()
    STATE["last_query_stdout"] = stdout
    STATE["last_query_stderr"] = stderr
    STATE["last_query_cmd"] = " ".join(shlex.quote(x) for x in cmd)
    if process.returncode != 0:
        raise RuntimeError(stderr or stdout or "Metadata query failed.")

    if progress_callback:
        progress_callback("Finalizing query outputs...")

    return read_table(Path(f"{output_prefix}.metadata_curated.tsv"))


def build_review_widgets(df: pd.DataFrame, refresh_callback: Callable[[], None] | None = None) -> widgets.Widget:
    items = []
    controls = []
    for idx, row in df.iterrows():
        header = widgets.HTML(
            value=(
                f"<b>{row['run_accession']}</b> | <code>{row['sample_title']}</code><br>"
                f"<span style='color:#444'>strategy={row['library_strategy']} | "
                f"inferred={row['inferred_type']} | suggested_group={row['suggested_group']}</span><br>"
                f"<span style='color:#666'>type evidence: {row['type_evidence']} | "
                f"group evidence: {row['group_evidence']}</span>"
            )
        )
        type_widget = widgets.Dropdown(
            options=TYPE_OPTIONS,
            value=row["inferred_type"] if row["inferred_type"] in TYPE_OPTIONS else "unknown",
            description="type",
            layout=widgets.Layout(width="260px"),
        )
        group_widget = widgets.Text(
            value=row["suggested_group"],
            description="group",
            layout=widgets.Layout(width="420px"),
        )
        items.append(widgets.VBox([header, widgets.HBox([type_widget, group_widget])]))
        controls.append((idx, type_widget, group_widget))

    apply_btn = widgets.Button(description="Apply review edits", button_style="success")
    export_path = widgets.Text(
        value=str(OUTPUT_DIR / "reviewed_metadata.csv"),
        description="CSV path",
        layout=widgets.Layout(width="720px"),
    )
    export_btn = widgets.Button(description="Export review CSV", button_style="primary")
    out = widgets.Output()

    def collect_reviewed() -> pd.DataFrame:
        reviewed = df.copy()
        for idx, type_widget, group_widget in controls:
            reviewed.at[idx, "inferred_type"] = type_widget.value
            reviewed.at[idx, "suggested_group"] = group_widget.value.strip()
            reviewed.at[idx, "needs_manual_review"] = "true" if type_widget.value == "unknown" else "false"
        return reviewed

    def on_apply(_: widgets.Button) -> None:
        reviewed = collect_reviewed()
        clear_state(
            "download_manifest_df",
            "download_commands",
            "download_command_path",
            "samplesheet_df",
            "samplesheet_path",
            "nextflow_command",
        )
        STATE["reviewed_df"] = reviewed
        if refresh_callback:
            refresh_callback()
        with out:
            clear_output()
            render_message("Reviewed metadata saved in STATE['reviewed_df']", "ok")
            display(
                reviewed[
                    ["run_accession", "sample_title", "inferred_type", "suggested_group", "needs_manual_review"]
                ]
            )

    def on_export(_: widgets.Button) -> None:
        reviewed = collect_reviewed()
        path = save_text_file(Path(export_path.value), reviewed.to_csv(index=False))
        clear_state(
            "download_manifest_df",
            "download_commands",
            "download_command_path",
            "samplesheet_df",
            "samplesheet_path",
            "nextflow_command",
        )
        STATE["reviewed_df"] = reviewed
        STATE["review_export_path"] = path
        if refresh_callback:
            refresh_callback()
        with out:
            clear_output()
            render_message(f"Review CSV written to <code>{path}</code>", "ok")
            display(
                reviewed[
                    ["run_accession", "sample_title", "inferred_type", "suggested_group", "needs_manual_review"]
                ]
            )

    apply_btn.on_click(on_apply)
    export_btn.on_click(on_export)
    return widgets.VBox(items + [widgets.HBox([apply_btn, export_btn]), export_path, out])


def current_metadata_df() -> pd.DataFrame:
    if "reviewed_df" in STATE:
        return STATE["reviewed_df"].copy()  # type: ignore[return-value]
    if "metadata_df" in STATE:
        return STATE["metadata_df"].copy()  # type: ignore[return-value]
    raise RuntimeError("No metadata is loaded yet. Run Step 0 first.")


def build_download_manifest_from_metadata(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for _, row in df.iterrows():
        if row.get("fastq_ftp_1"):
            rows.append(
                {
                    "run_accession": row["run_accession"],
                    "file_role": "R1",
                    "file_name": Path(row["fastq_ftp_1"]).name,
                    "ftp_url": row["fastq_ftp_1"],
                    "md5": row["fastq_md5_1"],
                    "recommended_method": "ena-ftp",
                    "fallback_sra_accession": row["sra_download_accession"],
                }
            )
        if row.get("fastq_ftp_2"):
            rows.append(
                {
                    "run_accession": row["run_accession"],
                    "file_role": "R2",
                    "file_name": Path(row["fastq_ftp_2"]).name,
                    "ftp_url": row["fastq_ftp_2"],
                    "md5": row["fastq_md5_2"],
                    "recommended_method": "ena-ftp",
                    "fallback_sra_accession": row["sra_download_accession"],
                }
            )
        if not row.get("fastq_ftp_1") and not row.get("fastq_ftp_2"):
            rows.append(
                {
                    "run_accession": row["run_accession"],
                    "file_role": "SRA",
                    "file_name": f"{row['sra_download_accession']}.sra",
                    "ftp_url": "",
                    "md5": "",
                    "recommended_method": "ncbi-sra",
                    "fallback_sra_accession": row["sra_download_accession"],
                }
            )
    return pd.DataFrame(rows)


def download_commands(df: pd.DataFrame, method: str, outdir: str = "data") -> list[str]:
    commands = [f"mkdir -p {shlex.quote(outdir)}"]
    manifest = build_download_manifest_from_metadata(df)
    for _, row in manifest.iterrows():
        dest = f"{outdir}/{row['file_name']}"
        if method == "ena-ftp-wget" and row["ftp_url"]:
            commands.append(f"wget -c ftp://{row['ftp_url']} -O {shlex.quote(dest)}")
        elif method == "ena-ftp-curl" and row["ftp_url"]:
            commands.append(f"curl -L ftp://{row['ftp_url']} -o {shlex.quote(dest)}")
        elif method == "ena-ftp-ascp" and row["ftp_url"]:
            commands.append(
                "ascp -QT -l 300m -P33001 "
                f"era-fasp@fasp.sra.ebi.ac.uk:{row['ftp_url']} {shlex.quote(outdir)}/"
            )
        else:
            acc = row["fallback_sra_accession"]
            commands.append(
                f"prefetch {acc} && fasterq-dump --split-files --threads 8 {acc} && pigz -p 8 {acc}*.fastq"
            )
    return commands


def build_samplesheet(df: pd.DataFrame, strandedness: str = "auto") -> pd.DataFrame:
    if (df["inferred_type"] == "unknown").any():
        unknown_runs = ", ".join(df.loc[df["inferred_type"] == "unknown", "run_accession"].tolist())
        raise ValueError(f"Unknown sample type remains for: {unknown_runs}. Please review Step 0 first.")
    rows = []
    name_counts = Counter()
    for _, row in df.iterrows():
        sample = row["sample_title"] or row["run_accession"]
        safe = "".join(ch if ch.isalnum() or ch in "._-" else "_" for ch in sample)
        name_counts[safe] += 1
        if name_counts[safe] > 1:
            safe = f"{safe}_{name_counts[safe]}"
        rows.append(
            {
                "sample": safe,
                "fastq_1": row["fastq_ftp_1"] or row["sra_download_accession"],
                "fastq_2": row["fastq_ftp_2"],
                "strandedness": strandedness,
                "type": row["inferred_type"],
                "group": row["suggested_group"],
                "run_accession": row["run_accession"],
                "input_accession": row["input_accession"],
            }
        )
    return pd.DataFrame(rows)


def save_samplesheet(df: pd.DataFrame, path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)
    STATE["samplesheet_path"] = path
    STATE["samplesheet_df"] = df.copy()
    return path


def build_nextflow_command(
    samplesheet_path: Path,
    outdir: str,
    profile: str,
    aligner: str,
    merge_replicates: bool,
    skip_unify: bool,
    skip_classify: bool,
    input_mode: str = "fastq",
    resume: bool = False,
) -> str:
    cmd = ["nextflow", "run", ".", "-profile", profile, "--outdir", outdir]
    if input_mode == "fastq":
        cmd.extend(["--input", str(samplesheet_path)])
    else:
        cmd.extend(["--input", "samplesheet_bam.csv"])
    cmd.extend(["--aligner", aligner])
    if merge_replicates:
        cmd.append("--merge_replicates")
    if skip_unify:
        cmd.extend(["--skip_unify_orf_predictions", "true"])
    if skip_classify:
        cmd.extend(["--skip_orf_classification", "true"])
    if resume:
        cmd.append("-resume")
    return " ".join(shlex.quote(part) for part in cmd)


def plot_demo_orf_summary(df: pd.DataFrame) -> None:
    exploded_tools = []
    for cell in df["tools"]:
        exploded_tools.extend([tool.strip() for tool in str(cell).split(",") if tool.strip()])
    tool_counts = pd.Series(exploded_tools).value_counts().sort_values(ascending=False)
    plt.figure(figsize=(6, 3.5))
    tool_counts.plot(kind="bar", color=["#2f6db5", "#f28e2b", "#59a14f"])
    plt.title("Demo ORFs by tool / 示例 ORF 来源工具")
    plt.ylabel("Count")
    plt.tight_layout()
    plt.show()


def launch_tutorial() -> None:
    if MISSING_NOTEBOOK_DEPS:
        raise RuntimeError(
            "Notebook dependencies are missing. Please install at least: "
            + ", ".join(sorted(set(MISSING_NOTEBOOK_DEPS + ["pandas", "ipywidgets", "matplotlib"])))
        )

    render_message(f"Repo root detected at: <code>{REPO_ROOT}</code>", "ok")
    render_message("Run Step 0 to Step 4 in order. The notebook is optimized for teaching, review, and reproducible command generation.", "info")

    intro = widgets.HTML(
        value=(
            "<div style='padding:18px 20px;border:1px solid #dbe4f0;border-radius:14px;"
            "background:linear-gradient(135deg,#fffdf7 0%,#f4f8ff 100%)'>"
            "<h2 style='margin:0 0 8px 0'>Ribo-seq Guided Tutorial</h2>"
            "<p style='margin:0 0 12px 0;color:#475467'>"
            "This notebook is a guided launcher for new users: public-data discovery, metadata review, "
            "download preparation, samplesheet generation, Nextflow command assembly, and demo ORF orientation."
            "</p>"
            "<div style='display:flex;flex-wrap:wrap;gap:18px'>"
            "<div style='flex:1 1 240px'><b>What you will produce</b><br>"
            "<span style='color:#526074'>Reviewed metadata CSV, download shell script, samplesheet CSV, run command.</span></div>"
            "<div style='flex:1 1 240px'><b>Recommended rhythm</b><br>"
            "<span style='color:#526074'>Review first, export second, execute last. The notebook intentionally blocks risky steps when metadata is incomplete.</span></div>"
            "<div style='flex:1 1 240px'><b>Good follow-up reading</b><br>"
            "<span style='color:#526074'><code>docs/lab_tutorial.md</code>, <code>docs/usage.md</code>, <code>docs/output.md</code></span></div>"
            "</div></div>"
        )
    )
    state_panel = widgets.HTML()

    def refresh_state_panel() -> None:
        state_panel.value = build_state_summary_html()

    refresh_state_panel()

    accession_box = widgets.Textarea(
        value=(DEMO_DIR / "accessions.txt").read_text(encoding="utf-8").strip(),
        description="Accessions",
        layout=widgets.Layout(width="720px", height="120px"),
    )
    source_dropdown = widgets.Dropdown(
        options=["ena-first", "ncbi-first", "ena-only", "ncbi-only"],
        value="ena-first",
        description="Source",
    )
    use_demo_checkbox = widgets.Checkbox(value=True, description="Use bundled demo results if available")
    query_button = widgets.Button(description="Step 0: Query metadata", button_style="primary")
    step0_out = widgets.Output()

    def on_query(_: widgets.Button) -> None:
        with step0_out:
            clear_output()
            status_bar = widgets.IntProgress(
                value=0,
                min=0,
                max=100,
                description="Status",
                bar_style="info",
                layout=widgets.Layout(width="720px"),
            )
            status_text = widgets.HTML(
                "<span style='color:#1f5aa6'>Preparing metadata query...</span>"
            )
            display(widgets.VBox([status_bar, status_text]))

            def update_status(message: str) -> None:
                status_text.value = f"<span style='color:#1f5aa6'>{message}</span>"

            raw_accessions = [line for line in accession_box.value.splitlines() if line.strip()]
            accessions = normalize_accessions(accession_box.value)
            accession_box.value = "\n".join(accessions)
            if not accessions:
                status_bar.bar_style = "danger"
                status_bar.value = 100
                render_message("Please provide at least one accession.", "error")
                return
            if len(accessions) != len(raw_accessions):
                render_message("Accessions were normalized by removing duplicates, comments, or extra separators.", "info")
            demo_accessions = [
                line.strip()
                for line in (DEMO_DIR / "accessions.txt").read_text(encoding="utf-8").splitlines()
                if line.strip()
            ]
            if use_demo_checkbox.value and accessions == demo_accessions:
                update_status("Loading bundled demo metadata...")
                df = load_demo_metadata()
                status_bar.bar_style = "success"
                status_bar.value = 100
                status_text.value = "<span style='color:#1d6d2e'>Done. Bundled metadata loaded.</span>"
                render_message("Loaded bundled demo metadata from test_data/tutorial_demo_public_data/.", "ok")
            else:
                prefix = OUTPUT_DIR / "notebook_public_metadata"
                try:
                    update_status("Launching metadata query command...")
                    df = run_metadata_query(accessions, source_dropdown.value, prefix, progress_callback=update_status)
                    status_bar.bar_style = "success"
                    status_bar.value = 100
                    status_text.value = "<span style='color:#1d6d2e'>Done. Metadata query completed.</span>"
                    render_message("Metadata query completed successfully.", "ok")
                    print(STATE["last_query_cmd"])
                except Exception as exc:
                    status_bar.bar_style = "danger"
                    status_bar.value = 100
                    status_text.value = "<span style='color:#a61f2d'>Failed. Check error details below.</span>"
                    render_message(str(exc), "error")
                    if STATE.get("last_query_stderr"):
                        print(STATE["last_query_stderr"])
                    return
            clear_state(
                "reviewed_df",
                "review_export_path",
                "download_manifest_df",
                "download_commands",
                "download_command_path",
                "samplesheet_df",
                "samplesheet_path",
                "nextflow_command",
            )
            STATE["metadata_df"] = df.copy()
            STATE["requested_accessions"] = accessions
            refresh_state_panel()
            display(build_metadata_summary_widget(df))
            render_message("Review inferred type and suggested group below.", "info")
            display(build_review_widgets(df, refresh_callback=refresh_state_panel))

    query_button.on_click(on_query)

    download_method = widgets.Dropdown(options=DOWNLOAD_METHODS, value="ena-ftp-wget", description="Method")
    download_dir = widgets.Text(value="data", description="Out dir")
    show_downloads_btn = widgets.Button(description="Step 1: Generate download commands", button_style="primary")
    download_cmd_path = widgets.Text(
        value=str(OUTPUT_DIR / "download_commands.sh"),
        description="Cmd path",
        layout=widgets.Layout(width="720px"),
    )
    export_downloads_btn = widgets.Button(description="Export commands to file", button_style="success")
    allow_exec = widgets.Checkbox(value=False, description="I understand this may download real data")
    run_download_btn = widgets.Button(description="Run download commands", button_style="warning")
    download_out = widgets.Output()

    def on_generate_downloads(_: widgets.Button) -> None:
        with download_out:
            clear_output()
            try:
                df = current_metadata_df()
            except Exception as exc:
                render_message(str(exc), "error")
                return
            manifest = build_download_manifest_from_metadata(df)
            STATE["download_manifest_df"] = manifest.copy()
            outdir = download_dir.value.strip() or "data"
            commands = download_commands(df, download_method.value, outdir)
            STATE["download_commands"] = commands
            refresh_state_panel()
            render_message("Review the generated commands before executing them.", "warn")
            display(build_download_summary_widget(manifest, commands, download_method.value, outdir))
            print("\n".join(commands))

    def on_export_downloads(_: widgets.Button) -> None:
        with download_out:
            clear_output()
            try:
                df = current_metadata_df()
            except Exception as exc:
                render_message(str(exc), "error")
                return
            outdir = download_dir.value.strip() or "data"
            commands = download_commands(df, download_method.value, outdir)
            command_text = "#!/usr/bin/env bash\nset -euo pipefail\n\n" + "\n".join(commands) + "\n"
            path = save_text_file(Path(download_cmd_path.value), command_text)
            STATE["download_commands"] = commands
            STATE["download_command_path"] = path
            refresh_state_panel()
            render_message(f"Download commands written to <code>{path}</code>", "ok")
            print("\n".join(commands))

    def on_run_downloads(_: widgets.Button) -> None:
        with download_out:
            if not allow_exec.value:
                render_message("Enable the confirmation checkbox before running downloads.", "error")
                return
            commands = STATE.get("download_commands")
            if not commands:
                render_message("Generate download commands first.", "error")
                return
            render_message("Executing download commands one by one. Interrupt the kernel if you need to stop.", "warn")
            for command in commands:  # type: ignore[assignment]
                print("$", command)
                subprocess.run(command, cwd=REPO_ROOT, shell=True, check=False)

    show_downloads_btn.on_click(on_generate_downloads)
    export_downloads_btn.on_click(on_export_downloads)
    run_download_btn.on_click(on_run_downloads)

    strandedness_widget = widgets.Dropdown(
        options=["auto", "forward", "reverse", "unstranded"],
        value="auto",
        description="Strand",
    )
    samplesheet_path_widget = widgets.Text(
        value=str(OUTPUT_DIR / "public_samplesheet.csv"),
        description="Save path",
        layout=widgets.Layout(width="720px"),
    )
    build_samplesheet_btn = widgets.Button(description="Step 2: Build samplesheet", button_style="primary")
    samplesheet_out = widgets.Output()

    def on_build_samplesheet(_: widgets.Button) -> None:
        with samplesheet_out:
            clear_output()
            try:
                df = current_metadata_df()
                samplesheet_df = build_samplesheet(df, strandedness_widget.value)
                path = save_samplesheet(samplesheet_df, Path(samplesheet_path_widget.value))
            except Exception as exc:
                render_message(str(exc), "error")
                return
            clear_state("nextflow_command")
            refresh_state_panel()
            render_message(f"Samplesheet written to <code>{path}</code>", "ok")
            display(build_samplesheet_summary_widget(samplesheet_df, strandedness_widget.value))

    build_samplesheet_btn.on_click(on_build_samplesheet)

    profile_widget = widgets.Text(value="docker", description="Profile")
    aligner_widget = widgets.Dropdown(options=["star", "hisat2"], value="star", description="Aligner")
    outdir_widget = widgets.Text(
        value="results/tutorial_run",
        description="Outdir",
        layout=widgets.Layout(width="360px"),
    )
    merge_widget = widgets.Checkbox(value=False, description="merge replicates")
    skip_unify_widget = widgets.Checkbox(value=False, description="skip unify")
    skip_classify_widget = widgets.Checkbox(value=False, description="skip classify")
    resume_widget = widgets.Checkbox(value=False, description="append -resume")
    input_mode_widget = widgets.Dropdown(options=["fastq", "bam"], value="fastq", description="Input mode")
    build_cmd_btn = widgets.Button(description="Step 3: Build nextflow command", button_style="primary")
    command_out = widgets.Output()

    def on_build_command(_: widgets.Button) -> None:
        with command_out:
            clear_output()
            samplesheet_path = STATE.get("samplesheet_path")
            if input_mode_widget.value == "fastq" and not samplesheet_path:
                render_message("Please generate a samplesheet in Step 2 first.", "error")
                return
            command = build_nextflow_command(
                samplesheet_path=samplesheet_path or Path("samplesheet.csv"),  # type: ignore[arg-type]
                outdir=outdir_widget.value,
                profile=profile_widget.value,
                aligner=aligner_widget.value,
                merge_replicates=merge_widget.value,
                skip_unify=skip_unify_widget.value,
                skip_classify=skip_classify_widget.value,
                input_mode=input_mode_widget.value,
                resume=resume_widget.value,
            )
            STATE["nextflow_command"] = command
            refresh_state_panel()
            render_message("Suggested command generated below.", "ok")
            print(command)
            current_df = current_state_df()
            if current_df is not None:
                riboseq_rows = int(current_df["inferred_type"].astype(str).eq("riboseq").sum())
                if riboseq_rows == 0:
                    render_message("Current reviewed metadata does not include any `riboseq` rows. Double-check that this run setup matches your teaching goal.", "warn")
                if merge_widget.value:
                    mergeable_groups = int(current_df["suggested_group"].astype(str).value_counts().gt(1).sum())
                    if mergeable_groups == 0:
                        render_message("`merge_replicates` is enabled, but no repeated groups were found in the current metadata.", "warn")
            print("Replace docker with singularity or another local profile if needed.")
            print("For BAM input mode, remember to provide a compatible BAM samplesheet and explicit strandedness.")

    build_cmd_btn.on_click(on_build_command)

    demo_orfs_btn = widgets.Button(description="Step 4: Show demo ORF outputs", button_style="primary")
    demo_out = widgets.Output()

    def on_show_demo(_: widgets.Button) -> None:
        with demo_out:
            clear_output()
            demo_orfs = load_demo_orfs()
            display(build_demo_orf_summary_widget(demo_orfs))
            plot_demo_orf_summary(demo_orfs)
            render_message(
                "Tip: tools tells you which ORF caller supported the candidate. "
                "unique_psites and pN are useful support metrics to compare ORFs.",
                "info",
            )
            stats_lines = (DEMO_DIR / "unified_orfs.stats.txt").read_text(encoding="utf-8").splitlines()
            for line in stats_lines[:25]:
                print(line)

    demo_orfs_btn.on_click(on_show_demo)

    layout = widgets.VBox(
        [
            intro,
            state_panel,
            widgets.HTML("<h3>Step 0: Find public data</h3>"),
            widgets.HTML("<p style='color:#526074;margin-top:-6px'>Paste SRR/ERR/DRR/SRP/PRJNA-style accessions. Comments and duplicates are cleaned automatically.</p>"),
            accession_box,
            widgets.HBox([source_dropdown, use_demo_checkbox, query_button]),
            step0_out,
            widgets.HTML("<h3>Step 1: Generate download commands</h3>"),
            widgets.HTML("<p style='color:#526074;margin-top:-6px'>This step is review-first: preview the manifest, export a script, then run only if you truly want to fetch data from public services.</p>"),
            widgets.HBox([download_method, download_dir, show_downloads_btn]),
            download_cmd_path,
            export_downloads_btn,
            widgets.HBox([allow_exec, run_download_btn]),
            download_out,
            widgets.HTML("<h3>Step 2: Build samplesheet</h3>"),
            widgets.HTML("<p style='color:#526074;margin-top:-6px'>The notebook blocks this step when any row is still `unknown`, so that downstream analyses are not silently misconfigured.</p>"),
            widgets.HBox([strandedness_widget, build_samplesheet_btn]),
            samplesheet_path_widget,
            samplesheet_out,
            widgets.HTML("<h3>Step 3: Build nextflow command</h3>"),
            widgets.HBox([input_mode_widget, profile_widget, aligner_widget]),
            outdir_widget,
            widgets.HBox([merge_widget, skip_unify_widget, skip_classify_widget, resume_widget]),
            build_cmd_btn,
            command_out,
            widgets.HTML("<h3>Step 4: Explore demo outputs</h3>"),
            widgets.HTML("<p style='color:#526074;margin-top:-6px'>Use the bundled ORF demo to learn what a unified ORF table looks like before you wait for a full pipeline run.</p>"),
            demo_orfs_btn,
            demo_out,
        ]
    )
    display(layout)
