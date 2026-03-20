from __future__ import annotations

import shlex
import subprocess
from collections import Counter
from pathlib import Path

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


def load_demo_metadata() -> pd.DataFrame:
    return read_table(DEMO_DIR / "metadata_curated.tsv")


def load_demo_orfs() -> pd.DataFrame:
    return read_table(DEMO_DIR / "unified_orfs_demo.tsv")


def run_metadata_query(accessions: list[str], source_strategy: str, output_prefix: Path) -> pd.DataFrame:
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
    result = subprocess.run(cmd, cwd=REPO_ROOT, capture_output=True, text=True)
    STATE["last_query_stdout"] = result.stdout
    STATE["last_query_stderr"] = result.stderr
    STATE["last_query_cmd"] = " ".join(shlex.quote(x) for x in cmd)
    if result.returncode != 0:
        raise RuntimeError(result.stderr or result.stdout or "Metadata query failed.")
    return read_table(Path(f"{output_prefix}.metadata_curated.tsv"))


def build_review_widgets(df: pd.DataFrame) -> widgets.Widget:
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
    out = widgets.Output()

    def on_apply(_: widgets.Button) -> None:
        reviewed = df.copy()
        for idx, type_widget, group_widget in controls:
            reviewed.at[idx, "inferred_type"] = type_widget.value
            reviewed.at[idx, "suggested_group"] = group_widget.value.strip()
            reviewed.at[idx, "needs_manual_review"] = (
                "true" if type_widget.value == "unknown" else reviewed.at[idx, "needs_manual_review"]
            )
        STATE["reviewed_df"] = reviewed
        with out:
            clear_output()
            render_message("Reviewed metadata saved in STATE['reviewed_df']", "ok")
            display(
                reviewed[
                    ["run_accession", "sample_title", "inferred_type", "suggested_group", "needs_manual_review"]
                ]
            )

    apply_btn.on_click(on_apply)
    return widgets.VBox(items + [apply_btn, out])


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
    commands = []
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
    render_message("Run the widgets from top to bottom. This notebook is a guided launcher, not a full GUI.", "info")

    intro = widgets.HTML(
        value=(
            "<h2>Ribo-seq Guided Tutorial</h2>"
            "<p>Step 0 查询公开数据, Step 1 生成下载命令, Step 2 生成 samplesheet, "
            "Step 3 生成 nextflow 命令, Step 4 浏览 demo ORF 结果。</p>"
        )
    )

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
            accessions = [line.strip() for line in accession_box.value.splitlines() if line.strip()]
            if not accessions:
                render_message("Please provide at least one accession.", "error")
                return
            demo_accessions = [
                line.strip()
                for line in (DEMO_DIR / "accessions.txt").read_text(encoding="utf-8").splitlines()
                if line.strip()
            ]
            if use_demo_checkbox.value and accessions == demo_accessions:
                df = load_demo_metadata()
                render_message("Loaded bundled demo metadata from test_data/tutorial_demo_public_data/.", "ok")
            else:
                prefix = OUTPUT_DIR / "notebook_public_metadata"
                try:
                    df = run_metadata_query(accessions, source_dropdown.value, prefix)
                    render_message("Metadata query completed successfully.", "ok")
                    print(STATE["last_query_cmd"])
                except Exception as exc:
                    render_message(str(exc), "error")
                    if STATE.get("last_query_stderr"):
                        print(STATE["last_query_stderr"])
                    return
            STATE["metadata_df"] = df.copy()
            display(df)
            render_message("Review inferred type and suggested group below.", "info")
            display(build_review_widgets(df))

    query_button.on_click(on_query)

    download_method = widgets.Dropdown(options=DOWNLOAD_METHODS, value="ena-ftp-wget", description="Method")
    download_dir = widgets.Text(value="data", description="Out dir")
    show_downloads_btn = widgets.Button(description="Step 1: Generate download commands", button_style="primary")
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
            commands = download_commands(df, download_method.value, download_dir.value.strip() or "data")
            STATE["download_commands"] = commands
            render_message("Review the generated commands before executing them.", "warn")
            display(manifest)
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
            render_message(f"Samplesheet written to <code>{path}</code>", "ok")
            display(samplesheet_df)

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
            )
            STATE["nextflow_command"] = command
            render_message("Suggested command generated below.", "ok")
            print(command)
            print("Add -resume to continue a failed run.")
            print("Replace docker with singularity or another local profile if needed.")

    build_cmd_btn.on_click(on_build_command)

    demo_orfs_btn = widgets.Button(description="Step 4: Show demo ORF outputs", button_style="primary")
    demo_out = widgets.Output()

    def on_show_demo(_: widgets.Button) -> None:
        with demo_out:
            clear_output()
            demo_orfs = load_demo_orfs()
            display(demo_orfs.head(10))
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
            widgets.HTML("<h3>Step 0: Find public data</h3>"),
            accession_box,
            widgets.HBox([source_dropdown, use_demo_checkbox, query_button]),
            step0_out,
            widgets.HTML("<h3>Step 1: Generate download commands</h3>"),
            widgets.HBox([download_method, download_dir, show_downloads_btn]),
            widgets.HBox([allow_exec, run_download_btn]),
            download_out,
            widgets.HTML("<h3>Step 2: Build samplesheet</h3>"),
            widgets.HBox([strandedness_widget, build_samplesheet_btn]),
            samplesheet_path_widget,
            samplesheet_out,
            widgets.HTML("<h3>Step 3: Build nextflow command</h3>"),
            widgets.HBox([input_mode_widget, profile_widget, aligner_widget]),
            outdir_widget,
            widgets.HBox([merge_widget, skip_unify_widget, skip_classify_widget]),
            build_cmd_btn,
            command_out,
            widgets.HTML("<h3>Step 4: Explore demo outputs</h3>"),
            demo_orfs_btn,
            demo_out,
        ]
    )
    display(layout)
