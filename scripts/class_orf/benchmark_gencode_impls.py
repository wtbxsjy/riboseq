#!/usr/bin/env python3
"""
Run original and fast GENCODE classification implementations on the same input
prefix and summarize runtime plus key output diffs.
"""

from __future__ import annotations

import argparse
import csv
import os
import shutil
import subprocess
import time
from pathlib import Path


def run_impl(
    repo_root: Path,
    impl: str,
    input_prefix: str,
    ensembl_dir: str,
    out_dir: Path,
    cpus: int,
) -> float:
    out_dir.mkdir(parents=True, exist_ok=True)
    cmd = [
        "python3",
        str(repo_root / "scripts" / "classify_orfs_wrapper.py"),
        "--mode",
        "gencode",
        "--input",
        input_prefix,
        "--output_dir",
        str(out_dir),
        "--ensembl_dir",
        ensembl_dir,
        "--gencode_impl",
        impl,
        "--cpus",
        str(cpus),
    ]
    t0 = time.perf_counter()
    env = os.environ.copy()
    env["GENCODE_PROFILE_PATH"] = str(out_dir / "gencode_mapper.profile.tsv")
    subprocess.run(cmd, check=True, cwd=repo_root, env=env)
    return time.perf_counter() - t0


def compare_outputs(repo_root: Path, left: Path, right: Path, summary: Path) -> None:
    cmd = [
        "python3",
        str(repo_root / "scripts" / "class_orf" / "compare_gencode_outputs.py"),
        "--left",
        str(left),
        "--right",
        str(right),
        "--summary",
        str(summary),
    ]
    subprocess.run(cmd, check=True, cwd=repo_root)


def main() -> None:
    parser = argparse.ArgumentParser(description="Benchmark original vs fast GENCODE classification")
    parser.add_argument("--repo_root", default=".", help="Repository root")
    parser.add_argument("--input_prefix", required=True, help="Unified ORF input prefix")
    parser.add_argument("--ensembl_dir", required=True, help="Prepared Ensembl directory")
    parser.add_argument("--workdir", required=True, help="Benchmark output directory")
    parser.add_argument("--cpus", type=int, default=4, help="CPU count for fast mode")
    parser.add_argument("--keep_existing", action="store_true", help="Do not clear previous benchmark directory")
    args = parser.parse_args()

    repo_root = Path(args.repo_root).resolve()
    workdir = Path(args.workdir).resolve()
    if workdir.exists() and not args.keep_existing:
        shutil.rmtree(workdir)
    workdir.mkdir(parents=True, exist_ok=True)

    original_dir = workdir / "original"
    fast_dir = workdir / "fast"
    diff_summary = workdir / "compare.summary.tsv"
    runtime_summary = workdir / "runtime.summary.tsv"

    original_seconds = run_impl(
        repo_root=repo_root,
        impl="original",
        input_prefix=args.input_prefix,
        ensembl_dir=args.ensembl_dir,
        out_dir=original_dir,
        cpus=1,
    )
    fast_seconds = run_impl(
        repo_root=repo_root,
        impl="fast",
        input_prefix=args.input_prefix,
        ensembl_dir=args.ensembl_dir,
        out_dir=fast_dir,
        cpus=args.cpus,
    )

    compare_outputs(
        repo_root=repo_root,
        left=original_dir / "gencode_results.orfs.out",
        right=fast_dir / "gencode_results.orfs.out",
        summary=diff_summary,
    )

    with runtime_summary.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["impl", "seconds"])
        writer.writerow(["original", f"{original_seconds:.6f}"])
        writer.writerow(["fast", f"{fast_seconds:.6f}"])
        if fast_seconds > 0:
            writer.writerow(["speedup_vs_original", f"{original_seconds / fast_seconds:.6f}"])


if __name__ == "__main__":
    main()
