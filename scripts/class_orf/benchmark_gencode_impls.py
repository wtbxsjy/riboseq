#!/usr/bin/env python3
"""
Run original / fast / indexed_fast GENCODE classification implementations on
the same input prefix and summarize runtime plus key output diffs.
"""

from __future__ import annotations

import argparse
import csv
import os
import shutil
import subprocess
import time
from pathlib import Path


STANDARD_REF_MAP = {
    "TRANSCRIPTOME_FASTA": [
        "TRANSCRIPTOME_FASTA",
        "Mus_musculus.GRCm39.transcriptome.fa",
    ],
    "SORTED_TRANSCRIPTOME_GTF": [
        "SORTED_TRANSCRIPTOME_GTF",
        "Mus_musculus.GRCm39.110.sorted.chr.gtf",
        "Mus_musculus.GRCm39.110.sorted.gtf",
    ],
    "PROTEOME_FASTA": [
        "PROTEOME_FASTA",
        "Mus_musculus.GRCm39.pep.all.fa",
    ],
    "TRANSCRIPT_SUPPORT": [
        "TRANSCRIPT_SUPPORT",
        "transcript_support_level.txt",
    ],
    "PSITES_BED": [
        "PSITES_BED",
        "psites.chr.bed",
        "psites.bed",
    ],
}


def prepare_standard_ensembl_dir(source_dir: Path, target_dir: Path) -> Path:
    target_dir.mkdir(parents=True, exist_ok=True)
    missing = []
    for standard_name, candidates in STANDARD_REF_MAP.items():
        source_path = None
        for candidate in candidates:
            candidate_path = source_dir / candidate
            if candidate_path.exists():
                source_path = candidate_path
                break
        if source_path is None or not source_path.exists():
            missing.append(standard_name)
            continue
        link_path = target_dir / standard_name
        if link_path.exists() or link_path.is_symlink():
            link_path.unlink()
        link_path.symlink_to(source_path.resolve())
    if missing:
        raise FileNotFoundError(
            f"Missing required Ensembl files for benchmark: {', '.join(missing)}"
        )
    return target_dir


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
    parser = argparse.ArgumentParser(description="Benchmark GENCODE classification implementations")
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
    indexed_fast_dir = workdir / "indexed_fast"
    diff_fast_summary = workdir / "compare.fast.summary.tsv"
    diff_indexed_fast_summary = workdir / "compare.indexed_fast.summary.tsv"
    runtime_summary = workdir / "runtime.summary.tsv"
    standard_ensembl_dir = prepare_standard_ensembl_dir(
        Path(args.ensembl_dir).resolve(),
        workdir / "ensembl_standard",
    )

    original_seconds = run_impl(
        repo_root=repo_root,
        impl="original",
        input_prefix=args.input_prefix,
        ensembl_dir=str(standard_ensembl_dir),
        out_dir=original_dir,
        cpus=1,
    )
    fast_seconds = run_impl(
        repo_root=repo_root,
        impl="fast",
        input_prefix=args.input_prefix,
        ensembl_dir=str(standard_ensembl_dir),
        out_dir=fast_dir,
        cpus=args.cpus,
    )
    indexed_fast_seconds = run_impl(
        repo_root=repo_root,
        impl="indexed_fast",
        input_prefix=args.input_prefix,
        ensembl_dir=str(standard_ensembl_dir),
        out_dir=indexed_fast_dir,
        cpus=args.cpus,
    )

    compare_outputs(
        repo_root=repo_root,
        left=original_dir / "gencode_results.orfs.out",
        right=fast_dir / "gencode_results.orfs.out",
        summary=diff_fast_summary,
    )
    compare_outputs(
        repo_root=repo_root,
        left=original_dir / "gencode_results.orfs.out",
        right=indexed_fast_dir / "gencode_results.orfs.out",
        summary=diff_indexed_fast_summary,
    )

    with runtime_summary.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["impl", "seconds"])
        writer.writerow(["original", f"{original_seconds:.6f}"])
        writer.writerow(["fast", f"{fast_seconds:.6f}"])
        writer.writerow(["indexed_fast", f"{indexed_fast_seconds:.6f}"])
        if fast_seconds > 0:
            writer.writerow(["fast_speedup_vs_original", f"{original_seconds / fast_seconds:.6f}"])
        if indexed_fast_seconds > 0:
            writer.writerow(["indexed_fast_speedup_vs_original", f"{original_seconds / indexed_fast_seconds:.6f}"])
        if indexed_fast_seconds > 0 and fast_seconds > 0:
            writer.writerow(["indexed_fast_speedup_vs_fast", f"{fast_seconds / indexed_fast_seconds:.6f}"])


if __name__ == "__main__":
    main()
