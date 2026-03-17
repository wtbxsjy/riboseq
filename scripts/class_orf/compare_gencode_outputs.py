#!/usr/bin/env python3
"""
Compare key GENCODE classifier outputs and emit a compact regression summary.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Dict, List, Tuple


KEY_COLUMNS = ["orf_id", "trans", "gene", "gene_name", "orf_biotype", "gene_biotype"]


def load_orfs_out(path: Path) -> Dict[str, Dict[str, str]]:
    records: Dict[str, Dict[str, str]] = {}
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            orf_id = row.get("orf_id")
            if not orf_id or orf_id.startswith("#"):
                continue
            records[orf_id] = {key: row.get(key, "") for key in KEY_COLUMNS}
    return records


def compare_rows(left: Dict[str, Dict[str, str]], right: Dict[str, Dict[str, str]]) -> Tuple[int, int, List[str]]:
    left_ids = set(left)
    right_ids = set(right)
    shared = sorted(left_ids & right_ids)
    mismatches: List[str] = []
    for orf_id in shared:
        for key in KEY_COLUMNS[1:]:
            if left[orf_id].get(key, "") != right[orf_id].get(key, ""):
                mismatches.append(
                    f"{orf_id}\t{key}\t{left[orf_id].get(key, '')}\t{right[orf_id].get(key, '')}"
                )
    return len(left_ids - right_ids), len(right_ids - left_ids), mismatches


def main() -> None:
    parser = argparse.ArgumentParser(description="Compare key gencode_results.orfs.out fields")
    parser.add_argument("--left", required=True, help="Reference .orfs.out")
    parser.add_argument("--right", required=True, help="Candidate .orfs.out")
    parser.add_argument("--summary", required=True, help="Output summary TSV")
    args = parser.parse_args()

    left = load_orfs_out(Path(args.left))
    right = load_orfs_out(Path(args.right))
    missing_left, missing_right, mismatches = compare_rows(left, right)

    summary_path = Path(args.summary)
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    with summary_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["metric", "value"])
        writer.writerow(["left_records", len(left)])
        writer.writerow(["right_records", len(right)])
        writer.writerow(["missing_in_right", missing_left])
        writer.writerow(["missing_in_left", missing_right])
        writer.writerow(["field_mismatches", len(mismatches)])

    if mismatches:
        mismatch_path = summary_path.with_suffix(".mismatches.tsv")
        with mismatch_path.open("w") as handle:
            handle.write("orf_id\tfield\tleft\tright\n")
            handle.write("\n".join(mismatches) + "\n")


if __name__ == "__main__":
    main()
