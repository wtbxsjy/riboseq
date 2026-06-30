#!/usr/bin/env python3
"""
Prepare deterministic unified-ORF subsets for GENCODE classifier regression tests.
"""

from __future__ import annotations

import argparse
import csv
import random
import re
import sys
from pathlib import Path
from typing import List, Set


csv.field_size_limit(sys.maxsize)


def read_ids_from_bed(path: Path) -> List[str]:
    ids: List[str] = []
    with path.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 4:
                ids.append(parts[3])
    return ids


def choose_ids(ids: List[str], count: int, seed: int) -> List[str]:
    if count >= len(ids):
        return list(ids)
    rng = random.Random(seed)
    return sorted(rng.sample(ids, count))


def filter_text_lines(src: Path, dst: Path, keep_ids: Set[str], id_column: int) -> int:
    kept = 0
    with src.open() as in_handle, dst.open("w") as out_handle:
        for line in in_handle:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) > id_column and parts[id_column] in keep_ids:
                out_handle.write(line)
                kept += 1
    return kept


def filter_gtf(src: Path, dst: Path, keep_ids: Set[str]) -> int:
    kept = 0
    orf_pattern = re.compile(r'orf_id "([^"]+)"')
    tx_pattern = re.compile(r'transcript_id "([^"]+)"')
    with src.open() as in_handle, dst.open("w") as out_handle:
        for line in in_handle:
            if line.startswith("#"):
                out_handle.write(line)
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            match = orf_pattern.search(parts[8])
            if match is None:
                match = tx_pattern.search(parts[8])
            if match and match.group(1) in keep_ids:
                out_handle.write(line)
                kept += 1
    return kept


def filter_metadata(src: Path, dst: Path, keep_ids: Set[str]) -> int:
    kept = 0
    with src.open(newline="") as in_handle, dst.open("w", newline="") as out_handle:
        reader = csv.reader(in_handle, delimiter="\t")
        writer = csv.writer(out_handle, delimiter="\t", lineterminator="\n")
        header = next(reader)
        writer.writerow(header)
        orf_idx = header.index("orf_id")
        for row in reader:
            if len(row) > orf_idx and row[orf_idx] in keep_ids:
                writer.writerow(row)
                kept += 1
    return kept


def filter_fasta(src: Path, dst: Path, keep_ids: Set[str]) -> int:
    kept = 0
    write_record = False
    with src.open() as in_handle, dst.open("w") as out_handle:
        for line in in_handle:
            if line.startswith(">"):
                record_id = line[1:].strip().split()[0].split("--")[0]
                write_record = record_id in keep_ids
                if write_record:
                    kept += 1
                    out_handle.write(line)
            elif write_record:
                out_handle.write(line)
    return kept


def main() -> None:
    parser = argparse.ArgumentParser(description="Prepare unified ORF subset for regression tests")
    parser.add_argument("--input_prefix", required=True, help="Input unified ORF prefix")
    parser.add_argument("--output_prefix", required=True, help="Output unified ORF prefix")
    parser.add_argument("--count", type=int, required=True, help="Number of ORFs to keep")
    parser.add_argument("--seed", type=int, default=17, help="Sampling seed")
    parser.add_argument("--source_bed", help="Optional bed/bed6 file to sample IDs from")
    args = parser.parse_args()

    input_prefix = Path(args.input_prefix)
    output_prefix = Path(args.output_prefix)
    output_prefix.parent.mkdir(parents=True, exist_ok=True)

    src_bed = Path(args.source_bed) if args.source_bed else Path(f"{input_prefix}.bed")
    metadata = Path(f"{input_prefix}.metadata.tsv")
    bed = Path(f"{input_prefix}.bed")
    gtf = Path(f"{input_prefix}.gtf")
    fasta = Path(f"{input_prefix}.orfs.fa")

    all_ids = read_ids_from_bed(src_bed)
    keep_ids = set(choose_ids(all_ids, args.count, args.seed))

    summary_lines = []
    summary_lines.append(f"selected_orfs\t{len(keep_ids)}")

    if metadata.exists():
        kept = filter_metadata(metadata, Path(f"{output_prefix}.metadata.tsv"), keep_ids)
        summary_lines.append(f"metadata_rows\t{kept}")
    if bed.exists():
        kept = filter_text_lines(bed, Path(f"{output_prefix}.bed"), keep_ids, 3)
        summary_lines.append(f"bed_rows\t{kept}")
    if gtf.exists():
        kept = filter_gtf(gtf, Path(f"{output_prefix}.gtf"), keep_ids)
        summary_lines.append(f"gtf_rows\t{kept}")
    if fasta.exists():
        kept = filter_fasta(fasta, Path(f"{output_prefix}.orfs.fa"), keep_ids)
        summary_lines.append(f"fasta_records\t{kept}")

    with Path(f"{output_prefix}.subset_summary.tsv").open("w") as handle:
        handle.write("\n".join(summary_lines) + "\n")


if __name__ == "__main__":
    main()
