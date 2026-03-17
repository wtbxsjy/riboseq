#!/usr/bin/env python3
"""
Build a synthetic unified-ORF input set directly from an Ensembl directory.

This is useful for benchmarking GENCODE classification on a non-empty dataset
without depending on full upstream ORF-caller outputs.
"""

from __future__ import annotations

import argparse
import csv
import random
import re
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple


ATTR_RE = re.compile(r'(\S+) "([^"]+)"')


def parse_attrs(attr_text: str) -> Dict[str, str]:
    return {key: value for key, value in ATTR_RE.findall(attr_text)}


def load_proteome(proteome_fasta: Path) -> Dict[str, str]:
    seqs: Dict[str, str] = {}
    current_id: Optional[str] = None
    chunks: List[str] = []
    with proteome_fasta.open() as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if current_id is not None:
                    seqs[current_id] = "".join(chunks)
                current_id = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
    if current_id is not None:
        seqs[current_id] = "".join(chunks)
    return seqs


def collect_cds_entries(gtf_path: Path) -> List[Dict[str, object]]:
    transcripts: Dict[str, Dict[str, object]] = {}
    with gtf_path.open() as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "CDS":
                continue
            attrs = parse_attrs(parts[8])
            tx_id = attrs.get("transcript_id")
            protein_id = attrs.get("protein_id")
            if not tx_id or not protein_id:
                continue
            entry = transcripts.setdefault(
                tx_id,
                {
                    "chrom": parts[0],
                    "strand": parts[6],
                    "gene_id": attrs.get("gene_id", "NA"),
                    "gene_name": attrs.get("gene_name", attrs.get("gene_id", "NA")),
                    "protein_id": protein_id,
                    "starts": [],
                    "ends": [],
                },
            )
            entry["starts"].append(int(parts[3]))
            entry["ends"].append(int(parts[4]))

    rows: List[Dict[str, object]] = []
    for tx_id, entry in transcripts.items():
        starts = sorted(entry["starts"])
        ends = sorted(entry["ends"])
        rows.append(
            {
                "transcript_id": tx_id,
                "protein_id": entry["protein_id"],
                "chrom": entry["chrom"],
                "strand": entry["strand"],
                "gene_id": entry["gene_id"],
                "gene_name": entry["gene_name"],
                "start": min(starts),
                "end": max(ends),
                "exon_blocks": ";".join(f"{s}-{e}" for s, e in zip(starts, ends)),
            }
        )
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description="Build synthetic unified input from Ensembl references")
    parser.add_argument("--ensembl_dir", required=True, help="Prepared Ensembl directory")
    parser.add_argument("--output_prefix", required=True, help="Output prefix for .bed/.metadata.tsv/.orfs.fa")
    parser.add_argument("--count", type=int, default=100, help="Number of ORFs to sample")
    parser.add_argument("--seed", type=int, default=17, help="Sampling seed")
    parser.add_argument("--sample_name", default="synthetic", help="Sample/study ID to embed")
    args = parser.parse_args()

    ensembl_dir = Path(args.ensembl_dir)
    output_prefix = Path(args.output_prefix)
    output_prefix.parent.mkdir(parents=True, exist_ok=True)

    gtf_path = ensembl_dir / "SORTED_TRANSCRIPTOME_GTF"
    proteome_path = ensembl_dir / "PROTEOME_FASTA"
    if not gtf_path.exists() or not proteome_path.exists():
        raise FileNotFoundError("Ensembl dir must contain SORTED_TRANSCRIPTOME_GTF and PROTEOME_FASTA")

    proteome = load_proteome(proteome_path)
    cds_rows = [row for row in collect_cds_entries(gtf_path) if row["protein_id"] in proteome]
    if not cds_rows:
        raise RuntimeError("No CDS entries with matching protein sequences were found")

    rng = random.Random(args.seed)
    selected = cds_rows if args.count >= len(cds_rows) else rng.sample(cds_rows, args.count)
    selected = sorted(selected, key=lambda row: (str(row["chrom"]), int(row["start"]), str(row["transcript_id"])))

    metadata_path = Path(f"{output_prefix}.metadata.tsv")
    bed_path = Path(f"{output_prefix}.bed")
    fasta_path = Path(f"{output_prefix}.orfs.fa")

    with metadata_path.open("w", newline="") as meta_handle, bed_path.open("w") as bed_handle, fasta_path.open("w") as fasta_handle:
        writer = csv.writer(meta_handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "orf_id",
                "chrom",
                "strand",
                "start",
                "end",
                "length_aa",
                "exon_blocks",
                "gene_id",
                "transcript_id",
                "tools",
                "samples",
                "tool_scores",
                "tool_pvalues",
                "total_reads",
                "unique_reads",
                "total_psites",
                "unique_psites",
                "pN",
                "unique_pN",
                "sequence",
            ]
        )

        for idx, row in enumerate(selected, start=1):
            protein_seq = proteome[str(row["protein_id"])]
            orf_id = f"SYNTH_ORF_{idx}_{row['transcript_id']}"
            nt_len = max(0, len(protein_seq.replace("*", "")) * 3)
            writer.writerow(
                [
                    orf_id,
                    row["chrom"],
                    row["strand"],
                    row["start"],
                    row["end"],
                    len(protein_seq.replace("*", "")),
                    row["exon_blocks"],
                    row["gene_id"],
                    row["transcript_id"],
                    "Synthetic",
                    args.sample_name,
                    "Synthetic:1.0",
                    "Synthetic:0",
                    0,
                    0,
                    0,
                    0,
                    "0.000000",
                    "0.000000",
                    "ATG",
                ]
            )
            # BED file only needs first 6 columns for wrapper -> mapper conversion.
            bed_handle.write(
                f"{row['chrom']}\t{int(row['start']) - 1}\t{row['end']}\t{orf_id}\t0\t{row['strand']}\n"
            )
            fasta_handle.write(f">{orf_id}--{args.sample_name}\n{protein_seq}\n")


if __name__ == "__main__":
    main()
