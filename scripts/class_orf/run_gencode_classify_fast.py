#!/usr/bin/env python3
"""
Process-parallel entrypoint for GENCODE classification.

The fast path keeps the reference mapper as the source of truth for annotation
semantics, but parallelizes safely by splitting the ORF input into independent
chromosome bins and running one mapper subprocess per bin.
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Sequence, Tuple


CODON_TABLE = {
    "TTT": "F",
    "TTC": "F",
    "TTA": "L",
    "TTG": "L",
    "CTT": "L",
    "CTC": "L",
    "CTA": "L",
    "CTG": "L",
    "ATT": "I",
    "ATC": "I",
    "ATA": "I",
    "ATG": "M",
    "GTT": "V",
    "GTC": "V",
    "GTA": "V",
    "GTG": "V",
    "TCT": "S",
    "TCC": "S",
    "TCA": "S",
    "TCG": "S",
    "CCT": "P",
    "CCC": "P",
    "CCA": "P",
    "CCG": "P",
    "ACT": "T",
    "ACC": "T",
    "ACA": "T",
    "ACG": "T",
    "GCT": "A",
    "GCC": "A",
    "GCA": "A",
    "GCG": "A",
    "TAT": "Y",
    "TAC": "Y",
    "TAA": "*",
    "TAG": "*",
    "CAT": "H",
    "CAC": "H",
    "CAA": "Q",
    "CAG": "Q",
    "AAT": "N",
    "AAC": "N",
    "AAA": "K",
    "AAG": "K",
    "GAT": "D",
    "GAC": "D",
    "GAA": "E",
    "GAG": "E",
    "TGT": "C",
    "TGC": "C",
    "TGA": "*",
    "TGG": "W",
    "CGT": "R",
    "CGC": "R",
    "CGA": "R",
    "CGG": "R",
    "AGT": "S",
    "AGC": "S",
    "AGA": "R",
    "AGG": "R",
    "GGT": "G",
    "GGC": "G",
    "GGA": "G",
    "GGG": "G",
}

MAIN_OUTPUT_SUFFIXES = [
    ".orfs.out",
    ".orfs.gtf",
    ".orfs.bed",
    ".orfs.fa",
    ".orfs.frames.bed",
    ".orfs.allframes.bed",
    ".logs",
]
TRANSCRIPT_ID_RE = re.compile(r'transcript_id "([^"]+)"')
PROTEIN_ID_RE = re.compile(r'protein_id "([^"]+)"')


def eprint(message: str) -> None:
    print(message, file=sys.stderr)


def find_script(script_name: str) -> Path:
    here = Path(__file__).resolve().parent
    candidates = [
        here.parent / "gencode-riboseqORFs" / script_name,
        here / script_name,
        here.parent / script_name,
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    raise FileNotFoundError(f"Could not locate {script_name}")


def translate_nt(nt_seq: str) -> str:
    nt_seq = nt_seq.upper().replace("U", "T")
    aa: List[str] = []
    for i in range(0, len(nt_seq) - 2, 3):
        aa.append(CODON_TABLE.get(nt_seq[i : i + 3], "X"))
    return "".join(aa)


def load_orf_to_study_and_sequences(
    metadata_file: Path,
) -> Tuple[Dict[str, str], Dict[str, str]]:
    orf_to_study: Dict[str, str] = {}
    orf_to_nt: Dict[str, str] = {}

    with metadata_file.open(newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        header_map = {name: idx for idx, name in enumerate(header)}
        orf_idx = header_map["orf_id"]
        samples_idx = header_map.get("samples", -1)
        seq_idx = header_map.get("sequence", -1)
        for row in reader:
            if len(row) <= orf_idx:
                continue
            orf_id = row[orf_idx]
            study_id = "unified"
            if 0 <= samples_idx < len(row) and row[samples_idx]:
                study_id = row[samples_idx].split(",")[0]
            orf_to_study[orf_id] = study_id
            if 0 <= seq_idx < len(row):
                orf_to_nt[orf_id] = row[seq_idx]
    return orf_to_study, orf_to_nt


def write_bed6(bed12_file: Path, bed6_file: Path, orf_to_study: Dict[str, str]) -> int:
    count = 0
    with bed12_file.open() as src, bed6_file.open("w") as dst:
        for line in src:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue
            chrom, start, end, name, _, strand = parts[:6]
            dst.write(
                f"{chrom}\t{start}\t{end}\t{name}\t{orf_to_study.get(name, 'unified')}\t{strand}\n"
            )
            count += 1
    return count


def write_protein_fasta(
    fasta_file: Path,
    orf_to_study: Dict[str, str],
    orf_to_nt: Dict[str, str],
) -> int:
    count = 0
    with fasta_file.open("w") as handle:
        for orf_id, nt_seq in orf_to_nt.items():
            aa_seq = translate_nt(nt_seq)
            handle.write(f">{orf_id}--{orf_to_study.get(orf_id, 'unified')}\n{aa_seq}\n")
            count += 1
    return count


def iter_fasta_records(fasta_path: Path) -> Iterator[Tuple[str, str]]:
    header = None
    seq_lines: List[str] = []
    with fasta_path.open() as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_lines)
                header = line[1:]
                seq_lines = []
            else:
                seq_lines.append(line)
    if header is not None:
        yield header, "".join(seq_lines)


def load_fasta_by_orf_id(fasta_path: Path) -> Dict[str, Tuple[str, str]]:
    records: Dict[str, Tuple[str, str]] = {}
    for header, seq in iter_fasta_records(fasta_path):
        orf_id = header.split("--", 1)[0]
        records[orf_id] = (header, seq)
    return records


def load_bed6_groups(bed6_file: Path) -> Tuple[List[str], Dict[str, List[str]]]:
    chrom_order: List[str] = []
    groups: Dict[str, List[str]] = {}
    with bed6_file.open() as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if not line:
                continue
            chrom = line.split("\t", 1)[0]
            if chrom not in groups:
                chrom_order.append(chrom)
                groups[chrom] = []
            groups[chrom].append(line)
    return chrom_order, groups


def assign_bins(chrom_order: Sequence[str], groups: Dict[str, List[str]], max_bins: int) -> List[List[str]]:
    if not chrom_order:
        return []
    n_bins = max(1, min(max_bins, len(chrom_order)))
    bins: List[List[str]] = [[] for _ in range(n_bins)]
    loads = [0] * n_bins
    for chrom in sorted(chrom_order, key=lambda key: len(groups[key]), reverse=True):
        idx = min(range(n_bins), key=lambda i: loads[i])
        bins[idx].append(chrom)
        loads[idx] += len(groups[chrom])
    chrom_pos = {chrom: i for i, chrom in enumerate(chrom_order)}
    for chroms in bins:
        chroms.sort(key=lambda chrom: chrom_pos[chrom])
    return [chroms for chroms in bins if chroms]


def build_chunk_inputs(
    chunk_dir: Path,
    chunk_index: int,
    chroms: Sequence[str],
    bed_groups: Dict[str, List[str]],
    fasta_records: Dict[str, Tuple[str, str]],
) -> Tuple[Path, Path, int, int]:
    chunk_dir.mkdir(parents=True, exist_ok=True)
    bed6_path = chunk_dir / f"chunk_{chunk_index:02d}.orfs.bed6"
    fasta_path = chunk_dir / f"chunk_{chunk_index:02d}.orfs.fa"
    orf_ids: List[str] = []
    skipped = 0
    with bed6_path.open("w") as bed_handle:
        for chrom in chroms:
            for line in bed_groups[chrom]:
                fields = line.split("\t")
                if len(fields) >= 4:
                    orf_id = fields[3]
                    if orf_id not in fasta_records:
                        skipped += 1
                        continue
                    orf_ids.append(orf_id)
                    bed_handle.write(line + "\n")

    with fasta_path.open("w") as fasta_handle:
        for orf_id in orf_ids:
            record = fasta_records.get(orf_id)
            if record is None:
                continue
            header, seq = record
            fasta_handle.write(f">{header}\n{seq}\n")
    return bed6_path, fasta_path, len(orf_ids), skipped


def parse_reference_paths(ensembl_dir: Path) -> Dict[str, Path]:
    return {
        "TRANSCRIPTOME_FASTA": ensembl_dir / "TRANSCRIPTOME_FASTA",
        "SORTED_TRANSCRIPTOME_GTF": ensembl_dir / "SORTED_TRANSCRIPTOME_GTF",
        "PROTEOME_FASTA": ensembl_dir / "PROTEOME_FASTA",
        "TRANSCRIPT_SUPPORT": ensembl_dir / "TRANSCRIPT_SUPPORT",
        "PSITES_BED": ensembl_dir / "PSITES_BED",
    }


def filter_fasta_by_ids(input_path: Path, output_path: Path, keep_ids: Sequence[str]) -> None:
    keep = set(keep_ids)
    with output_path.open("w") as out_handle:
        for header, seq in iter_fasta_records(input_path):
            record_id = header.split(None, 1)[0]
            if record_id in keep:
                out_handle.write(f">{record_id}\n{seq}\n")


def build_filtered_reference(
    source_dir: Path,
    target_dir: Path,
    chroms: Sequence[str],
) -> Tuple[int, int]:
    target_dir.mkdir(parents=True, exist_ok=True)
    paths = parse_reference_paths(source_dir)
    chrom_set = set(chroms)
    transcript_ids: set[str] = set()
    protein_ids: set[str] = set()

    gtf_in = paths["SORTED_TRANSCRIPTOME_GTF"]
    gtf_out = target_dir / "SORTED_TRANSCRIPTOME_GTF"
    with gtf_in.open() as src, gtf_out.open("w") as dst:
        for line in src:
            if line.startswith("#"):
                dst.write(line)
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[0] not in chrom_set:
                continue
            dst.write(line)
            attrs = fields[8]
            tx_match = TRANSCRIPT_ID_RE.search(attrs)
            if tx_match:
                transcript_ids.add(tx_match.group(1))
            protein_match = PROTEIN_ID_RE.search(attrs)
            if protein_match:
                protein_ids.add(protein_match.group(1))

    filter_fasta_by_ids(
        paths["TRANSCRIPTOME_FASTA"],
        target_dir / "TRANSCRIPTOME_FASTA",
        sorted(transcript_ids),
    )
    filter_fasta_by_ids(
        paths["PROTEOME_FASTA"],
        target_dir / "PROTEOME_FASTA",
        sorted(protein_ids),
    )

    support_out = target_dir / "TRANSCRIPT_SUPPORT"
    with paths["TRANSCRIPT_SUPPORT"].open() as src, support_out.open("w") as dst:
        for line in src:
            record_id = line.split("\t", 1)[0]
            if record_id in transcript_ids:
                dst.write(line)

    psites_out = target_dir / "PSITES_BED"
    with paths["PSITES_BED"].open() as src, psites_out.open("w") as dst:
        for line in src:
            if not line.strip():
                continue
            chrom = line.split("\t", 1)[0]
            if chrom in chrom_set:
                dst.write(line)

    return len(transcript_ids), len(protein_ids)


def prepare_chunk(
    chunk_index: int,
    chroms: Sequence[str],
    chunk_root: Path,
    bed_groups: Dict[str, List[str]],
    fasta_records: Dict[str, Tuple[str, str]],
    ensembl_dir: Path,
) -> Dict[str, object]:
    chunk_dir = chunk_root / f"chunk_{chunk_index:02d}"
    build_inputs_started = time.perf_counter()
    chunk_bed6, chunk_fasta, n_orfs, n_skipped = build_chunk_inputs(
        chunk_dir=chunk_dir,
        chunk_index=chunk_index,
        chroms=chroms,
        bed_groups=bed_groups,
        fasta_records=fasta_records,
    )
    build_inputs_seconds = time.perf_counter() - build_inputs_started

    build_ref_started = time.perf_counter()
    ref_dir = chunk_dir / "ensembl_dir"
    n_tx, n_prot = build_filtered_reference(
        source_dir=ensembl_dir,
        target_dir=ref_dir,
        chroms=chroms,
    )
    build_ref_seconds = time.perf_counter() - build_ref_started

    return {
        "chunk_index": chunk_index,
        "chroms": list(chroms),
        "bed6": chunk_bed6,
        "fasta": chunk_fasta,
        "ensembl_dir": ref_dir,
        "output_prefix": chunk_dir / "gencode_results",
        "profile_path": chunk_dir / "gencode_mapper.profile.tsv",
        "n_orfs": n_orfs,
        "n_skipped": n_skipped,
        "n_tx": n_tx,
        "n_prot": n_prot,
        "build_inputs_seconds": build_inputs_seconds,
        "build_ref_seconds": build_ref_seconds,
    }


def run_reference_mapper(
    mapper_script: Path,
    ensembl_dir: Path,
    fasta_file: Path,
    bed6_file: Path,
    output_prefix: Path,
    profile_path: Path,
) -> None:
    cmd = [
        "python3",
        str(mapper_script),
        "-d",
        str(ensembl_dir),
        "-f",
        str(fasta_file),
        "-b",
        str(bed6_file),
        "-o",
        str(output_prefix),
    ]
    env = os.environ.copy()
    env.pop("GENCODE_FAST_CPUS", None)
    env["GENCODE_PROFILE_PATH"] = str(profile_path)
    subprocess.run(cmd, check=True, env=env)


def merge_files(files: Sequence[Path], destination: Path, keep_header: bool) -> None:
    with destination.open("w") as out_handle:
        wrote_header = False
        for src in files:
            if not src.exists():
                continue
            with src.open() as in_handle:
                for line_number, line in enumerate(in_handle):
                    if keep_header and line_number == 0:
                        if wrote_header:
                            continue
                        wrote_header = True
                    out_handle.write(line)


def merge_profiles(chunk_rows: Sequence[Tuple[int, str, float, str]], profile_file: Path) -> None:
    with profile_file.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["chunk", "stage", "seconds", "details"])
        for row in chunk_rows:
            writer.writerow(row)


def append_chunk_profile(chunk_id: int, profile_path: Path, rows: List[Tuple[int, str, float, str]]) -> None:
    if not profile_path.exists():
        return
    with profile_path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            try:
                seconds = float(row.get("seconds", "0") or 0)
            except ValueError:
                seconds = 0.0
            rows.append((chunk_id, row.get("stage", ""), seconds, row.get("details", "")))


def write_profile(profile_file: Path, rows: Iterable[Tuple[str, float, str]]) -> None:
    with profile_file.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["stage", "seconds", "details"])
        for row in rows:
            writer.writerow(row)


def main() -> None:
    parser = argparse.ArgumentParser(description="Fast GENCODE classifier entrypoint")
    parser.add_argument("--input_prefix", required=True, help="Unified ORF prefix")
    parser.add_argument("--ensembl_dir", required=True, help="Prepared Ensembl directory")
    parser.add_argument("--output_dir", required=True, help="Output directory")
    parser.add_argument("--cpus", type=int, default=1, help="Number of mapper subprocesses")
    args = parser.parse_args()

    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    input_prefix = Path(args.input_prefix)
    metadata_file = Path(f"{input_prefix}.metadata.tsv")
    bed12_file = Path(f"{input_prefix}.bed")
    fasta_file = Path(f"{input_prefix}.orfs.fa")
    bed6_file = Path(f"{input_prefix}.orfs.bed6")
    output_prefix = output_dir / "gencode_results"
    profile_file = output_dir / "gencode_fast.profile.tsv"
    chunk_root = output_dir / "chunk_runs"

    for required in (metadata_file, bed12_file):
        if not required.exists():
            raise FileNotFoundError(f"Required input not found: {required}")

    timings: List[Tuple[str, float, str]] = []

    t0 = time.perf_counter()
    orf_to_study, orf_to_nt = load_orf_to_study_and_sequences(metadata_file)
    timings.append(("load_metadata", time.perf_counter() - t0, f"orfs={len(orf_to_study)}"))

    if not bed6_file.exists():
        t0 = time.perf_counter()
        n_bed = write_bed6(bed12_file, bed6_file, orf_to_study)
        timings.append(("write_bed6", time.perf_counter() - t0, f"records={n_bed}"))

    if not fasta_file.exists():
        t0 = time.perf_counter()
        n_fa = write_protein_fasta(fasta_file, orf_to_study, orf_to_nt)
        timings.append(("write_protein_fasta", time.perf_counter() - t0, f"records={n_fa}"))

    t0 = time.perf_counter()
    mapper_script = find_script("ORF_mapper_to_GENCODE_v1.1.py")
    fasta_records = load_fasta_by_orf_id(fasta_file)
    chrom_order, bed_groups = load_bed6_groups(bed6_file)
    bins = assign_bins(chrom_order, bed_groups, max(1, args.cpus))
    timings.append(
        (
            "prepare_chunks",
            time.perf_counter() - t0,
            f"chromosomes={len(chrom_order)} bins={len(bins)}",
        )
    )

    source_ensembl_dir = Path(args.ensembl_dir).resolve()
    chunk_specs: List[Dict[str, object]] = []
    t0 = time.perf_counter()
    with ThreadPoolExecutor(max_workers=max(1, min(args.cpus, len(bins)))) as executor:
        futures = {
            executor.submit(
                prepare_chunk,
                chunk_index,
                chroms,
                chunk_root,
                bed_groups,
                fasta_records,
                source_ensembl_dir,
            ): (chunk_index, chroms)
            for chunk_index, chroms in enumerate(bins, start=1)
        }
        for future in as_completed(futures):
            spec = future.result()
            chunk_specs.append(spec)
    timings.append(("prepare_chunk_artifacts", time.perf_counter() - t0, f"chunks={len(chunk_specs)}"))
    for spec in sorted(chunk_specs, key=lambda item: item["chunk_index"]):
        timings.append(
            (
                "build_chunk_inputs",
                float(spec["build_inputs_seconds"]),
                f"chunk={spec['chunk_index']} chroms={','.join(spec['chroms'])} orfs={spec['n_orfs']} skipped={spec['n_skipped']}",
            )
        )
        timings.append(
            (
                "build_chunk_reference",
                float(spec["build_ref_seconds"]),
                f"chunk={spec['chunk_index']} transcripts={spec['n_tx']} proteins={spec['n_prot']}",
            )
        )

    t0 = time.perf_counter()
    with ThreadPoolExecutor(max_workers=max(1, min(args.cpus, len(chunk_specs)))) as executor:
        futures = {
            executor.submit(
                run_reference_mapper,
                mapper_script,
                spec["ensembl_dir"],
                spec["fasta"],
                spec["bed6"],
                spec["output_prefix"],
                spec["profile_path"],
            ): spec
            for spec in chunk_specs
        }
        for future in as_completed(futures):
            spec = futures[future]
            future.result()
            eprint(
                "Completed GENCODE chunk "
                f"{spec['chunk_index']} ({','.join(spec['chroms'])})"
            )
    timings.append(("parallel_reference_mapper", time.perf_counter() - t0, f"chunks={len(chunk_specs)}"))

    t0 = time.perf_counter()
    ordered_specs = sorted(chunk_specs, key=lambda spec: spec["chunk_index"])
    merge_files(
        [Path(f"{spec['output_prefix']}.orfs.out") for spec in ordered_specs],
        Path(f"{output_prefix}.orfs.out"),
        keep_header=True,
    )
    for suffix in MAIN_OUTPUT_SUFFIXES[1:]:
        merge_files(
            [Path(f"{spec['output_prefix']}{suffix}") for spec in ordered_specs],
            Path(f"{output_prefix}{suffix}"),
            keep_header=(suffix == ".logs"),
        )
    chunk_profile_rows: List[Tuple[int, str, float, str]] = []
    for spec in ordered_specs:
        append_chunk_profile(spec["chunk_index"], spec["profile_path"], chunk_profile_rows)
    merge_profiles(chunk_profile_rows, output_dir / "gencode_mapper.profile.tsv")
    timings.append(("merge_outputs", time.perf_counter() - t0, f"chunks={len(ordered_specs)}"))

    write_profile(profile_file, timings)


if __name__ == "__main__":
    main()
