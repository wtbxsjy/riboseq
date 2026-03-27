#!/usr/bin/env python3
"""
Indexed GENCODE classification entrypoint.

This implementation avoids rebuilding filtered Ensembl reference directories
for each parallel chunk. Instead, it builds a once-per-run chromosome-sharded
reference cache and lets worker processes consume the cached transcript state
directly.
"""

from __future__ import annotations

import argparse
import csv
import os
import pickle
import random
import re
import string
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Sequence, Tuple

from Bio.SeqIO.FastaIO import SimpleFastaParser

HERE = Path(__file__).resolve().parent
REPO_SCRIPTS = HERE.parent
GENCODE_DIR = REPO_SCRIPTS / "gencode-riboseqORFs"
sys.path.insert(0, str(GENCODE_DIR))
import functions  # noqa: E402


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

TRANSCRIPT_ID_RE = re.compile(r'transcript_id "([^"]+)"')
PROTEIN_ID_RE = re.compile(r'protein_id "([^"]+)"')
GENE_ID_RE = re.compile(r'gene_id "([^"]+)"')
GENE_NAME_RE = re.compile(r'gene_name "([^"]+)"')
GENE_BIOTYPE_RE = re.compile(r'gene_biotype "([^"]+)"')
TRANSCRIPT_BIOTYPE_RE = re.compile(r'transcript_biotype "([^"]+)"')


def eprint(message: str) -> None:
    print(message, file=sys.stderr, flush=True)


def chrom_aliases(chrom: str) -> List[str]:
    if not chrom:
        return []

    aliases = [chrom]
    if chrom.startswith("chr"):
        base = chrom[3:]
        aliases.append(base)
        if base == "M":
            aliases.append("MT")
        elif base == "MT":
            aliases.append("M")
    else:
        aliases.append(f"chr{chrom}")
        if chrom == "MT":
            aliases.extend(["chrM", "M"])
        elif chrom == "M":
            aliases.extend(["chrM", "MT"])

    out = []
    seen = set()
    for alias in aliases:
        if alias and alias not in seen:
            seen.add(alias)
            out.append(alias)
    return out


@dataclass
class SimpleSeqRecord:
    seq: str
    description: str = ""


def build_reference_chrom_map(gtf_path: Path) -> Dict[str, str]:
    ref_chroms: set[str] = set()
    with gtf_path.open() as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t", 1)
            if parts and parts[0]:
                ref_chroms.add(parts[0])

    chrom_map = {chrom: chrom for chrom in ref_chroms}
    for chrom in ref_chroms:
        for alias in chrom_aliases(chrom):
            chrom_map.setdefault(alias, chrom)
    return chrom_map


def write_profile(profile_file: Path, rows: Iterable[Tuple[str, float, str]]) -> None:
    with profile_file.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["stage", "seconds", "details"])
        for row in rows:
            writer.writerow(row)


def load_orf_to_study_and_sequences(metadata_file: Path) -> Tuple[Dict[str, str], Dict[str, str]]:
    try:
        csv.field_size_limit(sys.maxsize)
    except OverflowError:
        csv.field_size_limit(2147483647)
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


def write_bed6(
    bed12_file: Path,
    bed6_file: Path,
    orf_to_study: Dict[str, str],
    chrom_map: Dict[str, str],
) -> Tuple[int, int]:
    count = 0
    remapped = 0
    with bed12_file.open() as src, bed6_file.open("w") as dst:
        for line in src:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue
            chrom, start, end, name, _, strand = parts[:6]
            mapped_chrom = chrom_map.get(chrom, chrom)
            if mapped_chrom != chrom:
                remapped += 1
            dst.write(
                f"{mapped_chrom}\t{start}\t{end}\t{name}\t{orf_to_study.get(name, 'unified')}\t{strand}\n"
            )
            count += 1
    return count, remapped


def translate_nt(nt_seq: str) -> str:
    nt_seq = nt_seq.upper().replace("U", "T")
    aa: List[str] = []
    for i in range(0, len(nt_seq) - 2, 3):
        aa.append(CODON_TABLE.get(nt_seq[i : i + 3], "X"))
    return "".join(aa)


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


def parse_reference_paths(ensembl_dir: Path) -> Dict[str, Path]:
    return {
        "TRANSCRIPTOME_FASTA": ensembl_dir / "TRANSCRIPTOME_FASTA",
        "SORTED_TRANSCRIPTOME_GTF": ensembl_dir / "SORTED_TRANSCRIPTOME_GTF",
        "PROTEOME_FASTA": ensembl_dir / "PROTEOME_FASTA",
        "TRANSCRIPT_SUPPORT": ensembl_dir / "TRANSCRIPT_SUPPORT",
        "PSITES_BED": ensembl_dir / "PSITES_BED",
    }


def load_bed6_groups(bed6_file: Path) -> Tuple[List[str], Dict[str, List[str]], Dict[str, str]]:
    chrom_order: List[str] = []
    groups: Dict[str, List[str]] = {}
    orf_to_chrom: Dict[str, str] = {}
    with bed6_file.open() as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 4:
                continue
            chrom = parts[0]
            if chrom not in groups:
                chrom_order.append(chrom)
                groups[chrom] = []
            groups[chrom].append(line)
            if len(parts) >= 5:
                orf_to_chrom[parts[3] + "--" + parts[4]] = chrom
    return chrom_order, groups, orf_to_chrom


def write_feature_filtered_gtf(
    gtf_path: Path,
    output_path: Path,
    needed_chroms: set[str],
) -> int:
    kept = 0
    keep_features = {"transcript", "exon", "CDS"}
    with gtf_path.open() as src, output_path.open("w") as dst:
        for line in src:
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t", 8)
            if len(parts) < 9:
                continue
            if parts[0] not in needed_chroms or parts[2] not in keep_features:
                continue
            dst.write(line)
            kept += 1
    return kept


def assign_bins(
    chrom_order: Sequence[str],
    groups: Dict[str, List[str]],
    max_bins: int,
) -> List[List[str]]:
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


def build_tx_to_protein(gtf_path: Path) -> Dict[str, str]:
    tx_to_protein: Dict[str, str] = {}
    with gtf_path.open() as handle:
        for line in handle:
            if "\tCDS\t" not in line:
                continue
            tx_match = TRANSCRIPT_ID_RE.search(line)
            prot_match = PROTEIN_ID_RE.search(line)
            if tx_match and prot_match:
                tx_to_protein[tx_match.group(1)] = prot_match.group(1)
    return tx_to_protein


def load_fasta_records(
    fasta_path: Path,
    needed_ids: set[str],
) -> Dict[str, SimpleSeqRecord]:
    records: Dict[str, SimpleSeqRecord] = {}
    if not needed_ids:
        return records
    with fasta_path.open() as handle:
        for title, sequence in SimpleFastaParser(handle):
            seq_id = title.split(None, 1)[0]
            if seq_id in needed_ids:
                records[seq_id] = SimpleSeqRecord(sequence, title)
    return records


def load_protein_sequences(
    fasta_path: Path,
    needed_ids: set[str],
) -> Dict[str, str]:
    sequences: Dict[str, str] = {}
    if not needed_ids:
        return sequences
    with fasta_path.open() as handle:
        for title, sequence in SimpleFastaParser(handle):
            seq_id = title.split(None, 1)[0]
            if seq_id in needed_ids:
                sequences[seq_id] = sequence.replace("X", "")
    return sequences


def load_support_filtered(
    t_support_path: Path,
    needed_tx_ids: set[str] | None = None,
) -> Tuple[Dict[str, str], Dict[str, str]]:
    appris: Dict[str, str] = {}
    supp: Dict[str, str] = {}
    with t_support_path.open() as handle:
        for line in handle:
            if line.startswith("Transcript"):
                continue
            tx_id = line.split("\t", 1)[0]
            if needed_tx_ids is not None and tx_id not in needed_tx_ids:
                continue
            if "principal" in line or "alternative" in line:
                appris[tx_id] = line.split("\t")[2].rstrip("\n")
            else:
                appris[tx_id] = "tslNA"

            if "ts" in line:
                supp[tx_id] = line.split("\t")[1].split(" (")[0].rstrip("\n")
            else:
                supp[tx_id] = "tslNA"
    return appris, supp


def parse_gtf_filtered(
    gtf_path: Path,
    needed_tx_ids: set[str],
) -> Tuple[Dict[str, object], Dict[str, str]]:
    trans: Dict[str, object] = {}
    tx_to_protein: Dict[str, str] = {}
    with gtf_path.open() as handle:
        for line in handle:
            if "\t" not in line or line.startswith("#"):
                continue
            tx_match = TRANSCRIPT_ID_RE.search(line)
            if tx_match is None:
                continue
            tx_id = tx_match.group(1)
            if tx_id not in needed_tx_ids:
                continue
            if "\texon\t" in line:
                g_name = GENE_ID_RE.search(line)
                gene_id = g_name.group(1) if g_name else "unknown"
                if "gene_biotype" in line:
                    biot = GENE_BIOTYPE_RE.search(line)
                    biotype = biot.group(1) if biot else "unknown"
                else:
                    biotype = "unknown"
                fields = line.split("\t")
                trans.setdefault(
                    tx_id,
                    functions.trans_object(
                        fields[0],
                        gene_id,
                        fields[6],
                        [],
                        [],
                        biotype,
                    ),
                )
                trans[tx_id].start.append(int(fields[3]))
                trans[tx_id].end.append(int(fields[4]))
            elif "\tCDS\t" in line:
                prot_match = PROTEIN_ID_RE.search(line)
                if prot_match:
                    tx_to_protein[tx_id] = prot_match.group(1)

    for tx_id in trans:
        trans[tx_id].start.sort()
        trans[tx_id].end.sort()
    return trans, tx_to_protein


def build_cum_offsets(exon_starts: Sequence[int], exon_ends: Sequence[int]) -> List[int]:
    cum = []
    total = 0
    for exon_start, exon_end in zip(exon_starts, exon_ends):
        cum.append(total)
        total += exon_end - exon_start + 1
    return cum


def find_cds_in_frames(
    frames: Sequence[str],
    protein_seq: str,
) -> Tuple[int | None, int | None, int | None]:
    if not protein_seq:
        return None, None, None
    hits = [frame.find(protein_seq) for frame in frames]
    valid_hits = [(hit, idx) for idx, hit in enumerate(hits) if hit >= 0]
    if not valid_hits:
        return None, None, None
    best_hit, frame = min(valid_hits, key=lambda item: item[0] * 3 + item[1])
    return best_hit * 3 + frame, frame, len(protein_seq)


def build_reference_index(
    cache_dir: Path,
    transcriptome_records: Dict[str, SimpleSeqRecord],
    protein_sequences: Dict[str, str],
    gtf: Dict[str, object],
    appris: Dict[str, str],
    supp: Dict[str, str],
    tx_to_protein: Dict[str, str],
) -> Dict[str, str]:
    cache_dir.mkdir(parents=True, exist_ok=True)
    shards: Dict[str, Dict[str, dict]] = {}
    for tx_id, tx_obj in gtf.items():
        record = transcriptome_records.get(tx_id)
        if record is None:
            continue
        tx_seq = record.seq
        protein_id = tx_to_protein.get(tx_id)
        protein_seq = protein_sequences.get(protein_id, "") if protein_id else ""
        chrom = tx_obj.chrm
        shards.setdefault(chrom, {})[tx_id] = {
            "chrom": chrom,
            "gene": tx_obj.gene,
            "strand": tx_obj.strand,
            "gene_biotype": tx_obj.biotype,
            "exon_starts": list(tx_obj.start),
            "exon_ends": list(tx_obj.end),
            "cum_offsets": build_cum_offsets(tx_obj.start, tx_obj.end),
            "tx_seq": tx_seq,
            "protein_seq": protein_seq,
            "tx_len": len(tx_seq),
            "appris": appris.get(tx_id, "none"),
            "supp": supp.get(tx_id, "0"),
        }

    manifest: Dict[str, str] = {}
    for chrom, shard in shards.items():
        shard_path = cache_dir / f"{chrom}.pkl"
        with shard_path.open("wb") as handle:
            pickle.dump(shard, handle, protocol=pickle.HIGHEST_PROTOCOL)
        manifest[chrom] = str(shard_path)
    return manifest


def run_intersect(orfs_bed_file: Path, transcriptome_gtf_file: Path, output_path: Path) -> None:
    with output_path.open("w") as out_handle:
        subprocess.run(
            [
                "intersectBed",
                "-s",
                "-a",
                str(orfs_bed_file),
                "-b",
                str(transcriptome_gtf_file),
                "-wa",
                "-wb",
            ],
            check=True,
            stdout=out_handle,
        )


def parse_overlap_file(
    overlap_path: Path,
) -> Tuple[Dict[str, List[List[str]]], Dict[str, List[str]], List[str], set[str]]:
    overlaps: Dict[str, List[List[str]]] = {}
    other_overlaps: Dict[str, List[str]] = {}
    overlap_seen: Dict[str, set] = {}
    total_studies = set()
    needed_tx_ids: set[str] = set()
    with overlap_path.open() as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 15:
                continue
            feature = parts[8]
            attrs = parts[14]
            name = parts[3] + "--" + parts[4]
            other_overlaps.setdefault(name, ["0", "0", "0"])
            total_studies.add(parts[4])

            if feature == "exon" and "pseudogene" in attrs:
                other_overlaps[name][0] = "2" if "unitary" in attrs else "1"
            if feature == "CDS":
                other_overlaps[name][1] = "1"

            if feature != "transcript":
                continue

            tx_match = TRANSCRIPT_ID_RE.search(attrs)
            gene_match = GENE_ID_RE.search(attrs)
            if tx_match is None or gene_match is None:
                continue
            gene_name_match = GENE_NAME_RE.search(attrs)
            gene_biotype_match = GENE_BIOTYPE_RE.search(attrs)
            transcript_biotype_match = TRANSCRIPT_BIOTYPE_RE.search(attrs)
            tx_id = tx_match.group(1)
            needed_tx_ids.add(tx_id)
            gene_id = gene_match.group(1)
            gene_name = gene_name_match.group(1) if gene_name_match else gene_id
            gene_biotype = gene_biotype_match.group(1) if gene_biotype_match else "unknown"
            transcript_biotype = (
                transcript_biotype_match.group(1) if transcript_biotype_match else "unknown"
            )
            overlap_key = (tx_id, gene_id, transcript_biotype, gene_biotype, gene_name)
            overlaps.setdefault(name, [])
            overlap_seen.setdefault(name, set())
            if overlap_key not in overlap_seen[name]:
                overlaps[name].append(
                    [tx_id, gene_id, transcript_biotype, gene_biotype, gene_name]
                )
                overlap_seen[name].add(overlap_key)
    return overlaps, other_overlaps, sorted(total_studies), needed_tx_ids


def load_shards(shard_paths: Sequence[str]) -> Dict[str, dict]:
    tx_index: Dict[str, dict] = {}
    for shard_path in shard_paths:
        with open(shard_path, "rb") as handle:
            tx_index.update(pickle.load(handle))
    return tx_index


def tx_to_genome(tx_data: dict, tx_pos: int) -> int | None:
    for idx in range(len(tx_data["cum_offsets"]) - 1, -1, -1):
        if tx_data["cum_offsets"][idx] <= tx_pos:
            return tx_data["exon_starts"][idx] + (tx_pos - tx_data["cum_offsets"][idx])
    return None


def resolve_tx_runtime(tx_data: dict, runtime_cache: Dict[str, dict], tx_id: str) -> dict:
    cached = runtime_cache.get(tx_id)
    if cached is not None:
        return cached
    tx_seq = tx_data["tx_seq"]
    frames = (
        translate_nt(tx_seq),
        translate_nt(tx_seq[1:]),
        translate_nt(tx_seq[2:]),
    )
    cds_start, cds_frame, cds_len = find_cds_in_frames(frames, tx_data["protein_seq"])
    cached = {
        "frames": frames,
        "cds_start": cds_start,
        "cds_frame": cds_frame,
        "cds_len": cds_len,
    }
    runtime_cache[tx_id] = cached
    return cached


def process_orf_chunk(chunk_payload: dict) -> Tuple[dict, dict, dict, dict, List[str], dict]:
    tx_index = load_shards(chunk_payload["shard_paths"])
    tx_runtime_cache: Dict[str, dict] = {}
    candidates: Dict[str, List[list]] = {}
    trans_orfs: Dict[str, List[list]] = {}
    coord_psites: Dict[str, set] = {}
    second_names: Dict[str, str] = {}
    atgstop_lines: List[str] = []
    stats = {"orfs": 0, "pairs": 0}

    for orf_id, orf_data in chunk_payload["orfs"].items():
        stats["orfs"] += 1
        orf_seq = orf_data["orf_seq"]
        if not orf_seq:
            continue
        aa_len = len(orf_seq.replace("*", ""))
        if aa_len < chunk_payload["len_cutoff"] or aa_len > chunk_payload["max_len_cutoff"]:
            continue
        try:
            second_names[orf_id] = orf_data["orf_desc"].split()[-1].split("--")[0]
        except Exception:
            second_names[orf_id] = orf_id.split("--")[0]
        orf_len_nt = len(orf_seq) * 3

        for trans in orf_data["overlaps"]:
            stats["pairs"] += 1
            tx_id = trans[0]
            tx_data = tx_index.get(tx_id)
            if tx_data is None:
                continue
            tx_runtime = resolve_tx_runtime(tx_data, tx_runtime_cache, tx_id)
            gene = trans[1]
            gene_name = trans[4]
            cat2 = "protein_coding" if trans[3] == "protein_coding" else "non-coding"

            frames = tx_runtime["frames"]
            f1 = frames[0].find(orf_seq)
            f2 = frames[1].find(orf_seq)
            f3 = frames[2].find(orf_seq)
            best_orf = max(f1, f2, f3)
            if best_orf < 0:
                continue
            fi = [f1, f2, f3].index(best_orf)
            f = best_orf * 3 + fi
            intersection = set()

            if tx_runtime["cds_start"] is not None:
                c = tx_runtime["cds_start"]
                ci = tx_runtime["cds_frame"]
                cc = tx_runtime["cds_len"]
                fr = range(f, f + (len(orf_seq) * 3))
                cr = range(c, c + (cc * 3))
                intersection = set(fr).intersection(cr)
                if len(intersection) == 0:
                    cat = "dORF" if f > c else "uORF"
                else:
                    if fi == ci:
                        cat = "CDS"
                    else:
                        if f > c:
                            if f + (len(orf_seq) * 3) > c + (cc * 3):
                                cat = "doORF"
                            else:
                                cat = "intORF"
                        else:
                            cat = "uoORF"
            else:
                cat = "lncRNA"

            candidates.setdefault(orf_id, []).append(
                [tx_id, gene, gene_name, cat, cat2, fi, f, orf_seq, len(intersection)]
            )
            trans_orfs.setdefault(tx_id, []).append([orf_id, fi, f, orf_seq])

            genomic_f = f
            if tx_data["strand"] == "-":
                genomic_f = tx_data["tx_len"] - (f + orf_len_nt)
            atg_gp = tx_to_genome(tx_data, genomic_f)
            if atg_gp is not None:
                atgstop_lines.append(
                    tx_data["chrom"]
                    + "\t"
                    + str(atg_gp)
                    + "\t"
                    + str(atg_gp)
                    + "\t"
                    + orf_id
                    + "\tboundaries\t"
                    + tx_data["strand"]
                    + "\n"
                )
            stop_gp = tx_to_genome(tx_data, genomic_f + orf_len_nt - 1)
            if stop_gp is not None:
                atgstop_lines.append(
                    tx_data["chrom"]
                    + "\t"
                    + str(stop_gp)
                    + "\t"
                    + str(stop_gp)
                    + "\t"
                    + orf_id
                    + "\tboundaries\t"
                    + tx_data["strand"]
                    + "\n"
                )
            coord_psites.setdefault(orf_id, set())
            if tx_data["strand"] == "+":
                psite_offsets = range(genomic_f + 2, genomic_f + orf_len_nt - 1, 3)
            else:
                psite_offsets = range(genomic_f, genomic_f + orf_len_nt, 3)
            for tx_pos in psite_offsets:
                gp = tx_to_genome(tx_data, tx_pos)
                if gp is not None:
                    coord_psites[orf_id].add(gp)

    return candidates, trans_orfs, coord_psites, second_names, atgstop_lines, stats


def build_worker_payloads(
    overlaps: Dict[str, List[List[str]]],
    orf_to_chrom: Dict[str, str],
    chrom_bins: Sequence[Sequence[str]],
    shard_manifest: Dict[str, str],
    orfs_fa,
    len_cutoff: int,
    max_len_cutoff: int,
) -> List[dict]:
    bin_index: Dict[str, int] = {}
    for idx, chroms in enumerate(chrom_bins, start=1):
        for chrom in chroms:
            bin_index[chrom] = idx
    payloads: Dict[int, dict] = {}
    for orf_id, overlap_rows in overlaps.items():
        chrom = orf_to_chrom.get(orf_id)
        if chrom is None:
            continue
        chunk_id = bin_index.get(chrom, 1)
        payload = payloads.setdefault(
            chunk_id,
            {
                "chunk_id": chunk_id,
                "chroms": list(chrom_bins[chunk_id - 1]),
                "shard_paths": [shard_manifest[c] for c in chrom_bins[chunk_id - 1] if c in shard_manifest],
                "orfs": {},
                "len_cutoff": len_cutoff,
                "max_len_cutoff": max_len_cutoff,
            },
        )
        if orf_id not in orfs_fa:
            continue
        payload["orfs"][orf_id] = {
            "orf_seq": str(orfs_fa[orf_id].seq),
            "orf_desc": str(orfs_fa[orf_id].description),
            "overlaps": overlap_rows,
        }
    return sorted(payloads.values(), key=lambda item: item["chunk_id"])


def merge_results(
    results: Sequence[Tuple[dict, dict, dict, dict, List[str], dict]],
    overlaps: Dict[str, List[List[str]]],
) -> Tuple[dict, dict, dict, dict, List[str]]:
    candidates: Dict[str, List[list]] = {}
    trans_orfs: Dict[str, List[list]] = {}
    coord_psites: Dict[str, set] = {}
    second_names: Dict[str, str] = {}
    atgstop_lines: List[str] = []

    for cand_local, trans_local, psites_local, names_local, atg_lines, _ in results:
        for orf, rows in cand_local.items():
            candidates.setdefault(orf, []).extend(rows)
        for tx_id, rows in trans_local.items():
            trans_orfs.setdefault(tx_id, []).extend(rows)
        for orf, coords in psites_local.items():
            coord_psites.setdefault(orf, set()).update(coords)
        second_names.update(names_local)
        atgstop_lines.extend(atg_lines)

    orf_order_map = {orf_id: idx for idx, orf_id in enumerate(overlaps.keys())}
    trans_orfs = functions._canonicalize_trans_orfs(trans_orfs, orf_order_map)
    coord_psites = {orf: sorted(coords) for orf, coords in coord_psites.items()}
    return candidates, trans_orfs, coord_psites, second_names, atgstop_lines


def main() -> None:
    parser = argparse.ArgumentParser(description="Indexed fast GENCODE classifier entrypoint")
    parser.add_argument("--input_prefix", required=True, help="Unified ORF prefix")
    parser.add_argument("--ensembl_dir", required=True, help="Prepared Ensembl directory")
    parser.add_argument("--output_dir", required=True, help="Output directory")
    parser.add_argument("--cpus", type=int, default=1, help="Number of worker processes")
    parser.add_argument("--len_cutoff", type=int, default=16, help="Minimum ORF length threshold")
    parser.add_argument("--max_len_cutoff", type=int, default=999999999999, help="Maximum ORF length threshold")
    parser.add_argument("--collapse_cutoff", type=float, default=0.9, help="Variant collapse threshold")
    parser.add_argument("--collapse_method", default="longest_string", choices=["longest_string", "psite_overlap"], help="Variant collapse method")
    args = parser.parse_args()

    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    input_prefix = Path(args.input_prefix)
    metadata_file = Path(f"{input_prefix}.metadata.tsv")
    bed12_file = Path(f"{input_prefix}.bed")
    fasta_file = Path(f"{input_prefix}.orfs.fa")
    bed6_file = Path(f"{input_prefix}.orfs.bed6")
    output_prefix = output_dir / "gencode_results"
    profile_file = output_dir / "gencode_indexed_fast.profile.tsv"
    cache_dir = output_dir / "indexed_fast_cache"
    overlap_file = output_dir / "indexed_fast.orfs_to_gtf.ov"
    filtered_gtf_file = output_dir / "indexed_fast.filtered.gtf"
    temp_seed = "".join(random.choice(string.ascii_letters) for _ in range(10))
    psites_overlap_seed = temp_seed
    temp_dir = output_dir / "tmp"
    temp_dir.mkdir(parents=True, exist_ok=True)
    atgstop_file = temp_dir / f"{psites_overlap_seed}atg_to_stop.bed"

    paths = parse_reference_paths(Path(args.ensembl_dir).resolve())
    for required in (metadata_file, bed12_file, paths["TRANSCRIPTOME_FASTA"], paths["SORTED_TRANSCRIPTOME_GTF"], paths["PROTEOME_FASTA"], paths["TRANSCRIPT_SUPPORT"], paths["PSITES_BED"]):
        if not Path(required).exists():
            raise FileNotFoundError(f"Required input not found: {required}")
    chrom_map = build_reference_chrom_map(paths["SORTED_TRANSCRIPTOME_GTF"])

    timings: List[Tuple[str, float, str]] = []

    t0 = time.perf_counter()
    orf_to_study, orf_to_nt = load_orf_to_study_and_sequences(metadata_file)
    timings.append(("load_metadata", time.perf_counter() - t0, f"orfs={len(orf_to_study)}"))

    t0 = time.perf_counter()
    n_bed, remapped = write_bed6(bed12_file, bed6_file, orf_to_study, chrom_map)
    timings.append(("write_bed6", time.perf_counter() - t0, f"records={n_bed} remapped={remapped}"))
    if remapped:
        eprint(f"Remapped {remapped} BED records to reference chromosome names")

    if not fasta_file.exists():
        t0 = time.perf_counter()
        n_fa = write_protein_fasta(fasta_file, orf_to_study, orf_to_nt)
        timings.append(("write_protein_fasta", time.perf_counter() - t0, f"records={n_fa}"))

    t0 = time.perf_counter()
    chrom_order, bed_groups, orf_to_chrom = load_bed6_groups(bed6_file)
    write_feature_filtered_gtf(
        paths["SORTED_TRANSCRIPTOME_GTF"],
        filtered_gtf_file,
        set(chrom_order),
    )
    timings.append(("filter_gtf_for_overlap", time.perf_counter() - t0, f"chromosomes={len(chrom_order)}"))

    t0 = time.perf_counter()
    run_intersect(bed6_file, filtered_gtf_file, overlap_file)
    overlaps, other_overlaps, total_studies, needed_tx_ids = parse_overlap_file(overlap_file)
    timings.append(("intersect_and_parse", time.perf_counter() - t0, f"orfs={len(overlaps)} studies={len(total_studies)} transcripts={len(needed_tx_ids)}"))

    t0 = time.perf_counter()
    if needed_tx_ids:
        gtf, tx_to_protein = parse_gtf_filtered(filtered_gtf_file, needed_tx_ids)
        appris, supp = load_support_filtered(paths["TRANSCRIPT_SUPPORT"], needed_tx_ids)
        transcriptome_fa = load_fasta_records(paths["TRANSCRIPTOME_FASTA"], needed_tx_ids)
        proteome_fa = load_protein_sequences(
            paths["PROTEOME_FASTA"],
            set(tx_to_protein.values()),
        )
        orfs_fa = load_fasta_records(fasta_file, set(overlaps.keys()))
        shard_manifest = build_reference_index(
            cache_dir=cache_dir,
            transcriptome_records=transcriptome_fa,
            protein_sequences=proteome_fa,
            gtf=gtf,
            appris=appris,
            supp=supp,
            tx_to_protein=tx_to_protein,
        )
    else:
        transcriptome_fa = {}
        proteome_fa = {}
        orfs_fa = {}
        gtf = {}
        appris = {}
        supp = {}
        shard_manifest = {}
    timings.append(("build_reference_index", time.perf_counter() - t0, f"chromosomes={len(shard_manifest)} transcripts={len(gtf)}"))

    t0 = time.perf_counter()
    chrom_bins = assign_bins(chrom_order, bed_groups, max(1, args.cpus))
    worker_payloads = build_worker_payloads(
        overlaps=overlaps,
        orf_to_chrom=orf_to_chrom,
        chrom_bins=chrom_bins,
        shard_manifest=shard_manifest,
        orfs_fa=orfs_fa,
        len_cutoff=args.len_cutoff,
        max_len_cutoff=args.max_len_cutoff,
    )
    timings.append(("prepare_worker_payloads", time.perf_counter() - t0, f"chunks={len(worker_payloads)}"))

    t0 = time.perf_counter()
    results = []
    with ProcessPoolExecutor(max_workers=max(1, min(args.cpus, len(worker_payloads) or 1))) as executor:
        futures = {
            executor.submit(process_orf_chunk, payload): payload["chunk_id"]
            for payload in worker_payloads
        }
        for future in as_completed(futures):
            chunk_id = futures[future]
            result = future.result()
            results.append(result)
            stats = result[-1]
            eprint(
                f"Completed indexed_fast chunk {chunk_id} "
                f"(orfs={stats['orfs']} pairs={stats['pairs']})"
            )
    timings.append(("parallel_orf_classification", time.perf_counter() - t0, f"chunks={len(worker_payloads)}"))

    t0 = time.perf_counter()
    candidates, trans_orfs, coord_psites, second_names, atgstop_lines = merge_results(results, overlaps)
    with atgstop_file.open("w") as handle:
        handle.writelines(atgstop_lines)
    timings.append(("merge_worker_results", time.perf_counter() - t0, f"candidates={len(candidates)}"))

    t0 = time.perf_counter()
    shutil_overlap = temp_dir / f"{psites_overlap_seed}atg_to_stop.ov"
    exc, variants, variants_names, datasets = functions.exclude_variants(
        trans_orfs,
        args.collapse_cutoff,
        candidates,
        args.collapse_method,
        coord_psites,
        str(output_dir),
        psites_overlap_seed,
    )
    timings.append(("exclude_variants", time.perf_counter() - t0, f"exc={len(exc)}"))
    if atgstop_file.exists():
        atgstop_file.unlink()
    if shutil_overlap.exists():
        shutil_overlap.unlink()

    t0 = time.perf_counter()
    functions.write_output(
        str(fasta_file),
        str(bed6_file),
        candidates,
        exc,
        variants,
        variants_names,
        datasets,
        appris,
        supp,
        gtf,
        transcriptome_fa,
        second_names,
        args.len_cutoff,
        args.max_len_cutoff,
        args.collapse_cutoff,
        total_studies,
        other_overlaps,
        str(paths["PSITES_BED"]),
        str(output_dir),
        str(output_prefix),
        args.collapse_method,
        psites_overlap_seed,
        "none",
        "none",
        "no",
    )
    timings.append(("write_output", time.perf_counter() - t0, ""))

    write_profile(profile_file, timings)


if __name__ == "__main__":
    main()
