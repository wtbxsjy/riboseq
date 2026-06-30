#!/usr/bin/env python
import sys
import string
import subprocess
import os
import random
import string
import time
from concurrent.futures import ThreadPoolExecutor
from Bio import SeqIO
from Bio.Seq import Seq
from datetime import datetime
from optparse import OptionParser

__author__ = "Jorge Ruiz-Orera"
__contributor__ = "..."
__copyright__ = ""
__credits__ = []
__license__ = ""
__version__ = "1.1.0"
__maintainer__ = "Jorge Ruiz-Orera"
__email__ = "jorruior@gmail.com"


### FUNCTIONS
class trans_object:
    def __init__(self, chrm, gene, strand, start, end, biotype):
        self.chrm = chrm
        self.gene = gene
        self.strand = strand
        self.start = start
        self.end = end
        self.biotype = biotype


class ProgressTicker:
    """Lightweight progress reporter for long-running terminal loops."""

    def __init__(self, label, total, min_interval=5.0, min_step=0.02):
        self.label = label
        self.total = max(int(total), 0)
        self.min_interval = float(min_interval)
        self.min_step = float(min_step)
        self.started_at = time.time()
        self.last_report_at = self.started_at
        self.last_fraction = -1.0

    def update(self, current, extra=""):
        if self.total <= 0:
            return
        current = min(max(int(current), 0), self.total)
        fraction = float(current) / float(self.total)
        now = time.time()
        is_done = current >= self.total
        should_report = (
            current == 1
            or is_done
            or (fraction - self.last_fraction) >= self.min_step
            or (now - self.last_report_at) >= self.min_interval
        )
        if not should_report:
            return

        elapsed = now - self.started_at
        message = (
            "[progress] "
            + self.label
            + ": "
            + str(current)
            + "/"
            + str(self.total)
            + " ("
            + "{:.1f}".format(fraction * 100.0)
            + "%, elapsed "
            + "{:.1f}".format(elapsed)
            + "s)"
        )
        if extra:
            message += " " + str(extra)
        print(message, flush=True)
        self.last_report_at = now
        self.last_fraction = fraction


def lcs(S, T):
    """Return longest common substring.

    Uses difflib.SequenceMatcher for fast C-level string comparison
    instead of O(|S|×|T|) pure-Python dynamic programming.
    Returns a set with a single string (the longest common substring)
    to preserve the original interface.
    """
    import difflib
    m = difflib.SequenceMatcher(None, S, T)
    block = m.find_longest_match(0, len(S), 0, len(T))
    if block.size == 0:
        return set()
    return {S[block.a : block.a + block.size]}


def get_index_positions(list_of_elems, element):
    """Returns the indexes of all occurrences of give element in
    the list- listOfElements"""
    index_pos_list = []
    index_pos = 0
    while True:
        try:
            # Search for item in list from indexPos to the end of list
            index_pos = list_of_elems.index(element, index_pos)
            # Add the index position in list
            index_pos_list.append(index_pos)
            index_pos += 1
        except ValueError as e:
            break
    return index_pos_list


def parse_gtf(gtf, field):
    """Read a gtf and create a dict with sorted transcript coordinates, chrm, strand, and gene"""
    trans = {}
    for line in open(gtf):
        if not "\t" + field + "\t" in line:
            continue
        t_name = line.split('transcript_id "')[1].split('"')[0]
        g_name = line.split('gene_id "')[1].split('"')[0]

        if "gene_biotype" in line:
            biot = line.split('gene_biotype "')[1].split('"')[0]
        else:
            biot = "unknown"

        trans.setdefault(
            t_name,
            trans_object(
                line.split("\t")[0], g_name, line.split("\t")[6], [], [], biot
            ),
        )
        trans[t_name].start.append(int(line.split("\t")[3]))
        trans[t_name].end.append(int(line.split("\t")[4]))

    [trans[x].start.sort() for x in trans]
    [trans[x].end.sort() for x in trans]

    return trans


def load_fasta(orfs_fa_file, transcriptome_fa_file, proteome_fa_file):
    orfs_fa = SeqIO.index(orfs_fa_file, "fasta")
    transcriptome_fa = SeqIO.index(transcriptome_fa_file, "fasta")
    proteome_fa = SeqIO.index(proteome_fa_file, "fasta")
    return orfs_fa, transcriptome_fa, proteome_fa


def read_support(t_support):
    appris = {}
    supp = {}
    for line in open(t_support):
        if line.startswith("Transcript"):
            continue
        if "principal" in line or "alternative" in line:
            appris[line.split("\t")[0]] = line.split("\t")[2].rstrip("\n")
        else:
            appris[line.split("\t")[0]] = "tslNA"

        if "ts" in line:
            supp[line.split("\t")[0]] = line.split("\t")[1].split(" (")[0].rstrip("\n")
        else:
            supp[line.split("\t")[0]] = "tslNA"

    return appris, supp


def project_tx_range_to_genome(exon_starts, exon_ends, tx_start, tx_len):
    """Project a transcript-space interval onto genomic exon blocks."""
    if tx_start < 0 or tx_len <= 0:
        return []

    blocks = []
    remaining_start = tx_start
    remaining_len = tx_len

    for exon_start, exon_end in zip(exon_starts, exon_ends):
        exon_len = exon_end - exon_start + 1
        if remaining_start >= exon_len:
            remaining_start -= exon_len
            continue

        seg_start = exon_start + remaining_start
        take_len = min(remaining_len, exon_len - remaining_start)
        seg_end = seg_start + take_len - 1
        blocks.append((seg_start, seg_end))

        remaining_len -= take_len
        if remaining_len <= 0:
            break
        remaining_start = 0

    return blocks


def _resolve_fast_cpus():
    """Return requested CPU count for fast mode, defaulting to sequential."""
    raw = os.environ.get("GENCODE_FAST_CPUS", "").strip()
    if not raw:
        return 1
    try:
        cpus = int(raw)
    except ValueError:
        return 1
    return max(1, cpus)


def _chunk_items(items, n_chunks):
    if n_chunks <= 1 or len(items) <= 1:
        return [items]
    chunk_size = max(1, (len(items) + n_chunks - 1) // n_chunks)
    return [items[i : i + chunk_size] for i in range(0, len(items), chunk_size)]


def _canonicalize_trans_orfs(trans_orfs, orf_order_map):
    """Rebuild trans_orfs with deterministic ordering matching ORF input order."""
    tx_with_order = []
    for tx_id, rows in trans_orfs.items():
        sorted_rows = sorted(
            rows,
            key=lambda row: (orf_order_map.get(row[0], 10**18), row[2], row[1], row[0]),
        )
        first_idx = orf_order_map.get(sorted_rows[0][0], 10**18) if sorted_rows else 10**18
        tx_with_order.append((first_idx, tx_id, sorted_rows))

    tx_with_order.sort(key=lambda item: (item[0], item[1]))
    ordered = {}
    for _, tx_id, rows in tx_with_order:
        ordered[tx_id] = rows
    return ordered


def _process_orf_chunk(
    orf_items,
    overlaps_cds,
    orf_seq_map,
    orf_desc_map,
    transcript_seq_map,
    protein_seq_map,
    gtf,
    len_cutoff,
    max_len_cutoff,
    progress=None,
):
    candidates = {}
    trans_orfs = {}
    coord_psites = {}
    second_names = {}
    atgstop_lines = []
    tx_trans_cache = {}

    def _get_tx_frames(tx_id):
        if tx_id not in tx_trans_cache:
            seq = transcript_seq_map[tx_id]
            tx_trans_cache[tx_id] = (
                str(Seq(seq).translate(cds=False)),
                str(Seq(seq[1:]).translate(cds=False)),
                str(Seq(seq[2:]).translate(cds=False)),
            )
        return tx_trans_cache[tx_id]

    for idx, (orf, overlap_rows) in enumerate(orf_items, start=1):
        if progress is not None:
            progress.update(idx, "orf=" + str(orf))
        cat2 = "non-coding"
        orf_seq = orf_seq_map.get(orf)
        if orf_seq is None:
            continue
        try:
            second_names[orf] = orf_desc_map[orf].split()[-1].split("--")[0]
        except Exception:
            second_names[orf] = orf.split("--")[0]
        aa_len = len(orf_seq.replace("*", ""))
        if aa_len < len_cutoff or aa_len > max_len_cutoff:
            continue

        for trans in overlap_rows:
            tx_id = trans[0]
            if tx_id not in transcript_seq_map:
                continue
            gene = trans[1]
            genename = trans[4]
            frames = _get_tx_frames(tx_id)
            f1 = frames[0].find(orf_seq)
            f2 = frames[1].find(orf_seq)
            f3 = frames[2].find(orf_seq)
            fi = [f1, f2, f3].index(max([f1, f2, f3]))
            f = max(f1, f2, f3) * 3 + fi
            if f < 0:
                continue
            if trans[3] == "protein_coding":
                cat2 = "protein_coding"

            intersection = []
            if tx_id in overlaps_cds:
                prot_seq = protein_seq_map.get(tx_id)
                if prot_seq is None:
                    continue
                c1 = frames[0].find(prot_seq)
                c2 = frames[1].find(prot_seq)
                c3 = frames[2].find(prot_seq)
                cc = len(prot_seq)
                ci = [c1, c2, c3].index(max([c1, c2, c3]))
                c = max(c1, c2, c3) * 3 + ci
                if c == -3:
                    continue

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

            candidates.setdefault(orf, [])
            candidates[orf].append([tx_id, gene, genename, cat, cat2, fi, f, orf_seq, len(intersection)])

            trans_orfs.setdefault(tx_id, [])
            trans_orfs[tx_id].append([orf, fi, f, orf_seq])

            strand = gtf[tx_id].strand
            genomic_f = f
            if strand == "-":
                genomic_f = len(transcript_seq_map[tx_id]) - (f + (len(orf_seq) * 3))
            orf_len_nt = len(orf_seq) * 3
            exon_starts = gtf[tx_id].start
            exon_ends = gtf[tx_id].end
            chrom = gtf[tx_id].chrm

            cum = []
            total = 0
            for exon_idx in range(len(exon_starts)):
                cum.append(total)
                total += exon_ends[exon_idx] - exon_starts[exon_idx] + 1

            def _tx2g(tx_pos):
                for idx in range(len(cum) - 1, -1, -1):
                    if cum[idx] <= tx_pos:
                        return exon_starts[idx] + (tx_pos - cum[idx])
                return None

            atg_gp = _tx2g(genomic_f)
            if atg_gp is not None:
                atgstop_lines.append(
                    chrom + "\t" + str(atg_gp) + "\t" + str(atg_gp) + "\t" + orf + "\tboundaries\t" + strand + "\n"
                )
            stop_gp = _tx2g(genomic_f + orf_len_nt - 1)
            if stop_gp is not None:
                atgstop_lines.append(
                    chrom + "\t" + str(stop_gp) + "\t" + str(stop_gp) + "\t" + orf + "\tboundaries\t" + strand + "\n"
                )

            coord_psites.setdefault(orf, [])
            if strand == "+":
                psite_offsets = range(genomic_f + 2, genomic_f + orf_len_nt - 1, 3)
            else:
                psite_offsets = range(genomic_f, genomic_f + orf_len_nt, 3)
            for tx_pos in psite_offsets:
                gp = _tx2g(tx_pos)
                if gp is not None:
                    coord_psites[orf].append(gp)

    return candidates, trans_orfs, second_names, coord_psites, atgstop_lines


def make_bed(
    prot,
    trans,
    gtf,
    len_cutoff,
    max_len_cutoff,
    calculate_coordinates,
    orfs_bed_file,
    out_name,
    mult,
    genomic,
    fgenomic,
):
    print(
        "Matching ORF to transcripts, it can take a while... STEP 1 - "
        + str(datetime.now())
    )
    nomap = open(out_name + "_unmapped", "w+")
    altmap = open(out_name + "_altmapped", "w+")
    lines = []
    multiple = []
    if fgenomic == "none":
        trans_sequences = []
        seqs = {}
        for t in trans:
            if not t in gtf:
                continue
            trans_sequences.append(str(trans[t].seq.translate(cds=False)))
            trans_sequences.append(str(trans[t].seq[1:].translate(cds=False)))
            trans_sequences.append(str(trans[t].seq[2:].translate(cds=False)))
            seqs[str(trans[t].seq.translate(cds=False))] = t
            seqs[str(trans[t].seq[1:].translate(cds=False))] = t
            seqs[str(trans[t].seq[2:].translate(cds=False))] = t
    if genomic != "none":
        fgenomic = genomic
    if fgenomic != "none":
        gen_bed = {}
        for line in open(fgenomic):
            name = line.split("\t")[3] + "--" + line.split("\t")[4]
            if not name in gen_bed:
                gen_bed[name] = [
                    line.split("\t")[0],
                    int(line.split("\t")[1]) - 10,
                    int(line.split("\t")[2]) + 10,
                ]
            if int(line.split("\t")[1]) - 10 < gen_bed[name][1]:
                gen_bed[name][1] = int(line.split("\t")[1]) - 10
            if int(line.split("\t")[2]) + 10 > gen_bed[name][2]:
                gen_bed[name][2] = int(line.split("\t")[2]) + 10
    n2 = 0
    for p in prot:
        if n2 % 5000 == 0:
            print(
                str(n2)
                + " out of "
                + str(len(prot))
                + " ORFs assigned to a transcript - "
                + str(datetime.now())
            )
        n2 += 1
        disc = 0
        orf_seq = str(prot[p].seq).replace("*", "") + "*"
        if len(orf_seq) < (len_cutoff + 1):
            nomap.write(p + "\tnanoORF\t" + orf_seq + "\n")
            continue
        elif len(orf_seq) > max_len_cutoff:
            nomap.write(p + "\tlongORF\t" + orf_seq + "\n")
            continue
        elif (calculate_coordinates == "ATG") and (orf_seq[0] != "M"):
            nomap.write(p + "\tNTG\t" + orf_seq + "\n")
            continue
        elif (calculate_coordinates == "NTG") and (orf_seq[0] == "M"):
            nomap.write(p + "\tATG\t" + orf_seq + "\n")
            continue
        det = 0
        if fgenomic != "none":
            trans_sequences2 = []
            seqs2 = {}
            if not p in gen_bed:
                if genomic == "none":
                    nomap.write(p + "\tunmapped\t" + orf_seq + "\n")
                    continue
            for t in trans:
                if not t in gtf:
                    if genomic == "none":
                        nomap.write(p + "\tunmapped\t" + orf_seq + "\n")
                        continue
                elif p in gen_bed:
                    if gtf[t].chrm != gen_bed[p][0]:
                        if genomic == "none":
                            nomap.write(p + "\tunmapped\t" + orf_seq + "\n")
                            continue
                    if (
                        gtf[t].start[0] >= gen_bed[p][1]
                        and gtf[t].start[0] <= gen_bed[p][2]
                    ) or (
                        gtf[t].end[-1] >= gen_bed[p][1]
                        and gtf[t].end[-1] <= gen_bed[p][2]
                    ):
                        det = 1
                        trans_sequences2.append(str(trans[t].seq.translate(cds=False)))
                        trans_sequences2.append(
                            str(trans[t].seq[1:].translate(cds=False))
                        )
                        trans_sequences2.append(
                            str(trans[t].seq[2:].translate(cds=False))
                        )
                        seqs2[str(trans[t].seq.translate(cds=False))] = t
                        seqs2[str(trans[t].seq[1:].translate(cds=False))] = t
                        seqs2[str(trans[t].seq[2:].translate(cds=False))] = t
            if len(trans_sequences2) == 0:
                if genomic == "none":
                    nomap.write(p + "\tunmapped\t" + orf_seq + "\n")
                    continue
        if det == 0:
            res = list(filter(lambda x: orf_seq in x, trans_sequences))
        else:
            res = list(filter(lambda x: orf_seq in x, trans_sequences2))
        done_coords = []
        if len(res) > 0:
            for r in res:
                if det == 0:
                    t = seqs[r]
                else:
                    t = seqs2[r]

                l = len(str(trans[t].seq))
                f1 = str(trans[t].seq.translate(cds=False)).find(str(orf_seq))
                f2 = str(trans[t].seq[1:].translate(cds=False)).find(str(orf_seq))
                f3 = str(trans[t].seq[2:].translate(cds=False)).find(str(orf_seq))
                fi = [f1, f2, f3].index(max([f1, f2, f3]))
                f = max(f1, f2, f3) * 3 + fi
                if f < 0:
                    continue
                o = [f + 1, f + (len(orf_seq) * 3)]

                if gtf[t].strand == "-":
                    o[0], o[1] = o[1], o[0]
                    o[0] = l - o[0] + 1
                    o[1] = l - o[1] + 1

                cumu = 0
                op = 0
                orf_coords = ([], [], gtf[t].chrm, gtf[t].strand)
                for n in range(len(gtf[t].start)):
                    if op == 0:
                        if (cumu + (gtf[t].end[n] - gtf[t].start[n] + 1)) >= o[0]:
                            orf_coords[0].append(gtf[t].start[n] + o[0] - cumu - 1)
                            op = 1
                    elif op == 1:
                        orf_coords[0].append(gtf[t].start[n])

                    if op == 1:
                        if (cumu + (gtf[t].end[n] - gtf[t].start[n] + 1)) >= o[1]:
                            orf_coords[1].append(gtf[t].start[n] + o[1] - cumu - 1)
                            op = -1
                        else:
                            orf_coords[1].append(gtf[t].end[n])

                    cumu = cumu + (gtf[t].end[n] - gtf[t].start[n] + 1)

                if (
                    len(done_coords) > 0
                ):  # Some ORFs can map to multiple genomic regions
                    if (orf_coords[0][0] != done_coords[0][0]) and (
                        orf_coords[1][-1] != done_coords[1][-1]
                    ):
                        altmap.write(
                            p
                            + "\taltmap\t"
                            + orf_seq
                            + "\t"
                            + gtf[t].chrm
                            + "\t"
                            + gtf[t].strand
                            + "\t"
                            + ";".join(map(str, orf_coords[0]))
                            + "\t"
                            + ";".join(map(str, orf_coords[1]))
                            + "\t"
                            + done_coords[2]
                            + "\t"
                            + done_coords[3]
                            + "\t"
                            + ";".join(map(str, done_coords[0]))
                            + "\t"
                            + ";".join(map(str, done_coords[1]))
                            + "\n"
                        )
                        multiple.append(p)
                        break
                    for n, coord in enumerate(orf_coords[0]):
                        if (orf_coords[0][n] != done_coords[0][n]) or (
                            orf_coords[1][n] != done_coords[1][n]
                        ):
                            altmap.write(
                                p
                                + "\taltexons\t"
                                + orf_seq
                                + "\t"
                                + gtf[t].chrm
                                + "\t"
                                + gtf[t].strand
                                + "\t"
                                + ";".join(map(str, orf_coords[0]))
                                + "\t"
                                + ";".join(map(str, orf_coords[1]))
                                + "\t"
                                + done_coords[2]
                                + "\t"
                                + done_coords[3]
                                + "\t"
                                + ";".join(map(str, done_coords[0]))
                                + "\t"
                                + ";".join(map(str, done_coords[1]))
                                + "\n"
                            )
                            multiple.append(p)
                            break
                else:
                    for n, coord in enumerate(orf_coords[0]):
                        lines.append(
                            gtf[t].chrm
                            + "\t"
                            + str(orf_coords[0][n])
                            + "\t"
                            + str(orf_coords[1][n])
                            + "\t"
                            + p.split("--")[0]
                            + "\t"
                            + p.split("--")[1]
                            + "\t"
                            + gtf[t].strand
                            + "\n"
                        )
                    done_coords = orf_coords
        else:
            nomap.write(p + "\tunmapped\t" + orf_seq + "\n")
    bedout = open(orfs_bed_file, "w+")
    done = set()
    for line in lines:
        if (mult == "no") and (
            line.split("\t")[3] + "--" + line.split("\t")[4] in multiple
        ):
            if not line.split("\t")[3] + "--" + line.split("\t")[4] in done:
                nomap.write(
                    line.split("\t")[3]
                    + "--"
                    + line.split("\t")[4]
                    + "\tmultiple_coords\t"
                    + str(
                        prot[line.split("\t")[3] + "--" + line.split("\t")[4]].seq
                    ).replace("*", "")
                    + "*"
                    + "\n"
                )
            done.append(line.split("\t")[3] + "--" + line.split("\t")[4])
        else:
            bedout.write(line)
    bedout.close()
    nomap.close()
    altmap.close()


def insersect_orf_gtf(orfs_bed_file, transcriptome_gtf_file, folder):
    """Intersect a set of exonic ORFs against a transcriptome and output all partial/global overlaps"""
    print("Intersecting ORFs with transcriptome")
    overlaps = {}
    overlaps_cds = {}
    other_overlaps = {}
    total_studies = set()
    overlap_seen = {}
    seed = "".join(random.choice(string.ascii_letters) for i in range(10))
    out = open(folder + "/tmp/" + seed + "orfs_to_gtf.ov", "w+")
    subprocess.call(
        ["intersectBed", "-a", orfs_bed_file, "-b", transcriptome_gtf_file, "-wo"],
        stdout=out,
    )
    out.close()
    for line in open(folder + "/tmp/" + seed + "orfs_to_gtf.ov"):
        if line.split("\t")[5] != line.split("\t")[12]:
            continue
        if "\ttranscript\t" in line:
            name = line.split("\t")[3] + "--" + line.split("\t")[4]
            total_studies.add(line.split("\t")[4])
            gene = line.split('gene_id "')[1].split('"')[0]
            if "gene_name" in line:
                gene_name = line.split('gene_name "')[1].split('"')[0]
            else:
                gene_name = line.split('gene_id "')[1].split('"')[0]
            trans = line.split('transcript_id "')[1].split('"')[0]
            g_biotype = line.split('gene_biotype "')[1].split('"')[0]
            t_biotype = line.split('transcript_biotype "')[1].split('"')[0]
            overlaps.setdefault(name, [])
            other_overlaps.setdefault(name, ["0", "0", "0"])
            overlap_seen.setdefault(name, set())
            overlap_key = (trans, gene, t_biotype, g_biotype, gene_name)
            if overlap_key not in overlap_seen[name]:
                overlaps[name].append([trans, gene, t_biotype, g_biotype, gene_name])
                overlap_seen[name].add(overlap_key)
    for line in open(transcriptome_gtf_file):
        if "\tCDS\t" in line:
            trans = line.split('transcript_id "')[1].split('"')[0]
            prot = line.split('protein_id "')[1].split('"')[0]
            overlaps_cds[trans] = prot
    total_studies = sorted(total_studies)
    return overlaps, overlaps_cds, other_overlaps, total_studies, seed


def pseudo_or_cds_ov(
    orfs_bed_file, transcriptome_gtf_file, other_overlaps, folder, seed
):
    """Intersect with pseudogenes or CDS in any strand and any frame"""
    print("Intersecting ORFs with transcriptome")
    for line in open(folder + "/tmp/" + seed + "orfs_to_gtf.ov"):
        if (
            ("\texon\t" in line)
            and ("pseudogene" in line)
            and (line.split("\t")[5] == line.split("\t")[12])
        ):
            name = line.split("\t")[3] + "--" + line.split("\t")[4]
            other_overlaps.setdefault(name, ["0", "0", "0"])
            other_overlaps[name][0] = "1"
            if "unitary" in line:
                other_overlaps[name][0] = "2"
        if "\tCDS\t" in line:
            name = line.split("\t")[3] + "--" + line.split("\t")[4]
            if ("\t+\t" in line) and ("\t-\t" in line):
                other_overlaps.setdefault(name, ["0", "0", "0"])
                other_overlaps[name][2] = "1"
            else:
                other_overlaps.setdefault(name, ["0", "0", "0"])
                other_overlaps[name][1] = "1"
    return other_overlaps


def orf_tags(
    overlaps,
    overlaps_cds,
    orfs_fa,
    transcriptome_fa,
    proteome_fa,
    gtf,
    len_cutoff,
    max_len_cutoff,
    folder,
    seed,
):
    """Check the relative overlap of the ORFs in transcript(s)"""
    print("Checking for ORF overlaps in transcriptome")
    candidates = {}
    trans_orfs = {}
    coord_psites = {}
    second_names = {}
    atgstop_path = folder + "/tmp/" + seed + "atg_to_stop.bed"

    needed_tx_ids = set()
    needed_orf_ids = set(overlaps.keys())
    for overlap_rows in overlaps.values():
        for trans in overlap_rows:
            needed_tx_ids.add(trans[0])
    orf_seq_map = {}
    orf_desc_map = {}
    for orf_id in needed_orf_ids:
        if orf_id in orfs_fa:
            orf_seq_map[orf_id] = str(orfs_fa[orf_id].seq)
            orf_desc_map[orf_id] = str(orfs_fa[orf_id].description)
    transcript_seq_map = {
        tx_id: str(transcriptome_fa[tx_id].seq)
        for tx_id in needed_tx_ids
        if tx_id in transcriptome_fa
    }
    protein_seq_map = {}
    for tx_id, protein_id in overlaps_cds.items():
        if tx_id in transcript_seq_map and protein_id in proteome_fa:
            protein_seq_map[tx_id] = str(proteome_fa[protein_id].seq).replace("X", "")

    orf_items = list(overlaps.items())
    orf_order_map = {orf_id: idx for idx, (orf_id, _) in enumerate(orf_items)}
    fast_cpus = _resolve_fast_cpus()
    chunks = _chunk_items(orf_items, fast_cpus)
    progress = ProgressTicker("orf_tags", len(orf_items), min_interval=10.0, min_step=0.05)

    if fast_cpus > 1 and len(chunks) > 1:
        print(
            "Fast mode enabled for orf_tags with "
            + str(fast_cpus)
            + " worker(s) across "
            + str(len(orf_items))
            + " ORFs"
        )
        results = []
        with ThreadPoolExecutor(max_workers=fast_cpus) as executor:
            futures = [
                executor.submit(
                    _process_orf_chunk,
                    chunk,
                    overlaps_cds,
                    orf_seq_map,
                    orf_desc_map,
                    transcript_seq_map,
                    protein_seq_map,
                    gtf,
                    len_cutoff,
                    max_len_cutoff,
                )
                for chunk in chunks
            ]
            completed = 0
            for chunk, future in zip(chunks, futures):
                results.append(future.result())
                completed += len(chunk)
                progress.update(completed, "completed chunk")
    else:
        results = [
            _process_orf_chunk(
                orf_items,
                overlaps_cds,
                orf_seq_map,
                orf_desc_map,
                transcript_seq_map,
                protein_seq_map,
                gtf,
                len_cutoff,
                max_len_cutoff,
                progress=progress,
            )
        ]

    with open(atgstop_path, "w+") as atgstop:
        for cand_local, trans_local, names_local, psites_local, atg_lines in results:
            for orf, rows in cand_local.items():
                candidates.setdefault(orf, []).extend(rows)
            for tx_id, rows in trans_local.items():
                trans_orfs.setdefault(tx_id, []).extend(rows)
            second_names.update(names_local)
            for orf, coords in psites_local.items():
                coord_psites.setdefault(orf, []).extend(coords)
            for line in atg_lines:
                atgstop.write(line)

    trans_orfs = _canonicalize_trans_orfs(trans_orfs, orf_order_map)
    progress.update(len(orf_items), "done")

    return candidates, trans_orfs, second_names, coord_psites


def exclude_variants(
    trans_orfs, col_thr, candidates, method, coord_psites, folder, seed
):
    """Cluster ORFs"""
    print("Collapsing shorter variants")
    out = open(folder + "/tmp/" + seed + "atg_to_stop.ov", "w+")
    subprocess.call(
        [
            "intersectBed",
            "-s",
            "-a",
            folder + "/tmp/" + seed + "atg_to_stop.bed",
            "-b",
            folder + "/tmp/" + seed + "atg_to_stop.bed",
            "-wo",
        ],
        stdout=out,
    )
    out.close()
    ovs = {}
    for line in open(folder + "/tmp/" + seed + "atg_to_stop.ov"):
        if "boundaries" in line:
            left = line.split("\t")[3]
            right = line.split("\t")[9]
            ovs.setdefault(left, set()).add(right)
            ovs.setdefault(right, set()).add(left)

    exc = set()
    variants = {}
    variants_names = {}
    datasets = {}
    processed_pairs = set()
    psite_sets = {}
    total_pairs_hint = sum(len(rows) for rows in trans_orfs.values())
    progress = ProgressTicker("exclude_variants", total_pairs_hint, min_interval=10.0, min_step=0.05)
    processed_orfs = 0
    for trans in trans_orfs:
        for orf in trans_orfs[trans]:
            processed_orfs += 1
            progress.update(processed_orfs, "transcript=" + str(trans))
            orf_name = orf[0]
            for orf2_name in ovs.get(orf_name, ()):
                pair_key = tuple(sorted((orf_name, orf2_name)))
                if pair_key in processed_pairs:
                    continue
                processed_pairs.add(pair_key)
                if orf_name in exc:
                    continue
                # if (orf[0] in exc) or (orf2_name in exc):
                # 	continue
                variants.setdefault(orf_name, set())
                variants_names.setdefault(orf_name, set())
                variants_names[orf_name].add(orf_name)
                datasets.setdefault(orf_name, set())
                datasets[orf_name].add(orf_name.split("--")[1])
                orf2 = candidates[orf2_name]
                if orf_name != orf2_name:
                    if method == "longest_string":
                        longest_shared = next(iter(lcs(orf[3], orf2[0][7])), "")
                        match = len(longest_shared)
                        if match > 0:
                            if (len(orf[3]) >= len(orf2[0][7])) and (
                                match > (float(len(orf2[0][7])) * col_thr)
                            ):  # Remove shorter variant if intersect more than thr%
                                if orf2[0][7] != orf[3]:
                                    variants[orf_name].add(orf2[0][7])
                                variants_names[orf_name].add(orf2_name)
                                datasets[orf_name].add(orf2_name.split("--")[1])
                                exc.add(orf2_name)
                    elif method == "psite_overlap":
                        psite_sets.setdefault(orf_name, set(coord_psites[orf_name]))
                        psite_sets.setdefault(orf2_name, set(coord_psites[orf2_name]))
                        match = len(psite_sets[orf_name].intersection(psite_sets[orf2_name]))
                        if match > 0:
                            if (len(orf[3]) >= len(orf2[0][7])) and (
                                match > (float(len(orf2[0][7])) * col_thr)
                            ):  # Remove shorter variant if intersect more than thr%
                                if orf2[0][7] != orf[3]:
                                    variants[orf_name].add(orf2[0][7])
                                variants_names[orf_name].add(orf2_name)
                                datasets[orf_name].add(orf2_name.split("--")[1])
                                exc.add(orf2_name)

    # for var in variants:
    # 	print(var + "\t" + str(variants[var]))

    exc = sorted(exc)
    variants = {k: sorted(v) for k, v in variants.items()}
    variants_names = {k: sorted(v) for k, v in variants_names.items()}
    datasets = {k: sorted(v) for k, v in datasets.items()}
    progress.update(total_pairs_hint, "done")
    return exc, variants, variants_names, datasets


def write_output(
    orfs_fa_file,
    orfs_bed_file,
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
    len_cutoff,
    max_len_cutoff,
    col_thr,
    total_studies,
    other_overlaps,
    psites_bed_file,
    folder,
    out_name,
    method,
    seed,
    genomic,
    fgenomic,
    cds_cases,
):
    """Select main transcript and write output"""
    print("Writing output")
    # Check Riboseq ORF annotations
    done = []
    riboseq_orfs = {}
    codes = {}
    # Modified to handle missing file (single-sample mode)
    riboseq_list_file = "list_riboseqs/list_riboseq_orfs.txt"
    if os.path.exists(riboseq_list_file):
        for line in open(riboseq_list_file):
            if "\tCDS\t" in line or "\tNMD\t" in line:
                continue
            if line.startswith("#"):
                continue
            seq = line.split("\t")[4]
            if not seq in riboseq_orfs:
                riboseq_orfs[seq] = []
            riboseq_orfs[seq].append(
                [line.split("\t")[0], line.split("\t")[1], line.split("\t")[3]]
            )
            if "norep" in line:
                if (
                    not line.split("\t")[0].split("norep")[0].replace("c", "") + "norep"
                    in codes
                ):
                    codes[
                        line.split("\t")[0].split("norep")[0].replace("c", "") + "norep"
                    ] = 0
                if (
                    int(line.split("\t")[0].split("norep")[1])
                    > codes[
                        line.split("\t")[0].split("norep")[0].replace("c", "") + "norep"
                    ]
                ):
                    codes[
                        line.split("\t")[0].split("norep")[0].replace("c", "") + "norep"
                    ] = int(line.split("\t")[0].split("norep")[1])
            elif "riboseqorf" in line:
                if (
                    not line.split("\t")[0].split("riboseqorf")[0].replace("c", "")
                    + "riboseqorf"
                    in codes
                ):
                    codes[
                        line.split("\t")[0].split("riboseqorf")[0].replace("c", "")
                        + "riboseqorf"
                    ] = 0
                if (
                    int(line.split("\t")[0].split("riboseqorf")[1])
                    > codes[
                        line.split("\t")[0].split("riboseqorf")[0].replace("c", "")
                        + "riboseqorf"
                    ]
                ):
                    codes[
                        line.split("\t")[0].split("riboseqorf")[0].replace("c", "")
                        + "riboseqorf"
                    ] = int(line.split("\t")[0].split("riboseqorf")[1])

            if line.split("\t")[5].rstrip("\n") != "none":
                if ";" in line.split("\t")[5].rstrip("\n"):
                    for s in line.split("\t")[5].rstrip("\n"):
                        if not s in riboseq_orfs:
                            riboseq_orfs[s] = []
                        riboseq_orfs[s].append(
                            [
                                line.split("\t")[0] + "_var",
                                line.split("\t")[1],
                                line.split("\t")[3],
                            ]
                        )
                else:
                    if not line.split("\t")[5].rstrip("\n") in riboseq_orfs:
                        riboseq_orfs[line.split("\t")[5].rstrip("\n")] = []
                    riboseq_orfs[line.split("\t")[5].rstrip("\n")].append(
                        [
                            line.split("\t")[0] + "_var",
                            line.split("\t")[1],
                            line.split("\t")[3],
                        ]
                    )

    outs = []
    outs2 = []
    new_cases = {}
    all_biotypes = {}
    same_prot = set()
    exc = set(exc)
    out = open(out_name + ".orfs.out", "w+")
    out.write(
        "orf_id\tversion\tchrm\tstarts\tends\tstrand\ttrans\tgene\tgene_name\torf_biotype\tgene_biotype\tpep\torf_length\tn_datasets\t"
    )
    out.write("\t".join(total_studies))
    out.write(
        "\tpseudogene_ov\tCDS_ov\tCDS_as_ov\tall_trans\tall_genes\tall_gene_names\tn_variants\tseq_variants\tall_orf_names\tphaseI_id\tphaseI_biotype\n"
    )
    out3 = open(out_name + ".orfs.bed", "w+")
    out3b = open(out_name + ".orfs.gtf", "w+")
    out4 = open(out_name + ".orfs.fa", "w+")
    outf = open(out_name + ".orfs.frames.bed", "w+")
    outf2 = open(out_name + ".orfs.allframes.bed", "w+")
    outlogs = open(out_name + ".logs", "w+")
    outlogs.write(
        "#input fasta: "
        + orfs_fa_file
        + "\n#input bed: "
        + orfs_bed_file
        + "\n#annotation folder:"
        + folder
        + "\n#min_length_cutoff: "
        + str(len_cutoff)
        + "\n#max_length_cutoff: "
        + str(max_len_cutoff).replace("999999999999", "none")
        + "\n#collapse_method: "
        + str(method)
        + "\n#collapse_cutoff: "
        + str(col_thr)
        + "\n#genomic_bed : "
        + str(genomic)
        + "\n#forced_genomic_bed: "
        + str(fgenomic)
        + "\n#total_studies: "
        + ";".join(total_studies)
        + "\n"
    )
    progress = ProgressTicker("write_output", len(candidates), min_interval=10.0, min_step=0.05)
    processed_orfs = 0
    for orf in candidates:
        processed_orfs += 1
        progress.update(processed_orfs, "orf=" + str(orf))
        if orf in exc:
            continue
        all_t = set()
        all_g = set()
        all_gn = set()
        ncod = set()
        cod = set()
        for cand in candidates[orf]:
            all_t.add(cand[0])
            all_g.add(cand[1])
            all_gn.add(cand[2])
            if cand[3] == "lncRNA":
                ncod.add(cand[0])
            else:
                cod.add(cand[0])

        all_t = sorted(all_t)  # Sort list, in case of several transcripts with equal evidence, the first one is selected
        all_g = sorted(all_g)
        all_gn = sorted(all_gn)
        ncod = sorted(ncod)
        cod = sorted(cod)

        # Evidence in APPRIS and TSL
        # if (len(ncod) > 0) and (len(cod) > 0): #If the ORF overlaps coding and non-coding transcripts, prioritize coding ones
        # 	appris_list = []
        # 	for trans in all_t:
        # 		if trans in ncod:
        # 			appris_list.append("none")
        # 		elif trans in appris:
        # 			appris_list.append(appris[trans])
        # 		else:
        # 			appris_list.append("none")
        # else:
        appris_list = []
        for trans in all_t:
            if trans in appris:
                appris_list.append(appris[trans])
            else:
                appris_list.append("none")
        i = 0
        if "principal1" in appris_list:
            i2 = get_index_positions(appris_list, "principal1")
        elif "principal2" in appris_list:
            i2 = get_index_positions(appris_list, "principal2")
        elif "principal3" in appris_list:
            i2 = get_index_positions(appris_list, "principal3")
        elif "principal4" in appris_list:
            i2 = get_index_positions(appris_list, "principal4")
        elif "principal5" in appris_list:
            i2 = get_index_positions(appris_list, "principal5")
        elif "alternative1" in appris_list:
            i2 = get_index_positions(appris_list, "alternative1")
        elif "alternative2" in appris_list:
            i2 = get_index_positions(appris_list, "alternative2")
        else:
            i2 = []

        # if (len(ncod) > 0) and (len(cod) > 0): #If the ORF overlaps coding and non-coding transcripts, prioritize coding ones
        # 	appris_list = []
        # 	for trans in all_t:
        # 		if trans in ncod:
        # 			appris_list.append("0")
        # 		elif trans in supp:
        # 			appris_list.append(supp[trans].replace("tsl",""))
        # 		else:
        # 			appris_list.append("0")
        # else:
        appris_list = []
        for trans in all_t:
            if trans in supp:
                appris_list.append(supp[trans].replace("tsl", ""))
            else:
                appris_list.append("0")

        if "1" in appris_list:
            i3 = get_index_positions(appris_list, "1")
        elif "2" in appris_list:
            i3 = get_index_positions(appris_list, "2")
        elif "3" in appris_list:
            i3 = get_index_positions(appris_list, "3")
        elif "4" in appris_list:
            i3 = get_index_positions(appris_list, "4")
        elif "5" in appris_list:
            i3 = get_index_positions(appris_list, "5")
        else:
            i3 = []

        inter = []
        if len(i2) == 1:
            i = i2[0]
        elif len(i3) == 1:
            i = i3[0]
        elif len(i2) == 0 and len(i3) > 1:
            i = i3[0]  # First element
        elif len(i2) > 1 and len(i3) == 0:
            i = i2[0]  # First element
        elif len(i2) > 1 and len(i3) > 0:
            inter = [value for value in i2 if value in i3]
            if len(inter) == 1:
                i = inter[0]
            elif len(inter) > 1:
                i = inter[0]  # First element
            elif len(inter) == 0:
                i = i2[0]  # First element
        elif len(i2) == 0 and len(i3) == 0:
            i = 0  # First element

        t = all_t[i]

        for trans in candidates[orf]:  # Prioritize CDS overlaps
            if "CDS" in trans[3]:
                t = trans[0]
                break

        for trans in candidates[orf]:
            if trans[0] == t:

                # Annotate gene as pseudogene
                if other_overlaps[orf][0] == "1":
                    trans[4] = "pseudogene"
                elif other_overlaps[orf][0] == "2":
                    trans[4] = "unitary_pseudogene"
                # Vector 0/1 studies
                dataset_set = set(datasets.get(orf, []))
                stu = []
                for study in total_studies:
                    if study in dataset_set:
                        stu.append("1")
                    else:
                        stu.append("0")

                # Compare with original Riboseq annotation
                id2 = "none"
                if trans[7] in riboseq_orfs:
                    if len(riboseq_orfs[trans[7]]) == 1:
                        id2 = riboseq_orfs[trans[7]][0][0]
                    else:
                        for elemento in riboseq_orfs[trans[7]]:
                            if elemento[2] == trans[1]:
                                id2 = elemento[0]

                # Write coordinate output
                all_coords = [[], [], gtf[t].chrm, gtf[t].strand, [], []]
                if gtf[t].strand == "-":
                    rev = len(str(transcriptome_fa[trans[0]].seq)) - (
                        trans[6] + (len(trans[7]) * 3)
                    )
                    trans[6] = rev

                genomic_blocks = project_tx_range_to_genome(
                    gtf[t].start, gtf[t].end, trans[6], len(trans[7]) * 3
                )
                n3 = 1
                for start, end in genomic_blocks:
                    all_coords[4].append(
                        gtf[t].chrm
                        + "\t"
                        + str(start)
                        + "\t"
                        + str(end)
                        + "\tP1_ID\t"
                        + trans[1]
                        + "\t"
                        + gtf[t].strand
                        + "\n"
                    )
                    # all_coords[5].append(gtf[t].chrm + "\tphaseI\tCDS\t" + str(start) + "\t" + str(end) + "\t.\t" + gtf[t].strand + "\t.\tgene_id \"" + trans[1] + "\"; gene_name \"" + trans[2] + "--" + trans[3] + "\"; transcript_id \"P1_ID\";\n")
                    all_coords[5].append(
                        gtf[t].chrm
                        + "\tphaseI\tCDS\t"
                        + str(start)
                        + "\t"
                        + str(end)
                        + "\t.\t"
                        + gtf[t].strand
                        + '\t.\tgene_id "'
                        + trans[1]
                        + '"; gene_name "'
                        + trans[2]
                        + '"; gene_biotype "'
                        + trans[4]
                        + '"; transcript_id "'
                        + t
                        + '"; orf_id "P1_ID"; orf_biotype "'
                        + trans[3]
                        + '"; phaseI_id "'
                        + id2
                        + '"; exon_number "'
                        + str(n3)
                        + '";\n'
                    )
                    all_coords[0].append(str(start))
                    all_coords[1].append(str(end))
                    n3 += 1

                # Convert variants to second name
                new_variant_names = []
                for variant in variants_names.get(orf, []):
                    new_variant_names.append(second_names[variant])

                # In-frame overlaps with annotated CDSs
                id = (
                    "P1_"
                    + all_coords[2]
                    + ":"
                    + all_coords[0][0]
                    + "_"
                    + all_coords[1][-1]
                    + ":"
                    + all_coords[3]
                    + ":"
                    + str(len(all_coords[1]))
                    + ":"
                    + str(len(trans[7]) * 3)
                )
                phase = -1
                for line in all_coords[4]:
                    out3.write(line.replace("P1_ID", id))
                    # Frames
                    if phase == -1:
                        if line.split("\t")[5].rstrip("\n") == "+":
                            phase = 0
                        elif line.split("\t")[5].rstrip("\n") == "-":
                            phase = 2
                    for coord in range(
                        int(line.split("\t")[1]), int(line.split("\t")[2]) + 1
                    ):
                        if line.split("\t")[5].rstrip("\n") == "+":
                            outf2.write(
                                line.split("\t")[0]
                                + "\t"
                                + str(coord)
                                + "\t"
                                + str(coord)
                                + "\t"
                                + id
                                + "\tp"
                                + str(phase % 3)
                                + "\t"
                                + line.split("\t")[5].rstrip("\n")
                                + "\n"
                            )
                        else:
                            outf2.write(
                                line.split("\t")[0]
                                + "\t"
                                + str(coord)
                                + "\t"
                                + str(coord)
                                + "\t"
                                + id
                                + "\tp"
                                + str((phase + 1) % 3)
                                + "\t"
                                + line.split("\t")[5].rstrip("\n")
                                + "\n"
                            )
                        if phase % 3 == 2:
                            outf.write(
                                line.split("\t")[0]
                                + "\t"
                                + str(coord)
                                + "\t"
                                + str(coord)
                                + "\t"
                                + id
                                + "\t"
                                + str(phase)
                                + "\t"
                                + line.split("\t")[5].rstrip("\n")
                                + "\n"
                            )
                        phase += 1

                for line in all_coords[5]:
                    outs2.append(line.replace("P1_ID", id))

                outs.append(
                    id
                    + "\tphaseI\t"
                    + all_coords[2]
                    + "\t"
                    + ";".join(all_coords[0])
                    + "\t"
                    + ";".join(all_coords[1])
                    + "\t"
                    + all_coords[3]
                    + "\t"
                    + trans[0]
                    + "\t"
                    + trans[1]
                    + "\t"
                    + trans[2]
                    + "\t"
                    + trans[3]
                    + "\t"
                    + trans[4]
                    + "\t"
                    + trans[7]
                    + "\t"
                    + str(len(trans[7]) * 3)
                    + "\t"
                    + str(len(dataset_set))
                    + "\t"
                    + "\t".join(stu)
                    + "\t"
                    + "\t".join(other_overlaps[orf])
                    + "\t"
                    + ";".join(all_t)
                    + "\t"
                    + ";".join(all_g)
                    + "\t"
                    + ";".join(all_gn)
                    + "\t"
                    + str(len(variants[orf]))
                    + "\t"
                    + ";".join(variants[orf])
                    + "\t"
                    + ";".join(new_variant_names)
                )
                variants[trans[7]] = variants[orf]

                # Write FASTA
                out4.write(">" + id + "\n" + trans[7] + "\n")
                if (len(trans[7]) * 3 - 3) == int(trans[8]):
                    same_prot.add(trans[7])

    outf.close()
    outf2.close()
    outov = open(folder + "/tmp/" + seed + "frame_overlaps.ov", "w+")
    subprocess.call(
        [
            "intersectBed",
            "-s",
            "-a",
            out_name + ".orfs.frames.bed",
            "-b",
            psites_bed_file,
            "-wo",
        ],
        stdout=outov,
    )
    outov.close()
    cdsinf = []
    annot = {}
    for line in open(folder + "/tmp/" + seed + "frame_overlaps.ov"):
        if line.split("\t")[1] == line.split("\t")[7]:
            cdsinf.append(line.split("\t")[3])
    cdsinf = list(set(cdsinf))
    for line in outs:
        id = line.split("\t")[0]
        line2 = line

        if id in cdsinf:
            if line.split("\t")[9] != "CDS":
                line2 = line2.replace("\t" + line.split("\t")[9] + "\t", "\tCDS\t")

        # Compare with original Riboseq annotation
        if line2.split("\t")[9] in ("CDS", "NMD"):
            bio2 = "annotated"
            id2 = id
            annot[line.split("\t")[11]] = line2.split("\t")[9]
        elif line2.split("\t")[10] == "pseudogene":
            bio2 = "annotated"
            id2 = id
            annot[line.split("\t")[11]] = line2.split("\t")[10]
        elif line.split("\t")[11] in riboseq_orfs:
            if len(riboseq_orfs[line.split("\t")[11]]) == 1:
                done.add(
                    riboseq_orfs[line.split("\t")[11]][0][0].replace("_var", "")
                )
                id2 = riboseq_orfs[line.split("\t")[11]][0][0]
                bio2 = riboseq_orfs[line.split("\t")[11]][0][1]
            else:
                for elemento in riboseq_orfs[
                    line.split("\t")[11]
                ]:  # In the event of two sequences being the same, the gene_id will be inspected, be aware in case the gene_id changes across versions
                    if elemento[2] == line.split("\t")[7]:
                        done.add(elemento[0].replace("_var", ""))
                        id2 = elemento[0]
                        bio2 = elemento[1]
        else:
            if int(line.split("\t")[13]) == 1:
                key = id.split("_")[1].split(":")[0] + "norep"
                codes.setdefault(key, 0)
                id2 = "c" + id.split("_")[1].split(":")[0] + "norep" + str(codes[key] + 1)
                codes[key] += 1
            else:
                key = id.split("_")[1].split(":")[0] + "riboseqorf"
                codes.setdefault(key, 0)
                id2 = "c" + id.split("_")[1].split(":")[0] + "riboseqorf" + str(codes[key] + 1)
                codes[key] += 1
            bio2 = "novel"
            new_cases[id2] = line.split("\t")[11]

        out.write(line2 + "\t" + id2 + "\t" + bio2 + "\n")
        all_biotypes[line.split("\t")[11]] = [
            line2.split("\t")[9],
            line2.split("\t")[7],
        ]

    for line in outs2:
        id = line.split('orf_id "')[1].split('"')[0]
        bio = line.split('orf_biotype "')[1].split('"')[0]
        biog = line.split('gene_biotype "')[1].split('"')[0]
        line2 = line
        if id in cdsinf:
            if bio != "CDS":
                line2 = line2.replace("\t" + bio + "\t", "\tCDS\t")
                bio = "CDS"
        if ((bio != "CDS") and (biog != "pseudogene")) or (
            cds_cases == "yes"
        ):  # Not include CDS and non-unitary pseudogenes in GTF
            out3b.write(line2)
    others = "none"  # default; overwritten per-orf inside riboseq_orfs loop
    for orf in riboseq_orfs:
        for orf2 in riboseq_orfs[orf]:
            if "_var" in orf2[0]:
                continue
            others = "none"
            if orf in variants:
                    others = ";".join(variants[orf])

            if not orf2[0] in done:
                if orf in annot:
                    if orf in same_prot:
                        outlogs.write(
                            orf2[0]
                            + "\t"
                            + annot[orf]
                            + "\tannotated_complete\tunknown\t"
                            + orf
                            + "\t"
                            + others
                            + "\n"
                        )
                    else:
                        outlogs.write(
                            orf2[0]
                            + "\t"
                            + annot[orf]
                            + "\tannotated_alt\tunknown\t"
                            + orf
                            + "\t"
                            + others
                            + "\n"
                        )
                else:
                    outlogs.write(
                        orf2[0]
                        + "\t"
                        + orf2[1]
                        + "\tretired\tunknown\t"
                        + orf
                        + "\t"
                        + others
                        + "\n"
                    )
            elif orf in all_biotypes:
                outlogs.write(
                    orf2[0]
                    + "\t"
                    + all_biotypes[orf][0]
                    + "\tmapped\t"
                    + all_biotypes[orf][1]
                    + "\t"
                    + orf
                    + "\t"
                    + others
                    + "\n"
                )
            else:
                outlogs.write(
                    orf2[0]
                    + "\tvariant\tmapped\tunknown\t"
                    + orf
                    + "\t"
                    + others
                    + "\n"
                )
    for orf in new_cases:
        outlogs.write(
            orf
            + "\t"
            + all_biotypes[new_cases[orf]][0]
            + "\tnew\t"
            + all_biotypes[new_cases[orf]][1]
            + "\t"
            + new_cases[orf]
            + "\t"
            + others
            + "\n"
        )
    progress.update(len(candidates), "done")
    out.close()
    out3.close()
    out3b.close()
    out4.close()
    outlogs.close()

    # for f in os.listdir(folder + '/tmp/'):
    # 	if seed in f:
    # 		os.remove(os.path.join(folder + '/tmp/', f))
