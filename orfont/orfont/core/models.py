"""Core data models for ORF unification and classification.

Extracted and modularised from riboseq pipeline's unify_orf_predictions.py.
"""

import sys
import os
import bisect
from collections import defaultdict
from typing import List, Dict, Tuple, Set, Optional

from orfont.core.utils import chrom_aliases


def calculate_frame(chrom, strand, start, blocks):
    """Calculate reading frame from genomic coordinates."""
    if strand == '+':
        return (start - 1) % 3
    else:
        last_block_end = blocks[-1][1]
        return (last_block_end) % 3


# ---------------------------------------------------------------------------
# UnionFind
# ---------------------------------------------------------------------------

class UnionFind:
    """Union-Find with path compression and union by rank."""

    def __init__(self, n):
        self.parent = list(range(n))
        self.rank = [0] * n

    def find(self, x):
        if self.parent[x] != x:
            self.parent[x] = self.find(self.parent[x])
        return self.parent[x]

    def union(self, x, y):
        px, py = self.find(x), self.find(y)
        if px == py:
            return
        if self.rank[px] < self.rank[py]:
            px, py = py, px
        self.parent[py] = px
        if self.rank[px] == self.rank[py]:
            self.rank[px] += 1


# ---------------------------------------------------------------------------
# GTFIndex
# ---------------------------------------------------------------------------

class GTFIndex:
    """Load GTF into memory and provide spatial queries for ORF unification.

    Indexes transcripts, exons, CDS, and genes for:
    - transcript-to-genomic coordinate conversion
    - CDS overlap detection (in-frame check)
    - gene overlap queries
    """

    def __init__(self, gtf_file):
        self.transcripts = {}
        self.gene_map = {}
        self.gene_names = {}
        self.cds_by_chrom = {}
        self.gene_by_chrom = {}
        self.chrom_names = set()
        self._load_gtf(gtf_file)

    def _load_gtf(self, gtf_file):
        print(f"Loading GTF: {gtf_file}...", file=sys.stderr)
        cds_temp = defaultdict(list)
        gene_temp = defaultdict(list)

        with open(gtf_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) < 9:
                    continue

                feature = parts[2]
                attrs = self._parse_attributes(parts[8])
                chrom = parts[0]
                self.chrom_names.add(chrom)
                strand = parts[6]
                start, end = int(parts[3]), int(parts[4])
                gid = attrs.get('gene_id', 'NA')

                if feature == 'gene':
                    gene_temp[chrom].append((start, end, strand, gid))
                    continue
                if 'transcript_id' not in attrs:
                    continue

                tid = attrs['transcript_id']
                gname = attrs.get('gene_name', gid)
                self.gene_map[tid] = gid
                self.gene_names[gid] = gname

                if tid not in self.transcripts:
                    self.transcripts[tid] = {
                        'chrom': chrom, 'strand': strand,
                        'exons': [], 'cds': []}

                if feature == 'exon':
                    self.transcripts[tid]['exons'].append((start, end))
                elif feature == 'CDS':
                    self.transcripts[tid]['cds'].append((start, end))
                    cds_temp[chrom].append((start, end, strand, gid))

        for tid in self.transcripts:
            self.transcripts[tid]['exons'].sort()
            self.transcripts[tid]['cds'].sort()
        print(f"Loaded {len(self.transcripts)} transcripts.", file=sys.stderr)

        for chrom, intervals in cds_temp.items():
            self.cds_by_chrom[chrom] = sorted(intervals, key=lambda x: x[0])
        for chrom, intervals in gene_temp.items():
            self.gene_by_chrom[chrom] = sorted(intervals, key=lambda x: x[0])

    @staticmethod
    def _parse_attributes(attr_str):
        attrs = {}
        for p in attr_str.split(';'):
            p = p.strip()
            if not p:
                continue
            parts = p.split(' ')
            if len(parts) >= 2:
                attrs[parts[0]] = ' '.join(parts[1:]).strip('"')
        return attrs

    def resolve_chrom(self, chrom):
        for alias in chrom_aliases(chrom):
            if alias in self.chrom_names:
                return alias
        return chrom

    def get_genomic_blocks(self, tid, start_rel, end_rel, feature_type='exon'):
        """Convert transcript-relative (0-based) to genomic blocks (1-based)."""
        if tid not in self.transcripts:
            return None, None, None
        tx = self.transcripts[tid]
        chrom = tx['chrom']
        strand = tx['strand']
        blocks = tx['exons'] if feature_type == 'exon' else tx['cds']
        if not blocks and feature_type == 'cds':
            blocks = tx['exons']

        lengths = [e - s + 1 for s, e in blocks]
        total_len = sum(lengths)
        if start_rel < 0 or end_rel > total_len:
            return None, None, None

        genomic_blocks = []
        ordered_blocks = blocks if strand == '+' else blocks[::-1]
        remaining_start = start_rel
        remaining_len = end_rel - start_rel

        for s, e in ordered_blocks:
            blk_len = e - s + 1
            if remaining_start < blk_len:
                if strand == '+':
                    seg_s = s + remaining_start
                else:
                    seg_s = e - remaining_start

                fit_len = min(remaining_len, blk_len - remaining_start)
                if strand == '+':
                    seg_e = seg_s + fit_len - 1
                    genomic_blocks.append((seg_s, seg_e))
                else:
                    seg_start_genomic = seg_s - fit_len + 1
                    genomic_blocks.append((seg_start_genomic, seg_s))

                remaining_start = 0
                remaining_len -= fit_len
                if remaining_len <= 0:
                    break
            else:
                remaining_start -= blk_len

        genomic_blocks.sort()
        return chrom, strand, genomic_blocks

    def find_cds_overlap_inframe(self, chrom, strand, orf_start, orf_end, orf_frame):
        """Check if ORF overlaps in-frame with any annotated CDS."""
        chrom = self.resolve_chrom(chrom)
        intervals = self.cds_by_chrom.get(chrom, [])
        if not intervals:
            return False
        starts = [iv[0] for iv in intervals]
        idx = bisect.bisect_right(starts, orf_end)
        for i in range(min(idx, len(intervals)) - 1, -1, -1):
            iv_start, iv_end, iv_strand, _ = intervals[i]
            if iv_start > orf_end:
                continue
            if iv_end < orf_start:
                break
            if iv_strand != strand:
                continue
            cds_frame = iv_start % 3 if strand == '+' else iv_end % 3
            if cds_frame == orf_frame:
                return True
        return False

    def find_overlapping_genes(self, chrom, strand, orf_start, orf_end):
        """Return gene_ids whose genomic interval overlaps the ORF."""
        chrom = self.resolve_chrom(chrom)
        intervals = self.gene_by_chrom.get(chrom, [])
        if not intervals:
            return []
        starts = [iv[0] for iv in intervals]
        idx = bisect.bisect_right(starts, orf_end)
        result = []
        seen = set()
        for i in range(min(idx, len(intervals)) - 1, -1, -1):
            iv_start, iv_end, iv_strand, iv_gid = intervals[i]
            if iv_start > orf_end:
                continue
            if iv_end < orf_start:
                break
            if iv_strand != strand:
                continue
            if iv_gid not in seen:
                seen.add(iv_gid)
                result.append(iv_gid)
        return result


# ---------------------------------------------------------------------------
# ORFCandidate
# ---------------------------------------------------------------------------

class ORFCandidate:
    """Represents a single ORF prediction from one tool/sample."""

    def __init__(self, chrom, strand, blocks, tid, gid, tool, sample,
                 score=None, pvalue=None, sequence=None):
        self.chrom = chrom
        self.strand = strand
        self.blocks = tuple(sorted(blocks))  # 1-based (start, end)
        self.tid = tid
        self.gid = gid
        self.sources = {(tool, sample)}
        self.score = score
        self.sequence = sequence
        self.aa_sequence = ""
        self.tool_scores = {tool: score} if score is not None else {}
        self.tool_pvalues = {tool: pvalue} if pvalue is not None else {}

        self.total_psites = 0
        self.unique_psites = 0
        self.total_reads = 0
        self.unique_reads = 0
        self.subset_orfs = []
        self.is_representative = True

        self.start = self.blocks[0][0]
        self.end = self.blocks[-1][1]
        self.length_nt = sum(e - s + 1 for s, e in self.blocks)
        self.length_aa = self.length_nt // 3
        self.frame = calculate_frame(chrom, strand, self.start, self.blocks)

        self.start_codon = ""
        self.is_cds_overlap = False
        self.overlapping_gene_ids = []

        self.id_key = (self.chrom, self.strand, self.blocks)

    def merge(self, other):
        """Merge another candidate's metadata into this one."""
        self.sources.update(other.sources)
        for tool, score in other.tool_scores.items():
            if tool not in self.tool_scores or (score is not None and self.tool_scores.get(tool) is None):
                self.tool_scores[tool] = score
        for tool, pval in other.tool_pvalues.items():
            if tool not in self.tool_pvalues or (pval is not None and self.tool_pvalues.get(tool) is None):
                self.tool_pvalues[tool] = pval
        self.total_psites += other.total_psites
        self.unique_psites += other.unique_psites
        self.total_reads += other.total_reads
        self.unique_reads += other.unique_reads
        if not self.sequence and other.sequence:
            self.sequence = other.sequence
        if other.is_cds_overlap:
            self.is_cds_overlap = True
        for gid in other.overlapping_gene_ids:
            if gid not in self.overlapping_gene_ids:
                self.overlapping_gene_ids.append(gid)
        if not self.start_codon and other.start_codon:
            self.start_codon = other.start_codon
        self.subset_orfs.extend(other.subset_orfs)

    def add_subset_orf(self, other):
        """Record another ORF as a subset of this representative."""
        blocks_str = ",".join(f"{s}-{e}" for s, e in other.blocks)
        tools = ",".join(sorted({t for t, s in other.sources}))
        samples = ",".join(sorted({s for t, s in other.sources}))
        self.subset_orfs.append({
            'blocks': blocks_str, 'start': other.start, 'end': other.end,
            'length_aa': other.length_aa, 'tools': tools, 'samples': samples,
            'tool_scores': dict(other.tool_scores),
            'tool_pvalues': dict(other.tool_pvalues),
        })
        self.sources.update(other.sources)
        if other.is_cds_overlap:
            self.is_cds_overlap = True
        for gid in other.overlapping_gene_ids:
            if gid not in self.overlapping_gene_ids:
                self.overlapping_gene_ids.append(gid)
        for tool, score in other.tool_scores.items():
            if tool not in self.tool_scores or (score is not None and self.tool_scores.get(tool) is None):
                self.tool_scores[tool] = score
        for tool, pval in other.tool_pvalues.items():
            if tool not in self.tool_pvalues or (pval is not None and self.tool_pvalues.get(tool) is None):
                self.tool_pvalues[tool] = pval

    @property
    def pN(self):
        return self.total_psites / self.length_nt if self.length_nt > 0 else 0

    @property
    def unique_pN(self):
        return self.unique_psites / self.length_nt if self.length_nt > 0 else 0


# ---------------------------------------------------------------------------
# BedgraphIndex
# ---------------------------------------------------------------------------

class BedgraphIndex:
    """Pre-loaded binned bedgraph for fast P-site queries."""

    BIN_SIZE = 10000

    def __init__(self, bedgraph_file):
        self.data = defaultdict(lambda: defaultdict(list))
        self.loaded = False
        if not bedgraph_file or not os.path.exists(bedgraph_file):
            return
        try:
            with open(bedgraph_file, 'r') as f:
                for line in f:
                    if line.startswith('track') or line.startswith('#'):
                        continue
                    parts = line.strip().split('\t')
                    if len(parts) < 4:
                        continue
                    chrom = parts[0]
                    start = int(parts[1])
                    end = int(parts[2])
                    value = float(parts[3])
                    start_bin = start // self.BIN_SIZE
                    end_bin = (end - 1) // self.BIN_SIZE
                    for bi in range(start_bin, end_bin + 1):
                        self.data[chrom][bi].append((start, end, value))
            self.loaded = True
        except Exception as e:
            print(f"Warning: Error loading bedgraph: {e}", file=sys.stderr)

    def count_in_region(self, chrom, start_1based, end_1based):
        """Count signal in genomic region (1-based coords)."""
        if not self.loaded:
            return 0
        start_0 = start_1based - 1
        end_0 = end_1based
        total = 0
        start_bin = start_0 // self.BIN_SIZE
        end_bin = (end_0 - 1) // self.BIN_SIZE
        seen = set()
        for bi in range(start_bin, end_bin + 1):
            for entry in self.data[chrom].get(bi, []):
                eid = id(entry)
                if eid in seen:
                    continue
                seen.add(eid)
                bg_start, bg_end, bg_value = entry
                overlap_start = max(bg_start, start_0)
                overlap_end = min(bg_end, end_0)
                if overlap_start < overlap_end:
                    total += bg_value * (overlap_end - overlap_start)
        return int(total)
