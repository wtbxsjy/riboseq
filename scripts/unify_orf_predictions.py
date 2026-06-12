#!/usr/bin/env python3
import argparse
import sys
import os
import csv
import re
import gc
import gzip
from typing import List, Dict, Tuple, Set, Optional
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
import multiprocessing
import threading


def _open(path, mode='r'):
    if path.endswith('.gz'):
        return gzip.open(path, mode + 't')
    return open(path, mode)


def _detect_file_format(file_path):
    """Detect file format by reading content, not filename extension.

    Returns one of: 'gtf', 'tsv', 'unknown'.
    """
    try:
        with _open(file_path, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('##'):
                    continue
                if line.startswith('#'):
                    continue
                parts = line.split('\t')
                if len(parts) == 9 and parts[2] in (
                    'CDS', 'exon', 'ORF', 'start_codon', 'stop_codon',
                    'gene', 'transcript', 'five_prime_utr', 'three_prime_utr',
                ):
                    return 'gtf'
                if parts[0] in ('Gene', 'Gid', 'ORF_ID'):
                    return 'tsv'
                return 'unknown'
    except Exception:
        pass
    return 'unknown'

# ---------------------------------------------------------------------------
# Module-level globals shared with parallel worker processes via fork (Linux).
# Set in main() before the Pool is created; inherited by child processes via
# copy-on-write fork — no serialisation overhead for large objects.
# ---------------------------------------------------------------------------
_shared_gtf_index  = None   # GTFIndex instance
_shared_fasta_path = None   # str path to genome FASTA (workers open their own handle)

# Maps the canonical tool display-name used in tasks to the lowercase key used
# in tool_stats / sample_stats dictionaries.
_TOOL_STATS_KEY = {
    'Ribo-TISH': 'ribotish',
    'Ribotricer': 'ribotricer',
    'RiboCode':   'ribocode',
    'ORFquant':   'orfquant',
    'PRICE':      'price',
}

# Try to import biopython and pyfaidx
try:
    from Bio import SeqIO
    from Bio.Seq import Seq
    _BIOPYTHON_IMPORT_ERROR = None
except ImportError as exc:
    SeqIO = None
    Seq = None
    _BIOPYTHON_IMPORT_ERROR = exc

try:
    from pyfaidx import Fasta
    _PYFAIDX_IMPORT_ERROR = None
except ImportError as exc:
    Fasta = None
    _PYFAIDX_IMPORT_ERROR = exc


def _chrom_aliases(chrom: Optional[str]) -> List[str]:
    """Generate common chromosome-name aliases used across mouse references."""
    if not chrom:
        return []

    aliases = [chrom]
    if chrom.startswith('chr'):
        base = chrom[3:]
        aliases.append(base)
        if base == 'M':
            aliases.append('MT')
        elif base == 'MT':
            aliases.append('M')
    else:
        aliases.append(f"chr{chrom}")
        if chrom == 'MT':
            aliases.extend(['chrM', 'M'])
        elif chrom == 'M':
            aliases.extend(['chrM', 'MT'])

    seen = set()
    result = []
    for alias in aliases:
        if alias and alias not in seen:
            seen.add(alias)
            result.append(alias)
    return result


def infer_sample_id_from_prediction_path(file_path: str, tool_suffix: str) -> str:
    """
    Recover the sample ID from an ORF-prediction filename without truncating
    dots that are part of the original sample name.
    """
    name = os.path.basename(file_path)
    if name.endswith(tool_suffix):
        return name[:-len(tool_suffix)]
    return os.path.splitext(name)[0]

# GTF Parsing Helper
class GTFIndex:
    def __init__(self, gtf_file):
        self.transcripts = {} # tid -> {chrom, strand, exons: [(s,e), ...], cds: [(s,e), ...]}
        self.gene_map = {} # tid -> gid
        self.gene_names = {} # gid -> gene_name
        self.cds_by_chrom = {}   # chrom -> sorted list of (start, end, strand, gid)
        self.gene_by_chrom = {}  # chrom -> sorted list of (start, end, strand, gid)
        self.chrom_names = set()
        self._load_gtf(gtf_file)

    def _load_gtf(self, gtf_file):
        print(f"Loading GTF: {gtf_file}...", file=sys.stderr)
        cds_temp = defaultdict(list)
        gene_temp = defaultdict(list)
        with _open(gtf_file, 'r') as f:
            for line in f:
                if line.startswith('#'): continue
                parts = line.strip().split('\t')
                if len(parts) < 9: continue

                feature = parts[2]
                attributes = self._parse_attributes(parts[8])
                chrom = parts[0]
                self.chrom_names.add(chrom)
                strand = parts[6]
                start, end = int(parts[3]), int(parts[4])
                gid = attributes.get('gene_id', 'NA')

                # gene features have no transcript_id – collect for spatial index
                if feature == 'gene':
                    gene_temp[chrom].append((start, end, strand, gid))
                    continue

                if 'transcript_id' not in attributes: continue
                tid = attributes['transcript_id']
                gname = attributes.get('gene_name', gid)

                self.gene_map[tid] = gid
                self.gene_names[gid] = gname

                if tid not in self.transcripts:
                    self.transcripts[tid] = {
                        'chrom': chrom,
                        'strand': strand,
                        'exons': [],
                        'cds': []
                    }

                if feature == 'exon':
                    self.transcripts[tid]['exons'].append((start, end))
                elif feature == 'CDS':
                    self.transcripts[tid]['cds'].append((start, end))
                    cds_temp[chrom].append((start, end, strand, gid))

        # Sort exons and CDS
        for tid in self.transcripts:
            self.transcripts[tid]['exons'].sort()
            self.transcripts[tid]['cds'].sort()
        print(f"Loaded {len(self.transcripts)} transcripts.", file=sys.stderr)

        # Build spatial indexes (sorted by start for bisect queries)
        import bisect
        for chrom, intervals in cds_temp.items():
            self.cds_by_chrom[chrom] = sorted(intervals, key=lambda x: x[0])
        for chrom, intervals in gene_temp.items():
            self.gene_by_chrom[chrom] = sorted(intervals, key=lambda x: x[0])
        print(f"CDS spatial index: {sum(len(v) for v in self.cds_by_chrom.values())} intervals on {len(self.cds_by_chrom)} chroms.", file=sys.stderr)

    def _parse_attributes(self, attr_str):
        attrs = {}
        for p in attr_str.split(';'):
            p = p.strip()
            if not p: continue
            parts = p.split(' ')
            if len(parts) >= 2:
                key = parts[0]
                val = ' '.join(parts[1:]).strip('"')
                attrs[key] = val
        return attrs

    def resolve_chrom(self, chrom):
        """Resolve common chr/non-chr/mitochondrial aliases against loaded GTF contigs."""
        for alias in _chrom_aliases(chrom):
            if alias in self.chrom_names:
                return alias
        return chrom

    def get_genomic_blocks(self, tid, start_rel, end_rel, feature_type='exon'):
        """
        Convert relative coordinates (0-based) on a transcript (or CDS) to genomic blocks (1-based).
        start_rel: 0-based start index on the spliced sequence
        end_rel: 0-based end index (exclusive) on the spliced sequence
        feature_type: 'exon' (relative to full transcript) or 'cds' (relative to CDS start)
        """
        if tid not in self.transcripts:
            return None, None, None
        
        tx = self.transcripts[tid]
        chrom = tx['chrom']
        strand = tx['strand']
        
        blocks = tx['exons'] if feature_type == 'exon' else tx['cds']
        if not blocks and feature_type == 'cds':
             blocks = tx['exons']

        # Calculate lengths
        lengths = [e - s + 1 for s, e in blocks]
        total_len = sum(lengths)
        
        if start_rel < 0 or end_rel > total_len:
            # Out of bounds
            return None, None, None

        genomic_blocks = []
        
        # Process blocks in transcriptional order
        ordered_blocks = blocks if strand == '+' else blocks[::-1]
        
        remaining_start = start_rel
        remaining_len = end_rel - start_rel
        
        for s, e in ordered_blocks:
            blk_len = e - s + 1
            
            # Check if our segment starts in this block or before
            if remaining_start < blk_len:
                # Our segment starts in this block
                if strand == '+':
                    seg_s = s + remaining_start
                else:
                    seg_s = e - remaining_start 
                
                # How much of the segment fits in this block?
                fit_len = min(remaining_len, blk_len - remaining_start)
                
                if strand == '+':
                    seg_e = seg_s + fit_len - 1
                    genomic_blocks.append((seg_s, seg_e))
                else:
                    seg_e = seg_s # seg_s is actually the higher coordinate
                    seg_start_genomic = seg_e - fit_len + 1
                    genomic_blocks.append((seg_start_genomic, seg_e))
                
                remaining_start = 0 
                remaining_len -= fit_len
                
                if remaining_len <= 0:
                    break
            else:
                remaining_start -= blk_len
        
        # Sort genomic blocks by coordinate
        genomic_blocks.sort()
        
        return chrom, strand, genomic_blocks

    def find_cds_overlap_inframe(self, chrom, strand, orf_start, orf_end, orf_frame):
        """Check if the ORF overlaps in-frame with any annotated CDS interval."""
        import bisect
        chrom = self.resolve_chrom(chrom)
        intervals = self.cds_by_chrom.get(chrom, [])
        if not intervals:
            return False
        starts = [iv[0] for iv in intervals]
        idx = bisect.bisect_right(starts, orf_end)
        for i in range(min(idx, len(intervals)) - 1, -1, -1):
            iv_start, iv_end, iv_strand, iv_gid = intervals[i]
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
        """Return list of gene_ids whose genomic interval overlaps the ORF (same strand)."""
        import bisect
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


def calculate_frame(chrom, strand, start, blocks):
    """
    Calculate the reading frame (0, 1, 2) of an ORF.
    Frame is determined by the start codon position modulo 3.
    For negative strand, we use the end position of the last block.
    """
    if strand == '+':
        return start % 3
    else:
        # For negative strand, frame is based on the end of the ORF
        end = blocks[-1][1] if blocks else start
        return end % 3


def get_frame_aware_key(cand, tolerance=3):
    """
    Generate a key for frame-aware grouping.
    ORFs with same chrom, strand, frame, and start within tolerance are grouped.
    
    Returns: (chrom, strand, frame, rounded_start, block_count)
    """
    frame = calculate_frame(cand.chrom, cand.strand, cand.start, cand.blocks)
    # Round start to nearest multiple of 3 (frame-aligned)
    rounded_start = (cand.start // tolerance) * tolerance
    return (cand.chrom, cand.strand, frame, rounded_start, len(cand.blocks))


def are_frame_compatible(cand1, cand2, tolerance=3):
    """
    Check if two ORF candidates are frame-compatible for merging.
    
    Conditions:
    1. Same chromosome and strand
    2. Same reading frame
    3. Start positions differ by <= tolerance bases
    4. End positions differ by <= tolerance bases  
    5. Same number of exon blocks
    """
    if cand1.chrom != cand2.chrom or cand1.strand != cand2.strand:
        return False
    
    if len(cand1.blocks) != len(cand2.blocks):
        return False
    
    frame1 = calculate_frame(cand1.chrom, cand1.strand, cand1.start, cand1.blocks)
    frame2 = calculate_frame(cand2.chrom, cand2.strand, cand2.start, cand2.blocks)
    
    if frame1 != frame2:
        return False
    
    # Check start and end within tolerance
    if abs(cand1.start - cand2.start) > tolerance:
        return False
    if abs(cand1.end - cand2.end) > tolerance:
        return False
    
    # Check each block boundary (for multi-exon ORFs)
    for (s1, e1), (s2, e2) in zip(cand1.blocks, cand2.blocks):
        if abs(s1 - s2) > tolerance or abs(e1 - e2) > tolerance:
            return False
    
    return True


def calculate_overlap(cand1, cand2):
    """
    Calculate overlap between two ORF candidates.
    Returns: (overlap_bp, overlap_fraction)
    overlap_fraction is relative to the shorter ORF.
    """
    if cand1.chrom != cand2.chrom or cand1.strand != cand2.strand:
        return 0, 0.0
    
    # Calculate genomic overlap
    overlap_start = max(cand1.start, cand2.start)
    overlap_end = min(cand1.end, cand2.end)
    
    if overlap_start >= overlap_end:
        return 0, 0.0
    
    # For single-exon ORFs, simple overlap calculation
    if len(cand1.blocks) == 1 and len(cand2.blocks) == 1:
        overlap_bp = overlap_end - overlap_start
        min_len = min(cand1.length_nt, cand2.length_nt)
        return overlap_bp, overlap_bp / min_len if min_len > 0 else 0.0
    
    # For multi-exon ORFs, use two-pointer merge for O(b1+b2) overlap calculation
    overlap_bp = 0
    b1 = sorted(cand1.blocks)
    b2 = sorted(cand2.blocks)
    i2, j2 = 0, 0
    while i2 < len(b1) and j2 < len(b2):
        s1, e1 = b1[i2]
        s2, e2 = b2[j2]
        os = max(s1, s2)
        oe = min(e1, e2)
        if os < oe:
            overlap_bp += oe - os
        if e1 <= e2:
            i2 += 1
        else:
            j2 += 1
    
    min_len = min(cand1.length_nt, cand2.length_nt)
    return overlap_bp, overlap_bp / min_len if min_len > 0 else 0.0


class UnionFind:
    """Union-Find data structure for grouping overlapping ORFs."""
    
    def __init__(self, n):
        self.parent = list(range(n))
        self.rank = [0] * n
    
    def find(self, x):
        if self.parent[x] != x:
            self.parent[x] = self.find(self.parent[x])  # Path compression
        return self.parent[x]
    
    def union(self, x, y):
        px, py = self.find(x), self.find(y)
        if px == py:
            return
        # Union by rank
        if self.rank[px] < self.rank[py]:
            px, py = py, px
        self.parent[py] = px
        if self.rank[px] == self.rank[py]:
            self.rank[px] += 1


def group_overlapping_orfs(candidates, min_overlap_fraction=0.5):
    """
    Group overlapping ORFs using Union-Find.
    
    Args:
        candidates: List of ORFCandidate
        min_overlap_fraction: Minimum overlap fraction to group ORFs
    
    Returns:
        List of groups, where each group is a list of ORFCandidate indices
    """
    n = len(candidates)
    if n == 0:
        return []
    
    uf = UnionFind(n)
    
    # Sort by chromosome and position for efficient comparison
    sorted_indices = sorted(range(n), key=lambda i: (candidates[i].chrom, candidates[i].strand, candidates[i].start))
    
    # Only compare nearby candidates (within max ORF length window)
    for i, idx1 in enumerate(sorted_indices):
        cand1 = candidates[idx1]
        
        # Compare with subsequent candidates on same chromosome/strand
        for j in range(i + 1, n):
            idx2 = sorted_indices[j]
            cand2 = candidates[idx2]
            
            # Stop if different chromosome or strand
            if cand2.chrom != cand1.chrom or cand2.strand != cand1.strand:
                break
            
            # Stop if too far apart (no possible overlap)
            if cand2.start > cand1.end:
                break
            
            # Check overlap
            _, overlap_frac = calculate_overlap(cand1, cand2)
            if overlap_frac >= min_overlap_fraction:
                uf.union(idx1, idx2)
    
    # Collect groups
    groups = defaultdict(list)
    for i in range(n):
        groups[uf.find(i)].append(i)
    
    return list(groups.values())


def merge_frame_compatible_orfs(candidates, min_overlap_fraction=0.9):
    """
    Merge ORFs that are frame-compatible with >= min_overlap_fraction positional overlap.

    Rules (Mudge et al. 2022-inspired):
    - ONLY single-exon ORFs are eligible for frame-aware merging
    - Multi-exon (spliced) ORFs pass through unchanged
    - Single-exon ORFs are merged if:
        1. Same chromosome and strand
        2. Same reading frame
        3. Overlap fraction >= min_overlap_fraction (relative to shorter ORF)
    - The longest ORF in a merged group becomes representative

    Args:
        candidates: List of ORFCandidate
        min_overlap_fraction: Minimum overlap of shorter ORF to merge (default 0.9 = 90%)

    Returns:
        List of ORFCandidate (multi-exon untouched; single-exon merged)
    """
    if not candidates:
        return []

    # Separate single-exon from multi-exon: multi-exon pass through
    single_exon = [c for c in candidates if len(c.blocks) == 1]
    multi_exon  = [c for c in candidates if len(c.blocks) > 1]

    if not single_exon:
        return multi_exon

    # Sort single-exon by chrom, strand, frame, start for efficient sweep
    single_exon.sort(key=lambda c: (c.chrom, c.strand, c.frame, c.start))

    merged = []
    used = set()

    for i, cand1 in enumerate(single_exon):
        if i in used:
            continue

        compatible_group = [cand1]
        max_possible_end = cand1.end  # tracks the furthest end in current group

        for j in range(i + 1, len(single_exon)):
            if j in used:
                continue
            cand2 = single_exon[j]

            # Different chrom/strand/frame → can never merge
            if cand2.chrom != cand1.chrom or cand2.strand != cand1.strand:
                break
            if cand2.frame != cand1.frame:
                break  # sorted by (chrom, strand, frame, start); no more matching frames
            # Early exit: cand2 starts beyond the furthest end of any group member
            if cand2.start > max_possible_end:
                break

            # Compute overlap fraction (relative to shorter ORF)
            overlap_start = max(cand1.start, cand2.start)
            overlap_end   = min(max_possible_end, cand2.end)
            if overlap_start >= overlap_end:
                continue
            overlap_bp = overlap_end - overlap_start
            min_len = min(cand1.length_nt, cand2.length_nt)
            if min_len > 0 and overlap_bp / min_len >= min_overlap_fraction:
                compatible_group.append(cand2)
                used.add(j)
                max_possible_end = max(max_possible_end, cand2.end)

        used.add(i)

        if len(compatible_group) == 1:
            merged.append(cand1)
        else:
            # Longest ORF is representative; others are merged into it
            representative = max(compatible_group, key=lambda c: c.length_nt)
            for cand in compatible_group:
                if cand is not representative:
                    representative.merge(cand)
            merged.append(representative)

    return merged + multi_exon


def select_representative_for_group(candidates, group_indices):
    """
    Select the representative ORF for a group of overlapping ORFs.
    The longest ORF becomes representative, shorter ones become subsets.
    
    Args:
        candidates: List of all ORFCandidate
        group_indices: List of indices belonging to this group
    
    Returns:
        The representative ORFCandidate with subset information recorded
    """
    if len(group_indices) == 1:
        return candidates[group_indices[0]]
    
    # Sort by length (longest first)
    group_cands = [candidates[i] for i in group_indices]
    group_cands.sort(key=lambda c: c.length_nt, reverse=True)
    
    representative = group_cands[0]
    
    # Add all other ORFs as subsets
    for cand in group_cands[1:]:
        representative.add_subset_orf(cand)
        cand.is_representative = False
    
    return representative


def process_overlap_groups(candidates, min_overlap_fraction=0.5):
    """
    Process all candidates: group overlapping ORFs and select representatives.
    
    Args:
        candidates: List of ORFCandidate
        min_overlap_fraction: Minimum overlap fraction to group ORFs
    
    Returns:
        List of representative ORFCandidate (subset info recorded in each)
    """
    # Group overlapping ORFs
    groups = group_overlapping_orfs(candidates, min_overlap_fraction)
    
    # Select representative for each group
    representatives = []
    for group_indices in groups:
        rep = select_representative_for_group(candidates, group_indices)
        representatives.append(rep)
    
    return representatives


def _seq_identity_from_shared_terminus(aa1, aa2):
    """
    Compute sequence identity aligned from the shared terminus.
    For two ORFs sharing the same stop codon, compare from the C-terminus;
    for same start, compare from the N-terminus.
    Returns fraction of identical positions / len(shorter).
    Caller is responsible for picking the correct direction.
    """
    if not aa1 or not aa2:
        return 0.0
    shorter, longer = (aa1, aa2) if len(aa1) <= len(aa2) else (aa2, aa1)
    n = len(shorter)
    # Compare last n characters of longer with shorter (shared stop → C-terminal)
    matches = sum(a == b for a, b in zip(shorter, longer[-n:]))
    return matches / n


def cluster_by_sequence(candidates, seq_identity=0.9, min_length_ratio=0.8,
                        terminus_tolerance=3):
    """
    Stage 3: Sequence-similarity-based clustering (UniRef90-inspired).

    Two ORFs are merged into one cluster when ALL THREE conditions hold:
      1. Same gene locus (gene_id match, or same chrom+strand for gene_id='NA')
      2. Share a start OR stop position within ±terminus_tolerance bp
      3. Length ratio (shorter/longer) >= min_length_ratio
         AND sequence identity (identical_aa/len_shorter) >= seq_identity

    Within each cluster the longest ORF is the representative; others are
    recorded as subset_orfs.
    """
    if not candidates:
        return []

    # ---- bucket by gene_id (fall back to chrom+strand for NA) ----
    from collections import defaultdict
    buckets = defaultdict(list)
    for idx, cand in enumerate(candidates):
        key = cand.gid if cand.gid and cand.gid != 'NA' else f"_locus_{cand.chrom}_{cand.strand}"
        buckets[key].append(idx)

    uf = UnionFind(len(candidates))

    for key, indices in buckets.items():
        n = len(indices)
        if n == 1:
            continue
        # Replace O(n²) pairwise scan with terminus-based hash bucketing.
        # Two ORFs can only merge if they share a start OR stop within
        # ±terminus_tolerance bp (Condition 2).  Instead of checking all
        # pairs, bucket by floor(pos / bin_width) and only compare within
        # the same or adjacent bins — reduces to O(n) average per gene.
        tol = terminus_tolerance
        bin_w = 2 * tol + 1

        def _check_pair(idx1, idx2):
            c1, c2 = candidates[idx1], candidates[idx2]
            len1, len2 = c1.length_aa, c2.length_aa
            if len1 == 0 or len2 == 0:
                return
            if min(len1, len2) / max(len1, len2) < min_length_ratio:
                return
            aa1, aa2 = c1.aa_sequence, c2.aa_sequence
            if aa1 and aa2:
                if _seq_identity_from_shared_terminus(aa1, aa2) < seq_identity:
                    return
            uf.union(idx1, idx2)

        # --- shared-start pairs ---
        start_buckets = defaultdict(list)
        for idx in indices:
            c = candidates[idx]
            start_buckets[c.start // bin_w].append(idx)
        for bin_key, bucket in list(start_buckets.items()):
            # check within this bin and against the next bin (boundary pairs)
            next_bucket = start_buckets.get(bin_key + 1, [])
            combined = bucket + next_bucket
            for ii in range(len(bucket)):
                for jj in range(ii + 1, len(combined)):
                    idx1, idx2 = bucket[ii], combined[jj]
                    if idx1 == idx2:
                        continue
                    if abs(candidates[idx1].start - candidates[idx2].start) <= tol:
                        _check_pair(idx1, idx2)

        # --- shared-stop pairs ---
        end_buckets = defaultdict(list)
        for idx in indices:
            c = candidates[idx]
            end_buckets[c.end // bin_w].append(idx)
        for bin_key, bucket in list(end_buckets.items()):
            next_bucket = end_buckets.get(bin_key + 1, [])
            combined = bucket + next_bucket
            for ii in range(len(bucket)):
                for jj in range(ii + 1, len(combined)):
                    idx1, idx2 = bucket[ii], combined[jj]
                    if idx1 == idx2:
                        continue
                    if abs(candidates[idx1].end - candidates[idx2].end) <= tol:
                        _check_pair(idx1, idx2)

    # Collect groups
    from collections import defaultdict
    group_map = defaultdict(list)
    for idx in range(len(candidates)):
        group_map[uf.find(idx)].append(idx)

    representatives = []
    for root, group_indices in group_map.items():
        rep = select_representative_for_group(candidates, group_indices)
        representatives.append(rep)

    return representatives


class ORFCandidate:
    def __init__(self, chrom, strand, blocks, tid, gid, tool, sample, score=None, pvalue=None, sequence=None):
        self.chrom = chrom
        self.strand = strand
        self.blocks = tuple(sorted(blocks)) # List of (start, end) tuples, 1-based
        self.tid = tid
        self.gid = gid
        self.sources = {(tool, sample)} # Set of (tool, sample) tuples
        self.score = score
        self.sequence = sequence # Extracted nucleotide sequence
        self.aa_sequence = ""   # Amino acid sequence (set by extract_sequence)
        self.tool_scores = {tool: score} if score is not None else {}  # Dict: tool -> score
        self.tool_pvalues = {tool: pvalue} if pvalue is not None else {}  # Dict: tool -> pvalue
        
        # Statistics from bedgraph (calculated later)
        self.total_psites = 0
        self.unique_psites = 0
        self.total_reads = 0
        self.unique_reads = 0
        
        # Subset ORFs (for representatives)
        self.subset_orfs = []  # List of (blocks_str, tools, samples) for subset ORFs
        self.is_representative = True  # Whether this is the representative of its group
        
        # Calculated fields
        self.start = self.blocks[0][0]
        self.end = self.blocks[-1][1]
        self.length_nt = sum(e - s + 1 for s, e in self.blocks)
        self.length_aa = self.length_nt // 3
        
        # Frame (0, 1, 2)
        self.frame = calculate_frame(chrom, strand, self.start, self.blocks)

        self.start_codon = ""         # First 3 nt of coding sequence (set after sequence extraction)
        self.is_cds_overlap = False   # Whether ORF overlaps in-frame with annotated CDS
        self.overlapping_gene_ids = [] # Gene IDs whose genomic region overlaps this ORF

        # ID for grouping (exact match)
        self.id_key = (self.chrom, self.strand, self.blocks)

    def merge(self, other):
        """Merge another candidate into this one."""
        self.sources.update(other.sources)
        # Merge tool scores
        for tool, score in other.tool_scores.items():
            if tool not in self.tool_scores or (score is not None and self.tool_scores.get(tool) is None):
                self.tool_scores[tool] = score
        # Merge tool pvalues
        for tool, pval in other.tool_pvalues.items():
            if tool not in self.tool_pvalues or (pval is not None and self.tool_pvalues.get(tool) is None):
                self.tool_pvalues[tool] = pval
        # Merge statistics (sum across samples)
        self.total_psites += other.total_psites
        self.unique_psites += other.unique_psites
        self.total_reads += other.total_reads
        self.unique_reads += other.unique_reads
        if not self.sequence and other.sequence:
            self.sequence = other.sequence
        # Propagate CDS overlap flag
        if other.is_cds_overlap:
            self.is_cds_overlap = True
        # Merge overlapping genes (union)
        for gid in other.overlapping_gene_ids:
            if gid not in self.overlapping_gene_ids:
                self.overlapping_gene_ids.append(gid)
        if not self.start_codon and other.start_codon:
            self.start_codon = other.start_codon
        # Merge subset ORFs
        self.subset_orfs.extend(other.subset_orfs)
    
    def add_subset_orf(self, other):
        """
        Add another ORF as a subset of this representative.
        Records the subset's coordinates and source information.
        """
        blocks_str = ",".join(f"{s}-{e}" for s, e in other.blocks)
        tools = ",".join(sorted(set(t for t, s in other.sources)))
        samples = ",".join(sorted(set(s for t, s in other.sources)))
        self.subset_orfs.append({
            'blocks': blocks_str,
            'start': other.start,
            'end': other.end,
            'length_aa': other.length_aa,
            'tools': tools,
            'samples': samples,
            'tool_scores': dict(other.tool_scores),
            'tool_pvalues': dict(other.tool_pvalues)
        })
        # Merge sources into representative
        self.sources.update(other.sources)
        if other.is_cds_overlap:
            self.is_cds_overlap = True
        for gid in other.overlapping_gene_ids:
            if gid not in self.overlapping_gene_ids:
                self.overlapping_gene_ids.append(gid)
        # Merge tool scores and pvalues
        for tool, score in other.tool_scores.items():
            if tool not in self.tool_scores or (score is not None and self.tool_scores.get(tool) is None):
                self.tool_scores[tool] = score
        for tool, pval in other.tool_pvalues.items():
            if tool not in self.tool_pvalues or (pval is not None and self.tool_pvalues.get(tool) is None):
                self.tool_pvalues[tool] = pval
    
    @property
    def pN(self):
        """P-sites per nucleotide"""
        return self.total_psites / self.length_nt if self.length_nt > 0 else 0
    
    @property
    def unique_pN(self):
        """Unique P-sites per nucleotide"""
        return self.unique_psites / self.length_nt if self.length_nt > 0 else 0

# Parsers
def parse_ribotish(file_path, gtf_index, sample_id, min_len=0,
                   exclude_tistypes=None, atg_only=False):
    candidates = []
    print(f"Parsing Ribo-TISH: {file_path}", file=sys.stderr)
    
    # Check for placeholder file
    try:
        with _open(file_path, 'r') as f:
            first_line = f.readline().strip()
            if first_line.startswith('#') and ('placeholder' in first_line.lower() or 'no orfs' in first_line.lower() or 'insufficient' in first_line.lower()):
                print(f"  -> Placeholder file detected, skipping", file=sys.stderr)
                return []
    except Exception:
        pass
    
    filtered_tis = 0
    filtered_atg = 0
    
    try:
        with _open(file_path, 'r') as f:
            header = f.readline().strip().split('\t')
            col_map = {name: i for i, name in enumerate(header)}
            
            # Check for required columns - support two formats:
            # Format 1: Start, Stop columns (separate)
            # Format 2: GenomePos column (chr:start-end:strand)
            has_start_stop = 'Tid' in col_map and 'Start' in col_map and 'Stop' in col_map
            has_genomepos = 'Tid' in col_map and 'GenomePos' in col_map
            
            if not has_start_stop and not has_genomepos:
                print("Warning: Ribo-TISH file missing required columns (need Start/Stop or GenomePos). Skipping.", file=sys.stderr)
                return []

            for line in f:
                parts = line.strip().split('\t')
                tid = parts[col_map['Tid']]
                
                # Stage 0: TisType filtering
                if exclude_tistypes and 'TisType' in col_map:
                    tis_type = parts[col_map['TisType']] if col_map['TisType'] < len(parts) else ''
                    if tis_type in exclude_tistypes:
                        filtered_tis += 1
                        continue
                
                # Stage 0: ATG-only filtering
                if atg_only and 'StartCodon' in col_map:
                    start_codon = parts[col_map['StartCodon']] if col_map['StartCodon'] < len(parts) else ''
                    if start_codon.upper() != 'ATG':
                        filtered_atg += 1
                        continue
                
                # Extract TisPvalue as score (lower is better, convert to -log10(p))
                score = None
                pvalue = None
                if 'TisPvalue' in col_map:
                    try:
                        pval = float(parts[col_map['TisPvalue']])
                        if pval > 0:
                            import math
                            score = -math.log10(pval)  # Convert to -log10(p), higher is better
                            pvalue = pval  # Keep original p-value
                    except (ValueError, OverflowError):
                        pass
                
                # Parse coordinates
                t_start = None
                t_stop = None
                chrom = None
                strand = None
                
                if has_start_stop:
                    # Format 1: Direct Start/Stop columns (transcriptomic coordinates)
                    try:
                        t_start = int(parts[col_map['Start']])
                        t_stop = int(parts[col_map['Stop']])
                    except ValueError:
                        continue
                    
                    # Convert to genomic coordinates via GTF
                    length_nt = t_stop - t_start
                    if length_nt // 3 < min_len: continue
                    
                    chrom, strand, blocks = gtf_index.get_genomic_blocks(tid, t_start, t_stop, feature_type='exon')
                    
                elif has_genomepos:
                    # Format 2: GenomePos (chr:start-end:strand) - ALREADY genomic coordinates
                    genome_pos = parts[col_map['GenomePos']]
                    import re
                    match = re.search(r'(.+):(\d+)-(\d+):([+-])', genome_pos)
                    if match:
                        chrom = gtf_index.resolve_chrom(match.group(1))
                        t_start = int(match.group(2))
                        t_stop = int(match.group(3))
                        strand = match.group(4)
                        
                        length_nt = t_stop - t_start + 1  # Inclusive
                        if length_nt // 3 < min_len: continue
                        
                        # GenomePos gives genomic coords directly - create blocks
                        blocks = [(t_start, t_stop)]
                    else:
                        continue
                else:
                    continue
                
                if blocks:
                    gid = gtf_index.gene_map.get(tid, 'NA')
                    cand = ORFCandidate(chrom, strand, blocks, tid, gid, 'Ribo-TISH', sample_id, score=score, pvalue=pvalue)
                    candidates.append(cand)
    except Exception as e:
        print(f"Error parsing Ribo-TISH file: {e}", file=sys.stderr)
    
    if filtered_tis or filtered_atg:
        print(f"  -> Filtered: {filtered_tis} by TisType, {filtered_atg} by start codon", file=sys.stderr)
    return candidates

def parse_ribotricer(file_path, gtf_index, sample_id, min_len=0,
                     exclude_tistypes=None, atg_only=False):
    candidates = []
    print(f"Parsing Ribotricer: {file_path}", file=sys.stderr)
    
    # Check for placeholder file
    try:
        with _open(file_path, 'r') as f:
            first_line = f.readline().strip()
            if first_line.startswith('#') and ('placeholder' in first_line.lower() or 'no orfs' in first_line.lower() or 'insufficient' in first_line.lower()):
                print(f"  -> Placeholder file detected, skipping", file=sys.stderr)
                return []
    except Exception:
        pass
    
    filtered_tis = 0
    filtered_atg = 0
    
    try:
        with _open(file_path, 'r') as f:
            header_line = f.readline()
            while header_line.startswith('#'): header_line = f.readline()
            header = header_line.strip().split('\t')
            col_map = {name: i for i, name in enumerate(header)}

            if 'transcript_id' not in col_map or 'ORF_ID' not in col_map:
                 print("Warning: Missing required columns in Ribotricer file.", file=sys.stderr)
                 return []
            
            for line in f:
                parts = line.strip().split('\t')
                tid = parts[col_map['transcript_id']]
                orf_id = parts[col_map['ORF_ID']]
                
                # Stage 0: ORF_type filtering (ribotricer uses lowercase)
                if exclude_tistypes and 'ORF_type' in col_map:
                    orf_type = parts[col_map['ORF_type']] if col_map['ORF_type'] < len(parts) else ''
                    if orf_type in exclude_tistypes:
                        filtered_tis += 1
                        continue
                
                # Stage 0: ATG-only filtering
                if atg_only and 'start_codon' in col_map:
                    start_codon = parts[col_map['start_codon']] if col_map['start_codon'] < len(parts) else ''
                    if start_codon.upper() != 'ATG':
                        filtered_atg += 1
                        continue
                
                # Extract phase_score (0-1, higher is better)
                score = None
                if 'phase_score' in col_map:
                    try:
                        score = float(parts[col_map['phase_score']])
                    except ValueError:
                        pass
                
                blocks = None
                chrom = parts[col_map.get('chrom', -1)] if 'chrom' in col_map else None
                if chrom:
                    chrom = gtf_index.resolve_chrom(chrom)
                strand = parts[col_map.get('strand', -1)] if 'strand' in col_map else None
                
                # Parsing logic for Ribotricer ORF_ID: tid_start_end or tid_start_end_length
                match = re.search(r'_(\d+)_(\d+)(?:_\d+)?$', orf_id)
                if match:
                    t_start_raw = int(match.group(1))
                    t_end_raw = int(match.group(2))
                    t_start = t_start_raw
                    t_end = t_end_raw
                    
                    # Convert 1-based to 0-based
                    t_start -= 1
                    
                    if t_end - t_start < min_len * 3: continue
                    
                    c, s, b = gtf_index.get_genomic_blocks(tid, t_start, t_end, feature_type='exon')
                    if b:
                        chrom, strand, blocks = c, s, b
                    elif chrom and strand and t_end_raw >= t_start_raw:
                        length_nt = t_end_raw - t_start_raw + 1
                        if length_nt // 3 >= min_len:
                            blocks = [(t_start_raw, t_end_raw)]
                
                if blocks:
                    gid = gtf_index.gene_map.get(tid, 'NA')
                    cand = ORFCandidate(chrom, strand, blocks, tid, gid, 'Ribotricer', sample_id, score=score)
                    candidates.append(cand)
    except Exception as e:
        print(f"Error parsing Ribotricer file: {e}", file=sys.stderr)

    if filtered_tis or filtered_atg:
        print(f"  -> Filtered: {filtered_tis} by ORF_type, {filtered_atg} by start codon", file=sys.stderr)
    return candidates

def parse_orfquant(file_path, gtf_index, sample_id, min_len=0,
                   exclude_tistypes=None, atg_only=False):
    candidates = []
    print(f"Parsing ORFquant: {file_path}", file=sys.stderr)
    
    # Check for placeholder file
    try:
        with _open(file_path, 'r') as f:
            first_line = f.readline().strip()
            if first_line.startswith('#') and ('placeholder' in first_line.lower() or 'no orfs' in first_line.lower() or 'insufficient' in first_line.lower()):
                print(f"  -> Placeholder file detected, skipping", file=sys.stderr)
                return []
    except Exception:
        pass
    
    try:
        current_orf = None
        current_blocks = []
        current_attrs = {}
        
        with _open(file_path, 'r') as f:
            for line in f:
                if line.startswith('#'): continue
                parts = line.strip().split('\t')
                if len(parts) < 9: continue
                
                feature = parts[2]
                if feature != 'CDS': continue

                attrs = {}
                for p in parts[8].split(';'):
                    p = p.strip()
                    if not p: continue
                    kv = p.split(' ')
                    attrs[kv[0]] = ' '.join(kv[1:]).strip('"')
                
                orf_id = attrs.get('ORF_id') or attrs.get('transcript_id')
                if not orf_id: continue
                
                if current_orf != orf_id:
                    if current_orf and current_blocks:
                        chrom = current_attrs.get('chrom')
                        strand = current_attrs.get('strand')
                        tid = current_attrs.get('transcript_id', 'NA')
                        gid = current_attrs.get('gene_id', 'NA')
                        
                        # Extract P_sites as score (integer, higher is better)
                        score = None
                        psites_count = 0
                        if 'P_sites' in current_attrs:
                            try:
                                psites_count = int(float(current_attrs['P_sites']))
                                score = float(current_attrs['P_sites'])
                            except ValueError:
                                pass
                        
                        # Extract unique P_sites if available
                        unique_psites_count = 0
                        if 'P_sites_uniq' in current_attrs:
                            try:
                                unique_psites_count = int(float(current_attrs['P_sites_uniq']))
                            except ValueError:
                                pass
                        
                        length_nt = sum(e-s+1 for s,e in current_blocks)
                        if length_nt // 3 >= min_len:
                            cand = ORFCandidate(chrom, strand, current_blocks, tid, gid, 'ORFquant', sample_id, score=score)
                            # Set P-site counts from ORFquant attributes
                            cand.total_psites = psites_count
                            cand.unique_psites = unique_psites_count if unique_psites_count > 0 else psites_count
                            candidates.append(cand)
                    
                    current_orf = orf_id
                    current_blocks = []
                    current_attrs = attrs
                    current_attrs['chrom'] = gtf_index.resolve_chrom(parts[0])
                    current_attrs['strand'] = parts[6]
                
                current_blocks.append((int(parts[3]), int(parts[4])))

        if current_orf and current_blocks:
            chrom = current_attrs.get('chrom')
            strand = current_attrs.get('strand')
            tid = current_attrs.get('transcript_id', 'NA')
            gid = current_attrs.get('gene_id', 'NA')
            
            # Extract P_sites as score
            score = None
            psites_count = 0
            if 'P_sites' in current_attrs:
                try:
                    psites_count = int(float(current_attrs['P_sites']))
                    score = float(current_attrs['P_sites'])
                except ValueError:
                    pass
            
            # Extract unique P_sites if available
            unique_psites_count = 0
            if 'P_sites_uniq' in current_attrs:
                try:
                    unique_psites_count = int(float(current_attrs['P_sites_uniq']))
                except ValueError:
                    pass
            
            length_nt = sum(e-s+1 for s,e in current_blocks)
            if length_nt // 3 >= min_len:
                cand = ORFCandidate(chrom, strand, current_blocks, tid, gid, 'ORFquant', sample_id, score=score)
                # Set P-site counts from ORFquant attributes
                cand.total_psites = psites_count
                cand.unique_psites = unique_psites_count if unique_psites_count > 0 else psites_count
                candidates.append(cand)

    except Exception as e:
        print(f"Error parsing ORFquant file: {e}", file=sys.stderr)
    
    return candidates

def parse_ribocode(file_path, gtf_index, sample_id, min_len=0,
                   exclude_tistypes=None, atg_only=False):
    candidates = []
    print(f"Parsing RiboCode: {file_path}", file=sys.stderr)

    try:
        with _open(file_path, 'r') as f:
            first_line = f.readline().strip()
            if first_line.startswith('#') and ('placeholder' in first_line.lower() or 'no orfs' in first_line.lower() or 'insufficient' in first_line.lower()):
                print(f"  -> Placeholder file detected, skipping", file=sys.stderr)
                return []
            if first_line.startswith('FAILED'):
                print(f"  -> Failure marker detected, skipping", file=sys.stderr)
                return []
    except Exception:
        pass

    def parse_attrs(attr_str):
        attrs = {}
        for item in attr_str.split(';'):
            item = item.strip()
            if not item:
                continue
            parts = item.split(' ', 1)
            if len(parts) == 2:
                attrs[parts[0]] = parts[1].strip().strip('"')
        return attrs

    def score_from_pvalue(pvalue):
        if pvalue is None:
            return None
        try:
            pval = float(pvalue)
            if pval > 0:
                import math
                return -math.log10(pval)
        except (TypeError, ValueError, OverflowError):
            return None
        return None

    def load_txt_metrics(path):
        metrics = {}
        if not path or not os.path.exists(path):
            return metrics
        try:
            with _open(path, 'r') as f:
                header = f.readline().strip().split('\t')
                col_map = {name: i for i, name in enumerate(header)}
                if 'ORF_ID' not in col_map:
                    return metrics
                for line in f:
                    if not line.strip() or line.startswith('#'):
                        continue
                    parts = line.rstrip('\n').split('\t')
                    if len(parts) <= col_map['ORF_ID']:
                        continue
                    orf_id = parts[col_map['ORF_ID']]
                    pvalue = None
                    for col in ('adjusted_pval', 'pval_combined'):
                        if col in col_map and col_map[col] < len(parts):
                            try:
                                pvalue = float(parts[col_map[col]])
                                break
                            except ValueError:
                                pass
                    metrics[orf_id] = {
                        'pvalue': pvalue,
                        'score': score_from_pvalue(pvalue),
                    }
        except Exception as e:
            print(f"Warning: could not read RiboCode sidecar metrics {path}: {e}", file=sys.stderr)
        return metrics

    def sidecar_txt_path(path):
        root, ext = os.path.splitext(path)
        if ext == '.gtf' or ext == '.bed':
            return root + '.txt'
        return path if ext == '.txt' else None

    metrics = load_txt_metrics(sidecar_txt_path(file_path))

    file_format = _detect_file_format(file_path)
    if file_format == 'gtf':
        grouped = {}
        try:
            with _open(file_path, 'r') as f:
                for line in f:
                    if line.startswith('#') or not line.strip():
                        continue
                    parts = line.rstrip('\n').split('\t')
                    if len(parts) < 9:
                        continue
                    attrs = parse_attrs(parts[8])
                    orf_id = attrs.get('orf_id') or attrs.get('ORF_ID')
                    if not orf_id:
                        continue
                    rec = grouped.setdefault(orf_id, {
                        'chrom': gtf_index.resolve_chrom(parts[0]),
                        'strand': parts[6],
                        'blocks': [],
                        'orf_span': None,
                        'attrs': attrs,
                    })
                    rec['attrs'].update(attrs)
                    feature = parts[2]
                    block = (int(parts[3]), int(parts[4]))
                    if feature == 'exon':
                        rec['blocks'].append(block)
                    elif feature == 'ORF':
                        rec['orf_span'] = block

            for orf_id, rec in grouped.items():
                blocks = rec['blocks'] or ([rec['orf_span']] if rec['orf_span'] else [])
                if not blocks:
                    continue
                length_nt = sum(e - s + 1 for s, e in blocks)
                if length_nt // 3 < min_len:
                    continue
                attrs = rec['attrs']
                tid = attrs.get('transcript_id', 'NA')
                gid = attrs.get('gene_id') or gtf_index.gene_map.get(tid, 'NA')
                metric = metrics.get(orf_id, {})
                cand = ORFCandidate(
                    rec['chrom'],
                    rec['strand'],
                    blocks,
                    tid,
                    gid,
                    'RiboCode',
                    sample_id,
                    score=metric.get('score'),
                    pvalue=metric.get('pvalue'),
                )
                candidates.append(cand)
        except Exception as e:
            print(f"Error parsing RiboCode GTF file: {e}", file=sys.stderr)
        return candidates

    try:
        with _open(file_path, 'r') as f:
            header = f.readline().strip().split('\t')
            col_map = {name: i for i, name in enumerate(header)}
            required = {'ORF_ID', 'chrom', 'strand', 'ORF_gstart', 'ORF_gstop'}
            if not required.issubset(col_map):
                # If GTF parsing already succeeded above, silently skip TSV fallback
                if candidates:
                    return candidates
                print("Warning: RiboCode file missing required columns. Skipping.", file=sys.stderr)
                return []
            for line in f:
                if not line.strip() or line.startswith('#'):
                    continue
                parts = line.rstrip('\n').split('\t')
                orf_id = parts[col_map['ORF_ID']]
                chrom = gtf_index.resolve_chrom(parts[col_map['chrom']])
                strand = parts[col_map['strand']]
                gstart = int(parts[col_map['ORF_gstart']])
                gstop = int(parts[col_map['ORF_gstop']])
                start, end = sorted((gstart, gstop))
                if (end - start + 1) // 3 < min_len:
                    continue
                tid = parts[col_map['transcript_id']] if 'transcript_id' in col_map and col_map['transcript_id'] < len(parts) else 'NA'
                gid = parts[col_map['gene_id']] if 'gene_id' in col_map and col_map['gene_id'] < len(parts) else gtf_index.gene_map.get(tid, 'NA')
                pvalue = None
                for col in ('adjusted_pval', 'pval_combined'):
                    if col in col_map and col_map[col] < len(parts):
                        try:
                            pvalue = float(parts[col_map[col]])
                            break
                        except ValueError:
                            pass
                cand = ORFCandidate(
                    chrom,
                    strand,
                    [(start, end)],
                    tid,
                    gid,
                    'RiboCode',
                    sample_id,
                    score=score_from_pvalue(pvalue),
                    pvalue=pvalue,
                )
                candidates.append(cand)
    except Exception as e:
        print(f"Error parsing RiboCode file: {e}", file=sys.stderr)

    return candidates

def extract_sequence(cand, genome_fasta):
    """Extract nucleotide + amino acid sequence for a candidate ORF.
    Sets cand.sequence (nt) and cand.aa_sequence (aa).
    """
    try:
        seq_parts = []
        for s, e in cand.blocks:
            seq = genome_fasta[cand.chrom][s-1:e]
            seq_parts.append(str(seq))
        full_seq = "".join(seq_parts)
        seq_obj = Seq(full_seq)
        if cand.strand == '-':
            seq_obj = seq_obj.reverse_complement()
        cand.sequence = str(seq_obj)
        # Translate to AA (stop codon may or may not be present)
        aa = str(seq_obj.translate(to_stop=False))
        cand.aa_sequence = aa
        # Set start codon from first 3 nt of coding sequence
        if len(cand.sequence) >= 3:
            cand.start_codon = cand.sequence[:3].upper()
    except Exception:
        cand.sequence = "N" * cand.length_nt
        cand.aa_sequence = ""
        cand.start_codon = "UNK"

def validate_sequence(cand, genome_fasta):
    """Legacy wrapper – kept for back-compat; sequences already extracted."""
    if not cand.sequence:
        extract_sequence(cand, genome_fasta)


def annotate_cds_overlap(candidates, gtf_index):
    """
    Annotate each candidate with CDS overlap status and overlapping gene IDs.
    Must be called after extract_sequence() so that frame info is available.

    Sets:
      cand.is_cds_overlap      - True if ORF overlaps in-frame with annotated CDS
      cand.overlapping_gene_ids - list of gene_ids overlapping this ORF (same strand)
    """
    for cand in candidates:
        cand.overlapping_gene_ids = gtf_index.find_overlapping_genes(
            cand.chrom, cand.strand, cand.start, cand.end
        )
        cand.is_cds_overlap = gtf_index.find_cds_overlap_inframe(
            cand.chrom, cand.strand, cand.start, cand.end, cand.frame
        )


class BedgraphIndex:
    """
    Pre-loaded and indexed bedgraph data for fast region queries.
    Uses a simple binning strategy for efficient overlap queries.
    """
    BIN_SIZE = 10000  # 10kb bins
    
    def __init__(self, bedgraph_file):
        """Load bedgraph file into memory with binned index."""
        self.data = defaultdict(lambda: defaultdict(list))  # chrom -> bin -> [(start, end, value), ...]
        self.loaded = False
        
        if not bedgraph_file or not os.path.exists(bedgraph_file):
            return
        
        try:
            with _open(bedgraph_file, 'r') as f:
                for line in f:
                    if line.startswith('track') or line.startswith('#'):
                        continue
                    parts = line.strip().split('\t')
                    if len(parts) < 4:
                        continue
                    
                    chrom = parts[0]
                    start = int(parts[1])  # 0-based
                    end = int(parts[2])    # 0-based, exclusive
                    value = float(parts[3])
                    
                    # Add to all overlapping bins
                    start_bin = start // self.BIN_SIZE
                    end_bin = (end - 1) // self.BIN_SIZE
                    for bin_idx in range(start_bin, end_bin + 1):
                        self.data[chrom][bin_idx].append((start, end, value))
            
            self.loaded = True
        except Exception as e:
            print(f"Warning: Error loading bedgraph: {e}", file=sys.stderr)
    
    def count_in_region(self, chrom, start_1based, end_1based):
        """Count signal in a genomic region (1-based coordinates)."""
        if not self.loaded:
            return 0
        
        # Convert to 0-based
        start_0 = start_1based - 1
        end_0 = end_1based
        
        total_count = 0
        start_bin = start_0 // self.BIN_SIZE
        end_bin = (end_0 - 1) // self.BIN_SIZE
        
        seen = set()  # Avoid double-counting entries in multiple bins
        
        for bin_idx in range(start_bin, end_bin + 1):
            for entry in self.data[chrom].get(bin_idx, []):
                entry_id = id(entry)
                if entry_id in seen:
                    continue
                seen.add(entry_id)
                
                bg_start, bg_end, bg_value = entry
                
                # Calculate overlap
                overlap_start = max(bg_start, start_0)
                overlap_end = min(bg_end, end_0)
                
                if overlap_start < overlap_end:
                    overlap_len = overlap_end - overlap_start
                    total_count += bg_value * overlap_len
        
        return int(total_count)


def load_bedgraph_indices(bedgraph_dir, sample_list, metric_types=None):
    """
    Pre-load all bedgraph files into indexed structures.
    Returns: dict of sample -> strand -> type -> BedgraphIndex
    """
    if not bedgraph_dir or not os.path.exists(bedgraph_dir):
        return None

    if metric_types is None:
        metric_types = ('psite', 'psite_uniq', 'coverage', 'coverage_uniq')
    
    print(f"Pre-loading bedgraph files for {len(sample_list)} samples...", file=sys.stderr)
    indices = {}
    
    for sample in sample_list:
        indices[sample] = {}
        for strand_suffix in ['plus', 'minus']:
            indices[sample][strand_suffix] = {}
            
            file_map = {
                'psite': os.path.join(bedgraph_dir, f"{sample}_P_sites_{strand_suffix}.bedgraph"),
                'psite_uniq': os.path.join(bedgraph_dir, f"{sample}_P_sites_uniq_{strand_suffix}.bedgraph"),
                'coverage': os.path.join(bedgraph_dir, f"{sample}_coverage_{strand_suffix}.bedgraph"),
                'coverage_uniq': os.path.join(bedgraph_dir, f"{sample}_coverage_uniq_{strand_suffix}.bedgraph"),
            }
            for metric in metric_types:
                indices[sample][strand_suffix][metric] = BedgraphIndex(file_map[metric])
    
    print(f"Bedgraph indices loaded.", file=sys.stderr)
    return indices


def count_psites_in_region(bedgraph_file, chrom, start, end):
    """
    Count P-sites in a genomic region from bedgraph file
    bedgraph format: chrom start end value
    Returns: total count
    
    NOTE: This function is kept for backward compatibility but is slow.
    Use BedgraphIndex for better performance.
    """
    if not os.path.exists(bedgraph_file):
        return 0
    
    total_count = 0
    try:
        with _open(bedgraph_file, 'r') as f:
            for line in f:
                if line.startswith('track') or line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) < 4:
                    continue
                
                bg_chrom = parts[0]
                bg_start = int(parts[1])  # 0-based
                bg_end = int(parts[2])    # 0-based, exclusive
                bg_value = float(parts[3])
                
                if bg_chrom != chrom:
                    continue
                
                # Convert ORF coords (1-based) to 0-based for overlap check
                orf_start_0 = start - 1
                orf_end_0 = end
                
                # Calculate overlap
                overlap_start = max(bg_start, orf_start_0)
                overlap_end = min(bg_end, orf_end_0)
                
                if overlap_start < overlap_end:
                    overlap_len = overlap_end - overlap_start
                    total_count += bg_value * overlap_len
    
    except Exception as e:
        print(f"Warning: Error reading bedgraph {bedgraph_file}: {e}", file=sys.stderr)
        return 0
    
    return int(total_count)

def calculate_statistics_from_bedgraphs(cand, bedgraph_dir, sample_list):
    """
    Calculate statistics from RiboseQC bedgraph files (slow legacy version)
    """
    if not bedgraph_dir or not os.path.exists(bedgraph_dir):
        return
    
    total_psites = 0
    unique_psites = 0
    total_reads = 0
    unique_reads = 0
    
    strand_suffix = 'plus' if cand.strand == '+' else 'minus'
    
    for sample in sample_list:
        # P-site bedgraphs
        psite_file = os.path.join(bedgraph_dir, f"{sample}_P_sites_{strand_suffix}.bedgraph")
        psite_uniq_file = os.path.join(bedgraph_dir, f"{sample}_P_sites_uniq_{strand_suffix}.bedgraph")
        
        # Coverage bedgraphs (optional)
        coverage_file = os.path.join(bedgraph_dir, f"{sample}_coverage_{strand_suffix}.bedgraph")
        coverage_uniq_file = os.path.join(bedgraph_dir, f"{sample}_coverage_uniq_{strand_suffix}.bedgraph")
        
        # Sum across all exon blocks
        for block_start, block_end in cand.blocks:
            total_psites += count_psites_in_region(psite_file, cand.chrom, block_start, block_end)
            unique_psites += count_psites_in_region(psite_uniq_file, cand.chrom, block_start, block_end)
            total_reads += count_psites_in_region(coverage_file, cand.chrom, block_start, block_end)
            unique_reads += count_psites_in_region(coverage_uniq_file, cand.chrom, block_start, block_end)
    
    cand.total_psites = total_psites
    cand.unique_psites = unique_psites
    cand.total_reads = total_reads
    cand.unique_reads = unique_reads


def calculate_statistics_from_indices(cand, bedgraph_indices, sample_list, metric_types=None):
    """
    Calculate statistics using pre-loaded bedgraph indices (fast version)
    """
    if not bedgraph_indices:
        return

    if metric_types is None:
        metric_types = ('psite', 'psite_uniq', 'coverage', 'coverage_uniq')
    
    total_psites = 0
    unique_psites = 0
    total_reads = 0
    unique_reads = 0
    
    strand_suffix = 'plus' if cand.strand == '+' else 'minus'
    
    for sample in sample_list:
        if sample not in bedgraph_indices:
            continue
        
        idx = bedgraph_indices[sample][strand_suffix]
        
        # Sum across all exon blocks
        for block_start, block_end in cand.blocks:
            if 'psite' in metric_types:
                total_psites += idx['psite'].count_in_region(cand.chrom, block_start, block_end)
            if 'psite_uniq' in metric_types:
                unique_psites += idx['psite_uniq'].count_in_region(cand.chrom, block_start, block_end)
            if 'coverage' in metric_types:
                total_reads += idx['coverage'].count_in_region(cand.chrom, block_start, block_end)
            if 'coverage_uniq' in metric_types:
                unique_reads += idx['coverage_uniq'].count_in_region(cand.chrom, block_start, block_end)
    
    cand.total_psites = total_psites
    cand.unique_psites = unique_psites
    cand.total_reads = total_reads
    cand.unique_reads = unique_reads


def process_candidate_batch(batch_data):
    """
    Process a batch of candidates for P-site statistics.
    Used for parallel processing.
    
    batch_data: tuple of (candidate_dicts, bedgraph_indices, sample_list)
    Returns: list of updated candidate dicts with statistics
    """
    candidates, bedgraph_indices, sample_list, metric_types = batch_data
    
    results = []
    for cand_dict in candidates:
        # Reconstruct minimal candidate info for counting
        chrom = cand_dict['chrom']
        strand = cand_dict['strand']
        blocks = cand_dict['blocks']
        strand_suffix = 'plus' if strand == '+' else 'minus'
        
        total_psites = 0
        unique_psites = 0
        total_reads = 0
        unique_reads = 0
        
        for sample in sample_list:
            if sample not in bedgraph_indices:
                continue
            
            idx = bedgraph_indices[sample][strand_suffix]
            
            for block_start, block_end in blocks:
                if 'psite' in metric_types:
                    total_psites += idx['psite'].count_in_region(chrom, block_start, block_end)
                if 'psite_uniq' in metric_types:
                    unique_psites += idx['psite_uniq'].count_in_region(chrom, block_start, block_end)
                if 'coverage' in metric_types:
                    total_reads += idx['coverage'].count_in_region(chrom, block_start, block_end)
                if 'coverage_uniq' in metric_types:
                    unique_reads += idx['coverage_uniq'].count_in_region(chrom, block_start, block_end)
        
        results.append({
            'id_key': cand_dict['id_key'],
            'total_psites': total_psites,
            'unique_psites': unique_psites,
            'total_reads': total_reads,
            'unique_reads': unique_reads
        })
    
    return results


def calculate_statistics_parallel(final_list, bedgraph_indices, sample_list, num_workers=None, metric_types=None):
    """
    Calculate P-site statistics for all candidates in parallel.
    """
    if not bedgraph_indices or not sample_list:
        return

    if metric_types is None:
        metric_types = ('psite', 'psite_uniq', 'coverage', 'coverage_uniq')
    
    if num_workers is None:
        num_workers = min(multiprocessing.cpu_count(), 8)
    
    # Convert candidates to serializable dicts
    cand_dicts = [
        {
            'id_key': cand.id_key,
            'chrom': cand.chrom,
            'strand': cand.strand,
            'blocks': cand.blocks
        }
        for cand in final_list
    ]
    
    # Create a lookup for updating results
    cand_lookup = {cand.id_key: cand for cand in final_list}
    
    # For small datasets or single worker, process sequentially
    if len(final_list) < 100 or num_workers <= 1:
        print(f"Processing {len(final_list)} candidates sequentially...", file=sys.stderr)
        for cand in final_list:
            calculate_statistics_from_indices(cand, bedgraph_indices, sample_list, metric_types=metric_types)
        return
    
    # Split into batches for parallel processing
    batch_size = max(50, len(cand_dicts) // num_workers)
    batches = [cand_dicts[i:i + batch_size] for i in range(0, len(cand_dicts), batch_size)]
    
    print(f"Processing {len(final_list)} candidates in {len(batches)} batches using {num_workers} workers...", file=sys.stderr)
    
    # Note: Since BedgraphIndex contains complex data structures that are hard to pickle,
    # we'll use threading instead of multiprocessing for simplicity
    from concurrent.futures import ThreadPoolExecutor
    
    with ThreadPoolExecutor(max_workers=num_workers) as executor:
        futures = []
        for batch in batches:
            futures.append(executor.submit(process_candidate_batch, (batch, bedgraph_indices, sample_list, metric_types)))
        
        processed = 0
        for future in as_completed(futures):
            results = future.result()
            for res in results:
                cand = cand_lookup.get(res['id_key'])
                if cand:
                    cand.total_psites = res['total_psites']
                    cand.unique_psites = res['unique_psites']
                    cand.total_reads = res['total_reads']
                    cand.unique_reads = res['unique_reads']
            processed += len(results)
            print(f"  Processed {processed}/{len(final_list)} candidates...", file=sys.stderr, end='\r')
    
    print(f"  Processed {len(final_list)}/{len(final_list)} candidates.", file=sys.stderr)


BEDGRAPH_METRIC_FILES = {
    'psite': 'P_sites',
    'psite_uniq': 'P_sites_uniq',
    'coverage': 'coverage',
    'coverage_uniq': 'coverage_uniq',
}


def _build_candidate_block_bins(final_list, bin_size=BedgraphIndex.BIN_SIZE):
    """
    Build a low-memory block index for streaming bedgraph updates.
    Returns: strand_suffix -> chrom -> bin -> [(cand_idx, block_start, block_end), ...]
    """
    block_bins = {
        'plus': defaultdict(lambda: defaultdict(list)),
        'minus': defaultdict(lambda: defaultdict(list)),
    }
    for cand_idx, cand in enumerate(final_list):
        strand_suffix = 'plus' if cand.strand == '+' else 'minus'
        for block_start, block_end in cand.blocks:
            record = (cand_idx, block_start, block_end)
            start_bin = (block_start - 1) // bin_size
            end_bin = (block_end - 1) // bin_size
            for bin_idx in range(start_bin, end_bin + 1):
                block_bins[strand_suffix][cand.chrom][bin_idx].append(record)
    return block_bins


def _stream_bedgraph_file_counts(task):
    """
    Stream a single bedgraph file and accumulate counts for candidate blocks.

    Returns: (counts, max_vals, val_sums, val_counts)
        counts      - defaultdict(int):  weighted sum (bg_value * overlap_bp)
        max_vals    - dict:              max bg_value per cand_idx (for pN)
        val_sums    - defaultdict(float): sum of bg_values per cand_idx (for pN)
        val_counts  - defaultdict(int):   count of overlapping intervals per cand_idx (for pN)
    """
    bedgraph_file, chrom_bins = task
    counts = defaultdict(int)
    max_vals = {}
    val_sums = defaultdict(float)
    val_counts = defaultdict(int)

    if not bedgraph_file or not os.path.exists(bedgraph_file):
        return (counts, max_vals, val_sums, val_counts)

    try:
        with _open(bedgraph_file, 'r') as handle:
            for line in handle:
                if line.startswith('track') or line.startswith('#'):
                    continue
                parts = line.rstrip('\n').split('\t')
                if len(parts) < 4:
                    continue

                chrom = parts[0]
                if chrom not in chrom_bins:
                    continue

                bg_start = int(parts[1])      # 0-based inclusive
                bg_end = int(parts[2])        # 0-based exclusive
                bg_value = float(parts[3])
                if bg_value == 0:
                    continue

                start_bin = bg_start // BedgraphIndex.BIN_SIZE
                end_bin = (bg_end - 1) // BedgraphIndex.BIN_SIZE
                seen = set()
                bins_for_chrom = chrom_bins[chrom]

                for bin_idx in range(start_bin, end_bin + 1):
                    for cand_idx, block_start, block_end in bins_for_chrom.get(bin_idx, []):
                        block_key = (cand_idx, block_start, block_end)
                        if block_key in seen:
                            continue
                        seen.add(block_key)

                        block_start0 = block_start - 1
                        overlap_start = max(bg_start, block_start0)
                        overlap_end = min(bg_end, block_end)
                        if overlap_start < overlap_end:
                            counts[cand_idx] += int(bg_value * (overlap_end - overlap_start))
                            # Per-interval stats for pN calculation
                            max_vals[cand_idx] = max(max_vals.get(cand_idx, bg_value), bg_value)
                            val_sums[cand_idx] += bg_value
                            val_counts[cand_idx] += 1
    except Exception as exc:
        print(f"Warning: Error streaming bedgraph {bedgraph_file}: {exc}", file=sys.stderr)

    return (counts, max_vals, val_sums, val_counts)


def _estimate_bedgraph_bytes(bedgraph_dir, sample_list, metric_types):
    total = 0
    for sample in sample_list:
        for strand_suffix in ('plus', 'minus'):
            for metric in metric_types:
                path = os.path.join(
                    bedgraph_dir,
                    f"{sample}_{BEDGRAPH_METRIC_FILES[metric]}_{strand_suffix}.bedgraph"
                )
                if os.path.exists(path):
                    total += os.path.getsize(path)
    return total


def calculate_statistics_streaming(final_list, bedgraph_dir, sample_list, num_workers=None, metric_types=None):
    """
    Low-memory statistics path: stream each bedgraph file once and update all ORFs.

    Accumulates both aggregate stats (on ORFCandidate objects) and per-sample
    breakdowns returned as a dict suitable for writing expression output files.
    """
    if not bedgraph_dir or not os.path.exists(bedgraph_dir) or not sample_list:
        return None

    if metric_types is None:
        metric_types = ('psite', 'psite_uniq', 'coverage', 'coverage_uniq')
    if num_workers is None:
        num_workers = min(multiprocessing.cpu_count(), 4)

    for cand in final_list:
        cand.total_psites = 0
        cand.unique_psites = 0
        cand.total_reads = 0
        cand.unique_reads = 0

    block_bins = _build_candidate_block_bins(final_list)
    attr_map = {
        'psite': 'total_psites',
        'psite_uniq': 'unique_psites',
        'coverage': 'total_reads',
        'coverage_uniq': 'unique_reads',
    }

    # --- Per-sample tracking setup ---
    n_orfs = len(final_list)
    n_samples = len(sample_list)
    sample_idx = {s: i for i, s in enumerate(sample_list)}
    per_sample_reads = [[0] * n_samples for _ in range(n_orfs)]
    per_sample_max = [[0.0] * n_samples for _ in range(n_orfs)]
    per_sample_val_sum = [[0.0] * n_samples for _ in range(n_orfs)]
    per_sample_val_cnt = [[0] * n_samples for _ in range(n_orfs)]
    per_sample_coverage = [[0.0] * n_samples for _ in range(n_orfs)]

    tasks = []
    for sample in sample_list:
        s_idx = sample_idx[sample]
        for strand_suffix in ('plus', 'minus'):
            chrom_bins = block_bins[strand_suffix]
            if not chrom_bins:
                continue
            for metric in metric_types:
                bedgraph_file = os.path.join(
                    bedgraph_dir,
                    f"{sample}_{BEDGRAPH_METRIC_FILES[metric]}_{strand_suffix}.bedgraph"
                )
                tasks.append((metric, bedgraph_file, chrom_bins, s_idx))

    if not tasks:
        return None

    max_workers = max(1, min(num_workers, len(tasks), 4))
    print(
        f"Streaming {len(tasks)} bedgraph files with {max_workers} worker(s)...",
        file=sys.stderr
    )

    completed = 0
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        future_map = {}
        for metric, bedgraph_file, chrom_bins, s_idx in tasks:
            future = executor.submit(_stream_bedgraph_file_counts, (bedgraph_file, chrom_bins))
            future_map[future] = (metric, s_idx)

        for future in as_completed(future_map):
            metric, s_idx = future_map[future]
            counts, max_vals, val_sums, val_counts = future.result()

            # Aggregate stats (unchanged behaviour)
            attr = attr_map[metric]
            for cand_idx, value in counts.items():
                setattr(final_list[cand_idx], attr, getattr(final_list[cand_idx], attr) + value)

            # Per-sample accumulation
            if metric in ('psite', 'psite_uniq'):
                for cand_idx, value in counts.items():
                    per_sample_reads[cand_idx][s_idx] += value
                for cand_idx, val in max_vals.items():
                    per_sample_max[cand_idx][s_idx] = max(per_sample_max[cand_idx][s_idx], val)
                for cand_idx, val in val_sums.items():
                    per_sample_val_sum[cand_idx][s_idx] += val
                for cand_idx, val in val_counts.items():
                    per_sample_val_cnt[cand_idx][s_idx] += val
            elif metric in ('coverage', 'coverage_uniq'):
                for cand_idx, value in counts.items():
                    per_sample_coverage[cand_idx][s_idx] += value

            completed += 1
            print(f"  Streamed {completed}/{len(tasks)} bedgraph files...", file=sys.stderr, end='\r')

    print(f"  Streamed {len(tasks)}/{len(tasks)} bedgraph files.", file=sys.stderr)

    # Compute derived per-sample values
    per_sample_pN = [[0.0] * n_samples for _ in range(n_orfs)]
    per_sample_rpkm = [[0.0] * n_samples for _ in range(n_orfs)]
    for i in range(n_orfs):
        length_kb = final_list[i].length_nt / 1000.0
        for s in range(n_samples):
            if per_sample_reads[i][s] > 0 and per_sample_val_sum[i][s] > 0:
                per_sample_pN[i][s] = round(
                    per_sample_max[i][s] * per_sample_val_cnt[i][s] / max(per_sample_val_sum[i][s], 1e-10), 4
                )
            if length_kb > 0:
                per_sample_rpkm[i][s] = per_sample_coverage[i][s] / length_kb

    # Compute TPM (per-sample normalization)
    per_sample_tpm = [[0.0] * n_samples for _ in range(n_orfs)]
    for s in range(n_samples):
        total_rpkm = sum(per_sample_rpkm[i][s] for i in range(n_orfs))
        if total_rpkm > 0:
            inv_total = 1e6 / total_rpkm
            for i in range(n_orfs):
                per_sample_tpm[i][s] = per_sample_rpkm[i][s] * inv_total

    return {
        'per_sample_reads': per_sample_reads,
        'per_sample_pN': per_sample_pN,
        'per_sample_coverage_rpm': per_sample_coverage,
        'per_sample_rpkm': per_sample_rpkm,
        'per_sample_tpm': per_sample_tpm,
        'sample_idx': sample_idx,
        'sample_list': sample_list,
    }


def parse_price(file_path, gtf_index, sample_id, min_len=0,
                exclude_tistypes=None, atg_only=False):
    """Parse PRICE (GEDI) ORF TSV output.

    PRICE TSV format:
      Gene  Id  Location  Candidate Location  Codon  Type  Start  Range  p_value  [conditions...]  Total

    Location column: aa_start-aa_stop:chr:genomic_start-genomic_end:strand
    Candidate Location: chr:start-end|start-end:strand (multi-exon if splice junctions)
    """
    candidates = []
    if exclude_tistypes is None:
        exclude_tistypes = set()
    elif isinstance(exclude_tistypes, str):
        exclude_tistypes = {t.strip() for t in exclude_tistypes.split(',')}

    # Route GTF files through the orfquant parser
    if _detect_file_format(file_path) == 'gtf':
        return parse_orfquant(file_path, gtf_index, sample_id,
                              min_len=min_len, exclude_tistypes=exclude_tistypes,
                              atg_only=atg_only)

    try:
        with open(file_path, 'r') as fh:
            header = None
            for line in fh:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                parts = line.split('\t')
                if header is None:
                    header = parts
                    continue
                if len(parts) < 9:
                    continue

                gene_id = parts[0] if parts[0] else 'unknown'
                orf_id = parts[1]
                loc_str = parts[2]
                cand_loc_str = parts[3] if len(parts) > 3 else ''
                start_codon = parts[4] if len(parts) > 4 else 'ATG'
                orf_type = parts[5] if len(parts) > 5 else 'unknown'
                range_score_str = parts[7] if len(parts) > 7 else '0'
                pvalue_str = parts[8] if len(parts) > 8 else '1'

                if exclude_tistypes and orf_type in exclude_tistypes:
                    continue
                if atg_only and start_codon.upper() != 'ATG':
                    continue

                if min_len > 1:
                    try:
                        aa_part = loc_str.split(':')[0]
                        aa_start, aa_stop = aa_part.split('-')
                        if int(aa_stop) - int(aa_start) + 1 < min_len:
                            continue
                    except (ValueError, IndexError):
                        pass

                chrom, strand, blocks = None, '+', []
                try:
                    colon_parts = loc_str.split(':')
                    if len(colon_parts) >= 3:
                        # Format A: aa_start-aa_stop:chr:start-end:strand (old)
                        chrom = colon_parts[1]
                        coords = colon_parts[2]
                        strand = colon_parts[3] if len(colon_parts) > 3 else '+'
                        blocks = [(int(coords.split('-')[0]), int(coords.split('-')[1]))]
                    elif len(colon_parts) == 2:
                        # Format B: chr+strand:start-end (actual PRICE output)
                        # chr+strand part: e.g. "1+" or "MT-"
                        chr_strand = colon_parts[0]
                        coords = colon_parts[1]
                        # Extract strand from last character
                        if chr_strand[-1] in ('+', '-'):
                            strand = chr_strand[-1]
                            chrom = chr_strand[:-1]
                        else:
                            strand = '+'
                            chrom = chr_strand
                        # Parse coordinates (may have | for multi-exon)
                        if '|' in coords:
                            exons = coords.split('|')
                            parsed = []
                            for e in exons:
                                se = e.split('-')
                                parsed.append((int(se[0]), int(se[1])))
                            if parsed:
                                blocks = sorted(parsed)
                        else:
                            blocks = [(int(coords.split('-')[0]), int(coords.split('-')[1]))]
                except (ValueError, IndexError):
                    continue

                if chrom is None:
                    continue

                # Also try Candidate Location for multi-exon blocks if not already multi
                if cand_loc_str and '|' in cand_loc_str and len(blocks) == 1:
                    try:
                        cp = cand_loc_str.split(':')
                        if len(cp) >= 2:
                            exons = cp[1].split('|')
                            parsed = []
                            for e in exons:
                                se = e.split('-')
                                parsed.append((int(se[0]), int(se[1])))
                            if parsed:
                                blocks = sorted(parsed)
                    except (ValueError, IndexError):
                        pass

                try:
                    pvalue = float(pvalue_str)
                except ValueError:
                    pvalue = 1.0
                try:
                    score = float(range_score_str)
                except ValueError:
                    score = None

                tid = f"PRICE_{gene_id}_{orf_id}"
                cand = ORFCandidate(chrom, strand, blocks, tid, gene_id,
                                  'PRICE', sample_id, score=score)
                candidates.append(cand)

    except Exception as e:
        print(f"Error parsing PRICE file {file_path}: {e}", file=sys.stderr)

    return candidates


# ---------------------------------------------------------------------------
# Top-level functions for parallel workers (must be at module level so that
# multiprocessing can pickle/reference them).
# ---------------------------------------------------------------------------

def _parse_file_task(task):
    """Parallel worker: parse one input file.

    Inherits *_shared_gtf_index* from the parent process via fork (Linux COW).
    No GTFIndex serialisation overhead.

    Args:
        task: (tool, file_path, sid, min_len, exclude_tistypes, atg_only)
    Returns:
        (tool, sid, orfs_list)
    """
    global _shared_gtf_index
    tool, file_path, sid, min_len, exclude_tistypes, atg_only = task
    _parse_fn = {
        'Ribo-TISH': parse_ribotish,
        'Ribotricer': parse_ribotricer,
        'RiboCode':   parse_ribocode,
        'ORFquant':   parse_orfquant,
        'PRICE':      parse_price,
    }[tool]
    orfs = _parse_fn(file_path, _shared_gtf_index, sid, min_len,
                     exclude_tistypes=exclude_tistypes, atg_only=atg_only)
    return tool, sid, orfs


def _extract_annotate_chunk(chunk):
    """Parallel worker: sequence extraction + CDS annotation for a batch of ORFs.

    Inherits *_shared_gtf_index* via fork; creates its own ``pyfaidx.Fasta``
    instance from *_shared_fasta_path* to avoid shared file-handle races.

    Args:
        chunk: list of (id_key, chrom, strand, blocks, frame, length_nt)
    Returns:
        list of (id_key, nt_seq, aa_seq, start_codon, is_cds_overlap, gene_ids)
    """
    global _shared_gtf_index, _shared_fasta_path
    try:
        local_fasta = Fasta(_shared_fasta_path)
    except Exception as exc:
        print(f"Worker: failed to open FASTA '{_shared_fasta_path}': {exc}",
              file=sys.stderr)
        return [(row[0], '', '', 'UNK', False, []) for row in chunk]

    results = []
    for id_key, chrom, strand, blocks, frame, length_nt in chunk:
        # --- Sequence extraction ---
        try:
            seq_parts = [str(local_fasta[chrom][s - 1:e]) for s, e in blocks]
            full_seq  = ''.join(seq_parts)
            seq_obj   = Seq(full_seq)
            if strand == '-':
                seq_obj = seq_obj.reverse_complement()
            nt_seq      = str(seq_obj)
            aa_seq      = str(seq_obj.translate(to_stop=False))
            start_codon = nt_seq[:3].upper() if len(nt_seq) >= 3 else ''
        except Exception:
            nt_seq      = 'N' * length_nt
            aa_seq      = ''
            start_codon = 'UNK'

        # --- CDS + gene overlap annotation ---
        try:
            orf_start = blocks[0][0]
            orf_end   = blocks[-1][1]
            gene_ids  = _shared_gtf_index.find_overlapping_genes(
                chrom, strand, orf_start, orf_end)
            is_cds    = _shared_gtf_index.find_cds_overlap_inframe(
                chrom, strand, orf_start, orf_end, frame)
        except Exception:
            gene_ids = []
            is_cds   = False

        results.append((id_key, nt_seq, aa_seq, start_codon, is_cds, gene_ids))
    return results


def annotate_sequences_and_cds(candidates, use_parallel, num_workers, genome_fasta, gtf_index):
    """
    Populate sequence/CDS-related annotations, streaming chunk results to avoid
    holding a second full copy of all sequence strings in memory.
    """
    if not candidates:
        return

    if use_parallel and len(candidates) > 200:
        chunk_size = max(200, len(candidates) // (num_workers * 4))
        chunk_data = [
            (c.id_key, c.chrom, c.strand, c.blocks, c.frame, c.length_nt)
            for c in candidates
        ]
        chunks = [chunk_data[i:i + chunk_size] for i in range(0, len(chunk_data), chunk_size)]
        cand_lookup = {c.id_key: c for c in candidates}

        print(
            f"Extracting sequences + annotating CDS overlap "
            f"({len(candidates)} candidates, {num_workers} workers, {len(chunks)} chunks)...",
            file=sys.stderr
        )
        try:
            ctx = multiprocessing.get_context('fork')
            processed = 0
            with ctx.Pool(num_workers) as pool:
                for chunk_results in pool.imap_unordered(_extract_annotate_chunk, chunks):
                    for id_key, nt_seq, aa_seq, start_codon, is_cds, gene_ids in chunk_results:
                        cand = cand_lookup.get(id_key)
                        if cand:
                            cand.sequence = nt_seq
                            cand.aa_sequence = aa_seq
                            cand.start_codon = start_codon
                            cand.is_cds_overlap = is_cds
                            cand.overlapping_gene_ids = gene_ids
                    processed += len(chunk_results)
                    print(f"  Annotated {processed}/{len(candidates)} candidates...", file=sys.stderr, end='\r')
            print(f"  Annotated {len(candidates)}/{len(candidates)} candidates.", file=sys.stderr)
            return
        except Exception as exc:
            print(
                f"Parallel extract/annotate failed ({exc}), retrying sequentially...",
                file=sys.stderr
            )

    print(f"Extracting sequences ({len(candidates)} candidates)...", file=sys.stderr)
    for cand in candidates:
        extract_sequence(cand, genome_fasta)
    print("Annotating CDS overlap and overlapping genes...", file=sys.stderr)
    annotate_cds_overlap(candidates, gtf_index)


def _write_expression_outputs(candidates, per_sample_data, prefix):
    """Write per-ORF expression summary and RPKM/TPM TSV files.

    Columns match the output of quantify_orf_expression.py + calc_orf_rpkm_tpm.py
    exactly so downstream consumers (MultiQC, user scripts) see no difference.
    """
    if not per_sample_data or not candidates:
        return

    sample_list = per_sample_data['sample_list']
    sample_idx_map = per_sample_data['sample_idx']
    sorted_samples = sorted(sample_list, key=lambda s: sample_idx_map[s])
    n_orfs = len(candidates)

    # --- expression_summary.tsv ---
    summary_path = f"{prefix}_expression_summary.tsv"
    with open(summary_path, 'w') as f:
        header = ["orf_id", "chrom", "start", "end", "strand"]
        for s in sorted_samples:
            header.append(f"{s}_reads")
            header.append(f"{s}_pN")
        header.extend(["total_reads", "n_expressed_samples"])
        f.write('\t'.join(header) + '\n')

        for i, cand in enumerate(candidates):
            orf_id = f"ORF_{i+1}_{cand.gid}"
            row = [orf_id, cand.chrom, str(cand.start), str(cand.end), cand.strand]
            for s in sorted_samples:
                s_idx = sample_idx_map[s]
                reads_val = per_sample_data['per_sample_reads'][i][s_idx]
                pN_val = per_sample_data['per_sample_pN'][i][s_idx]
                row.append(str(reads_val))
                row.append(f"{pN_val:.4f}")
            total_r = sum(per_sample_data['per_sample_reads'][i])
            n_exp = sum(1 for v in per_sample_data['per_sample_reads'][i] if v > 0)
            row.extend([str(total_r), str(n_exp)])
            f.write('\t'.join(row) + '\n')

    # --- expression_rpkm_tpm.tsv ---
    rpkm_path = f"{prefix}_expression_rpkm_tpm.tsv"
    with open(rpkm_path, 'w') as f:
        header = ["orf_id", "chrom", "start", "end", "strand", "orf_length", "orf_length_kb"]
        for s in sorted_samples:
            header.append(f"{s}_coverage_rpm")
            header.append(f"{s}_rpkm")
            header.append(f"{s}_tpm")
        f.write('\t'.join(header) + '\n')

        for i, cand in enumerate(candidates):
            orf_id = f"ORF_{i+1}_{cand.gid}"
            row = [orf_id, cand.chrom, str(cand.start), str(cand.end), cand.strand,
                   str(cand.length_nt), f"{cand.length_nt / 1000.0:.3f}"]
            for s in sorted_samples:
                s_idx = sample_idx_map[s]
                row.append(f"{per_sample_data['per_sample_coverage_rpm'][i][s_idx]:.4f}")
                row.append(f"{per_sample_data['per_sample_rpkm'][i][s_idx]:.4f}")
                row.append(f"{per_sample_data['per_sample_tpm'][i][s_idx]:.4f}")
            f.write('\t'.join(row) + '\n')

    print(f"Expression outputs written to {summary_path}, {rpkm_path}", file=sys.stderr)


def _write_orf_outputs(candidates: list, prefix: str) -> None:
    """Write BED12, metadata TSV, and GTF files for *candidates* to files named *prefix*.{bed,metadata.tsv,gtf}."""
    with open(f"{prefix}.metadata.tsv", 'w') as out:
        header = ["orf_id", "chrom", "strand", "start", "end", "length_aa", "exon_blocks",
                  "gene_id", "transcript_id", "tools", "samples",
                  "tool_scores", "tool_pvalues",
                  "total_reads", "unique_reads", "total_psites", "unique_psites", "pN", "unique_pN",
                  "num_subset_orfs", "subset_orfs",
                  "sequence", "start_codon", "aa_sequence", "is_cds_overlap", "overlapping_genes"]
        out.write('\t'.join(header) + '\n')

        for i, cand in enumerate(candidates):
            orf_id = f"ORF_{i+1}_{cand.gid}"
            tools = ",".join(sorted(set(t for t, s in cand.sources)))
            samples = ",".join(sorted(set(s for t, s in cand.sources)))
            blocks_str = ",".join(f"{s}-{e}" for s, e in cand.blocks)
            tool_scores_str = (",".join(f"{t}:{v:.3f}" if isinstance(v, float) else f"{t}:{v}"
                                        for t, v in sorted(cand.tool_scores.items()) if v is not None)
                               or "NA")
            tool_pvalues_str = (",".join(f"{t}:{p:.2e}" for t, p in sorted(cand.tool_pvalues.items()) if p is not None)
                                or "NA")
            if cand.subset_orfs:
                subset_orfs_str = ";".join(
                    f"{s['blocks']}|{s['tools']}|{s['samples']}" for s in cand.subset_orfs)
                num_subset_orfs = len(cand.subset_orfs)
            else:
                subset_orfs_str = "NA"
                num_subset_orfs = 0
            row = [orf_id, cand.chrom, cand.strand, str(cand.start), str(cand.end),
                   str(cand.length_aa), blocks_str, cand.gid, cand.tid, tools, samples,
                   tool_scores_str, tool_pvalues_str,
                   str(cand.total_reads), str(cand.unique_reads),
                   str(cand.total_psites), str(cand.unique_psites),
                   f"{cand.pN:.6f}", f"{cand.unique_pN:.6f}",
                   str(num_subset_orfs), subset_orfs_str,
                   cand.sequence, cand.start_codon, cand.aa_sequence,
                   "1" if cand.is_cds_overlap else "0",
                   ",".join(cand.overlapping_gene_ids) if cand.overlapping_gene_ids else "NA"]
            out.write('\t'.join(row) + '\n')

    with open(f"{prefix}.bed", 'w') as out:
        for i, cand in enumerate(candidates):
            orf_id = f"ORF_{i+1}_{cand.gid}"
            start0 = cand.start - 1
            end1 = cand.end
            block_sizes = []
            block_starts = []
            for s, e in cand.blocks:
                block_sizes.append(str(e - s + 1))
                block_starts.append(str((s - 1) - start0))
            out.write(f"{cand.chrom}\t{start0}\t{end1}\t{orf_id}\t0\t{cand.strand}\t"
                      f"{start0}\t{end1}\t0,0,0\t{len(cand.blocks)}\t"
                      f"{','.join(block_sizes)}\t{','.join(block_starts)}\n")

    with open(f"{prefix}.gtf", 'w') as out:
        out.write("##gff-version 2\n")
        source = "UnifiedRiboseq"
        for i, cand in enumerate(candidates):
            orf_id = f"ORF_{i+1}_{cand.gid}"
            tools = ",".join(sorted(set(t for t, s in cand.sources)))
            samples = ",".join(sorted(set(s for t, s in cand.sources)))
            num_tools = len(set(t for t, s in cand.sources))
            attr = (f'gene_id "{cand.gid}"; transcript_id "{cand.tid}"; orf_id "{orf_id}"; '
                    f'sources "{tools}"; samples "{samples}"; num_tools "{num_tools}";')
            for s, e in cand.blocks:
                out.write(f"{cand.chrom}\t{source}\texon\t{s}\t{e}\t.\t{cand.strand}\t.\t{attr}\n")
                out.write(f"{cand.chrom}\t{source}\tCDS\t{s}\t{e}\t.\t{cand.strand}\t.\t{attr}\n")


# Canonical tool names (as stored in ORFCandidate.sources) → file suffix
_TOOL_NAMES = [
    ("Ribo-TISH", "ribotish"),
    ("Ribotricer", "ribotricer"),
    ("RiboCode", "ribocode"),
    ("ORFquant", "orfquant"),
    ("PRICE", "price"),
]


def _write_per_tool_outputs(candidates: list, per_tool_prefix: str) -> None:
    """Write per-tool BED12/metadata/GTF files for each tool present in *candidates*."""
    for canonical, suffix in _TOOL_NAMES:
        tool_orfs = [c for c in candidates if canonical in {t for t, _s in c.sources}]
        if tool_orfs:
            _write_orf_outputs(tool_orfs, f"{per_tool_prefix}_{suffix}")
            print(f"  {canonical}: {len(tool_orfs)} ORFs → {per_tool_prefix}_{suffix}.*",
                  file=sys.stderr)
        else:
            print(f"  {canonical}: 0 ORFs (no output files written)", file=sys.stderr)


def main(argv=None):
    global _shared_gtf_index, _shared_fasta_path

    if _BIOPYTHON_IMPORT_ERROR is not None:
        print("Error: Biopython is required. Please install it with 'pip install biopython'", file=sys.stderr)
        sys.exit(1)
    if _PYFAIDX_IMPORT_ERROR is not None:
        print("Error: pyfaidx is required. Please install it with 'pip install pyfaidx'", file=sys.stderr)
        sys.exit(1)

    parser = argparse.ArgumentParser(description="Unify ORF predictions from multiple tools")
    parser.add_argument("--ribotish", nargs='+', help="Ribo-TISH output files")
    parser.add_argument("--ribotricer", nargs='+', help="Ribotricer output files")
    parser.add_argument("--ribocode", nargs='+', help="RiboCode GTF/TXT output files")
    parser.add_argument("--orfquant", nargs='+', help="ORFquant GTF output files")
    parser.add_argument("--price", nargs='+', help="PRICE ORF TSV output files (GEDI platform)")
    parser.add_argument("--gtf", required=True, help="Reference GTF file for coordinate mapping")
    parser.add_argument("--fasta", required=True, help="Genome FASTA file for validation")
    parser.add_argument("--output", required=True, help="Output prefix")
    parser.add_argument("--min_len", type=int, default=6, help="Minimum amino acid length")
    parser.add_argument("--bedgraph-dir", help="Directory containing RiboseQC bedgraph files (optional)")
    parser.add_argument("--sample-list", help="Comma-separated list of sample names for bedgraph stats (optional)")
    parser.add_argument("--threads", type=int, default=4, help="Number of threads for parallel processing (default: 4)")
    parser.add_argument("--stats-mode", choices=["auto", "preload", "stream"], default="auto",
                        help="Bedgraph statistics mode: preload for speed, stream for low memory, auto to choose based on input size (default: auto)")
    parser.add_argument("--skip-coverage-stats", action="store_true",
                        help="Skip coverage/coverage_uniq bedgraph aggregation and only compute P-site statistics")
    # Merging parameters
    parser.add_argument("--frame-merge-min-overlap", type=float, default=0.9,
                        help="Minimum overlap fraction of shorter ORF for single-exon frame-aware merging (default: 0.9)")
    parser.add_argument("--no-frame-merge", action="store_true",
                        help="Disable frame-aware merging (only use exact matches)")
    # Stage 0 filtering parameters (new)
    parser.add_argument("--exclude-tistypes", default="Annotated,annotated",
                        help="Comma-separated list of TisType/ORF_type values to exclude (default: Annotated,annotated)")
    parser.add_argument("--atg-only", action="store_true",
                        help="Only keep ORFs with ATG start codon (default: off)")
    # Stage 3 sequence clustering parameters (new)
    parser.add_argument("--seq-identity", type=float, default=0.9,
                        help="Minimum sequence identity for clustering (default: 0.9)")
    parser.add_argument("--min-length-ratio", type=float, default=0.8,
                        help="Minimum length ratio (shorter/longer) for clustering (default: 0.8)")
    parser.add_argument("--terminus-tolerance", type=int, default=3,
                        help="Start/stop position tolerance in bp for clustering (default: 3)")
    parser.add_argument("--seq-cluster", action="store_true",
                        help="Enable Stage 3 sequence-similarity clustering (disabled by default)")
    # Legacy option kept for back-compat (deprecated)
    parser.add_argument("--min-overlap", type=float, default=0.5,
                        help="[Deprecated] Minimum overlap fraction; use --seq-cluster to enable Stage 3")
    # Per-tool output mode: exact-match dedup per tool, no cross-tool merging
    parser.add_argument("--per-tool-output", type=str, default=None,
                        help="When set, output per-tool BED12/metadata/GTF files with this prefix "
                             "(e.g. 'ribotish_orfs', 'ribotricer_orfs', 'ribocode_orfs', 'orfquant_orfs') in addition "
                             "to the combined exact-dedup output at --output. "
                             "Skips Stage 2b (cross-tool frame-aware merge) and Stage 3 (seq clustering). "
                             "Use this for per-tool classification without cross-tool merging.")

    args = parser.parse_args(argv)
    
    # Build exclude_tistypes set
    exclude_tistypes = set()
    if args.exclude_tistypes:
        exclude_tistypes = {t.strip() for t in args.exclude_tistypes.split(',') if t.strip()}
    print(f"Stage 0 TisType exclusions: {exclude_tistypes or '(none)'}", file=sys.stderr)
    if args.atg_only:
        print("Stage 0: ATG-only mode enabled", file=sys.stderr)
    if args.skip_coverage_stats:
        print("Bedgraph stats: coverage aggregation disabled (--skip-coverage-stats)", file=sys.stderr)
    
    gtf_index = GTFIndex(args.gtf)
    
    print(f"Loading Genome FASTA: {args.fasta}...", file=sys.stderr)
    genome_fasta = Fasta(args.fasta)

    # ---------------------------------------------------------------------------
    # Determine parallelism.  Set module-level globals BEFORE any Pool.fork()
    # so child processes inherit them via Linux copy-on-write fork.
    # ---------------------------------------------------------------------------
    num_workers = min(args.threads, multiprocessing.cpu_count())
    _shared_gtf_index  = gtf_index
    _shared_fasta_path = args.fasta

    def _can_fork():
        try:
            multiprocessing.get_context('fork')
            return True
        except ValueError:
            return False

    use_parallel = num_workers > 1 and _can_fork()

    # Log which tools have input files
    print("\n=== Input Files Summary ===", file=sys.stderr)
    print(f"Ribo-TISH files: {len(args.ribotish) if args.ribotish else 0}", file=sys.stderr)
    if args.ribotish:
        for f in args.ribotish:
            print(f"  - {f}", file=sys.stderr)
    print(f"Ribotricer files: {len(args.ribotricer) if args.ribotricer else 0}", file=sys.stderr)
    if args.ribotricer:
        for f in args.ribotricer:
            print(f"  - {f}", file=sys.stderr)
    print(f"RiboCode files: {len(args.ribocode) if args.ribocode else 0}", file=sys.stderr)
    if args.ribocode:
        for f in args.ribocode:
            print(f"  - {f}", file=sys.stderr)
    print(f"ORFquant files: {len(args.orfquant) if args.orfquant else 0}", file=sys.stderr)
    if args.orfquant:
        for f in args.orfquant:
            print(f"  - {f}", file=sys.stderr)
    print(f"PRICE files: {len(args.price) if args.price else 0}", file=sys.stderr)
    if args.price:
        for f in args.price:
            print(f"  - {f}", file=sys.stderr)

    all_candidates = []

    # Track statistics by tool
    tool_stats   = defaultdict(lambda: {'count': 0, 'samples': set()})
    sample_stats = defaultdict(lambda: {'ribotish': 0, 'ribotricer': 0, 'ribocode': 0, 'orfquant': 0, 'price': 0})

    # Build task list: (tool_name, file_path, sample_id, min_len, excl, atg_only)
    parse_tasks = []
    if args.ribotish:
        for f in args.ribotish:
            sid = infer_sample_id_from_prediction_path(f, '_pred.txt')
            parse_tasks.append(('Ribo-TISH', f, sid, args.min_len, exclude_tistypes, args.atg_only))
    if args.ribotricer:
        for f in args.ribotricer:
            sid = infer_sample_id_from_prediction_path(f, '_translating_ORFs.tsv')
            parse_tasks.append(('Ribotricer', f, sid, args.min_len, exclude_tistypes, args.atg_only))
    if args.ribocode:
        for f in args.ribocode:
            # Detect suffix from filename for sample ID extraction
            if f.endswith('_collapsed.gtf.gz'):
                suffix = '_collapsed.gtf.gz'
            elif f.endswith('_collapsed.gtf'):
                suffix = '_collapsed.gtf'
            elif f.endswith('_collapsed.txt'):
                suffix = '_collapsed.txt'
            elif f.endswith('.gtf.gz'):
                suffix = '.gtf.gz'
            elif f.endswith('.gtf'):
                suffix = '.gtf'
            elif f.endswith('.txt'):
                suffix = '.txt'
            else:
                suffix = os.path.splitext(f)[1]  # Fallback
            sid = infer_sample_id_from_prediction_path(f, suffix)
            parse_tasks.append(('RiboCode', f, sid, args.min_len, exclude_tistypes, args.atg_only))
    if args.orfquant:
        for f in args.orfquant:
            sid = infer_sample_id_from_prediction_path(f, '_Detected_ORFs.gtf')
            parse_tasks.append(('ORFquant', f, sid, args.min_len, exclude_tistypes, args.atg_only))
    if args.price:
        for f in args.price:
            sid = infer_sample_id_from_prediction_path(f, '.orfs.tsv')
            parse_tasks.append(('PRICE', f, sid, args.min_len, exclude_tistypes, args.atg_only))

    # --- Parallel or sequential file parsing ---
    if use_parallel and len(parse_tasks) > 1:
        n_parse_workers = min(num_workers, len(parse_tasks))
        print(f"Parsing {len(parse_tasks)} input files using {n_parse_workers} workers (fork)...",
              file=sys.stderr)
        try:
            ctx = multiprocessing.get_context('fork')
            with ctx.Pool(n_parse_workers) as pool:
                for tool, sid, orfs in pool.imap_unordered(_parse_file_task, parse_tasks):
                    all_candidates.extend(orfs)
                    key = _TOOL_STATS_KEY[tool]
                    tool_stats[key]['count'] += len(orfs)
                    tool_stats[key]['samples'].add(sid)
                    sample_stats[sid][key] = len(orfs)
        except Exception as exc:
            print(f"Parallel parsing failed ({exc}), retrying sequentially...", file=sys.stderr)
            all_candidates.clear()
            for ts in [tool_stats, sample_stats]:
                ts.clear()
            use_parallel = False

    if not use_parallel or len(parse_tasks) <= 1:
        _parse_fn_map = {
            'Ribo-TISH': parse_ribotish,
            'Ribotricer': parse_ribotricer,
            'RiboCode':   parse_ribocode,
            'ORFquant':   parse_orfquant,
            'PRICE':      parse_price,
        }
        for tool, file_path, sid, min_len, excl, atg in parse_tasks:
            orfs = _parse_fn_map[tool](file_path, gtf_index, sid, min_len,
                                       exclude_tistypes=excl, atg_only=atg)
            all_candidates.extend(orfs)
            key = _TOOL_STATS_KEY[tool]
            tool_stats[key]['count'] += len(orfs)
            tool_stats[key]['samples'].add(sid)
            sample_stats[sid][key] = len(orfs)

    print(f"Total raw candidates: {len(all_candidates)}", file=sys.stderr)
    
    # Print statistics before merging
    print("\n=== Input Statistics (raw, per tool and per sample) ===", file=sys.stderr)
    print("By Tool:", file=sys.stderr)
    for tool in ['ribotish', 'ribotricer', 'ribocode', 'orfquant', 'price']:
        if tool in tool_stats:
            count = tool_stats[tool]['count']
            samples_ran = sorted(tool_stats[tool]['samples'])
            n = len(samples_ran)
            print(f"  {tool:12s}: {count:6d} ORFs from {n} sample(s): {', '.join(samples_ran)}", file=sys.stderr)
        else:
            print(f"  {tool:12s}: No input files provided", file=sys.stderr)
    print("By Sample (note: 0 means tool was not run on this sample):", file=sys.stderr)
    all_input_samples = sorted(set(sample_stats.keys()))
    for sample in all_input_samples:
        stats = sample_stats[sample]
        total = stats['ribotish'] + stats['ribotricer'] + stats['ribocode'] + stats['orfquant']
        if total > 0:
            print(f"  {sample:20s}: ribotish={stats['ribotish']:6d}, ribotricer={stats['ribotricer']:6d}, ribocode={stats['ribocode']:6d}, orfquant={stats['orfquant']:6d}, total={total:6d}", file=sys.stderr)
    
    # Stage 1: Exact match merging (same chrom, strand, and exact block coordinates).
    # Sort first so the same (gid, tid) is always chosen as the representative when
    # the same ORF appears in multiple samples (ensures deterministic output regardless
    # of parallel parse order).
    all_candidates.sort(key=lambda c: (c.chrom, c.strand, c.id_key, c.gid, c.tid))
    merged_candidates = {}
    for cand in all_candidates:
        if cand.id_key in merged_candidates:
            merged_candidates[cand.id_key].merge(cand)
        else:
            merged_candidates[cand.id_key] = cand
    
    print(f"After exact-match merging: {len(merged_candidates)}", file=sys.stderr)
    print(f"  Merged {len(all_candidates) - len(merged_candidates)} duplicates", file=sys.stderr)
    
    final_list = list(merged_candidates.values())
    del merged_candidates
    # Sort deterministically so ORF_IDs are consistent regardless of parse order
    final_list.sort(key=lambda c: (c.chrom, c.strand, c.start, c.end))

    del all_candidates
    gc.collect()

    skip_stage3 = not args.seq_cluster

    # ── Per-tool output mode ────────────────────────────────────────────────
    # When --per-tool-output is set, emit per-tool files based on the exact-match
    # deduplicated list (Stage 1 only) and skip cross-tool merging (Stage 2b) and
    # sequence clustering (Stage 3).  The combined exact-dedup set is also written
    # to args.output.* so downstream tools can consume the full set if needed.
    if args.per_tool_output:
        annotate_sequences_and_cds(final_list, use_parallel, num_workers, genome_fasta, gtf_index)
        _write_per_tool_outputs(final_list, args.per_tool_output)
        _write_orf_outputs(final_list, args.output)
        print(f"Per-tool outputs written to {args.per_tool_output}_{{ribotish,ribotricer,ribocode,orfquant}}.*",
              file=sys.stderr)
        print(f"Combined exact-dedup output written to {args.output}.*", file=sys.stderr)
        print(f"Done. Outputs written to {args.output}.*", file=sys.stderr)
        return

    # Stage 2: Frame-aware merging (single-exon only, overlap fraction threshold)
    if not args.no_frame_merge:
        frac = args.frame_merge_min_overlap
        print(f"Performing frame-aware merging (single-exon only, min_overlap_fraction={frac})...", file=sys.stderr)
        final_list = merge_frame_compatible_orfs(final_list, min_overlap_fraction=frac)
        print(f"After frame-aware merging: {len(final_list)}", file=sys.stderr)

    # Stage 3: Sequence-similarity clustering (disabled by default, opt-in with --seq-cluster)
    if not skip_stage3:
        annotate_sequences_and_cds(final_list, use_parallel, num_workers, genome_fasta, gtf_index)
        cds_overlap_count = sum(1 for c in final_list if c.is_cds_overlap)
        print(f"  CDS in-frame overlap: {cds_overlap_count} ORFs "
              f"({100*cds_overlap_count/max(len(final_list),1):.1f}%)", file=sys.stderr)
        print(f"Clustering by sequence similarity "
              f"(seq_identity={args.seq_identity}, min_length_ratio={args.min_length_ratio}, "
              f"terminus_tolerance={args.terminus_tolerance})...", file=sys.stderr)
        final_list = cluster_by_sequence(
            final_list,
            seq_identity=args.seq_identity,
            min_length_ratio=args.min_length_ratio,
            terminus_tolerance=args.terminus_tolerance,
        )
        print(f"After sequence clustering: {len(final_list)} representative ORFs", file=sys.stderr)
    
    # Calculate final statistics by tool
    print("\n=== Final Set Statistics (after cross-attribution via overlap grouping) ===", file=sys.stderr)
    print("Note: cross-attribution assigns subset ORF sources to representative, so all tools/samples may appear for each sample.", file=sys.stderr)
    final_tool_stats = defaultdict(int)
    final_sample_stats = defaultdict(lambda: defaultdict(int))
    
    for cand in final_list:
        tools_in_orf = set(t for t, s in cand.sources)
        for tool in tools_in_orf:
            final_tool_stats[tool] += 1
        
        samples_in_orf = set(s for t, s in cand.sources)
        for sample in samples_in_orf:
            for tool in tools_in_orf:
                final_sample_stats[sample][tool] += 1
    
    # Map canonical tool names used in ORFCandidate.sources to display names
    TOOL_KEYS = {
        'Ribo-TISH': 'ribotish',
        'Ribotricer': 'ribotricer',
        'RiboCode': 'ribocode',
        'ORFquant': 'orfquant',
    }

    print(f"Final unified ORFs: {len(final_list)}", file=sys.stderr)
    print("By Tool:", file=sys.stderr)
    for tool_key, tool_label in TOOL_KEYS.items():
        if tool_key in final_tool_stats:
            print(f"  {tool_label:12s}: {final_tool_stats[tool_key]:6d} ORFs", file=sys.stderr)

    print("By Sample:", file=sys.stderr)
    for sample in sorted(final_sample_stats.keys()):
        tool_counts = final_sample_stats[sample]
        ribotish_cnt   = tool_counts.get('Ribo-TISH', 0)
        ribotricer_cnt = tool_counts.get('Ribotricer', 0)
        ribocode_cnt   = tool_counts.get('RiboCode', 0)
        orfquant_cnt   = tool_counts.get('ORFquant', 0)
        total = ribotish_cnt + ribotricer_cnt + ribocode_cnt + orfquant_cnt
        if total > 0:
            print(f"  {sample:20s}: ribotish={ribotish_cnt:6d}, ribotricer={ribotricer_cnt:6d}, ribocode={ribocode_cnt:6d}, orfquant={orfquant_cnt:6d}, total={total:6d}", file=sys.stderr)
    
    # Calculate statistics from bedgraphs if provided (using optimized indexed version)
    per_sample_data = None
    if args.bedgraph_dir and args.sample_list:
        sample_list = args.sample_list.split(',')
        print(f"Calculating statistics from bedgraphs for {len(sample_list)} samples...", file=sys.stderr)

        metric_types = ['psite', 'psite_uniq']
        if not args.skip_coverage_stats:
            metric_types.extend(['coverage', 'coverage_uniq'])

        stats_mode = args.stats_mode
        if stats_mode == 'auto':
            total_bedgraph_bytes = _estimate_bedgraph_bytes(args.bedgraph_dir, sample_list, metric_types)
            if total_bedgraph_bytes <= 512 * 1024 * 1024 and len(final_list) <= 20000:
                stats_mode = 'preload'
            else:
                stats_mode = 'stream'
            print(
                f"Bedgraph stats mode auto-selected '{stats_mode}' "
                f"(files={total_bedgraph_bytes / (1024 * 1024):.1f} MiB, ORFs={len(final_list)})",
                file=sys.stderr
            )

        if stats_mode == 'preload':
            bedgraph_indices = load_bedgraph_indices(args.bedgraph_dir, sample_list, metric_types=metric_types)
            calculate_statistics_parallel(
                final_list,
                bedgraph_indices,
                sample_list,
                args.threads,
                metric_types=metric_types,
            )
            del bedgraph_indices
        else:
            per_sample_data = calculate_statistics_streaming(
                final_list,
                args.bedgraph_dir,
                sample_list,
                num_workers=args.threads,
                metric_types=metric_types,
            )
        gc.collect()

        # Write per-ORF expression outputs from per-sample data (streaming mode only)
        if per_sample_data:
            _write_expression_outputs(final_list, per_sample_data, args.output)
        else:
            print("Note: per-sample expression outputs not available (preload mode or no bedgraph data).",
                  file=sys.stderr)

    if skip_stage3:
        annotate_sequences_and_cds(final_list, use_parallel, num_workers, genome_fasta, gtf_index)
        cds_overlap_count = sum(1 for c in final_list if c.is_cds_overlap)
        print(f"  CDS in-frame overlap: {cds_overlap_count} ORFs "
              f"({100*cds_overlap_count/max(len(final_list),1):.1f}%)", file=sys.stderr)
    
    _write_orf_outputs(final_list, args.output)
    print(f"Done. Outputs written to {args.output}.*", file=sys.stderr)

if __name__ == "__main__":
    main(sys.argv[1:])
