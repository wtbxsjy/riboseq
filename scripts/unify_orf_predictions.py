#!/usr/bin/env python3
import argparse
import sys
import os
import csv
import re
from typing import List, Dict, Tuple, Set, Optional
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
import multiprocessing
import threading

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
    'ORFquant':   'orfquant',
}

# Try to import biopython and pyfaidx
try:
    from Bio import SeqIO
    from Bio.Seq import Seq
except ImportError:
    print("Error: Biopython is required. Please install it with 'pip install biopython'", file=sys.stderr)
    sys.exit(1)

try:
    from pyfaidx import Fasta
except ImportError:
    print("Error: pyfaidx is required. Please install it with 'pip install pyfaidx'", file=sys.stderr)
    sys.exit(1)

# GTF Parsing Helper
class GTFIndex:
    def __init__(self, gtf_file):
        self.transcripts = {} # tid -> {chrom, strand, exons: [(s,e), ...], cds: [(s,e), ...]}
        self.gene_map = {} # tid -> gid
        self.gene_names = {} # gid -> gene_name
        self.cds_by_chrom = {}   # chrom -> sorted list of (start, end, strand, gid)
        self.gene_by_chrom = {}  # chrom -> sorted list of (start, end, strand, gid)
        self._load_gtf(gtf_file)

    def _load_gtf(self, gtf_file):
        print(f"Loading GTF: {gtf_file}...", file=sys.stderr)
        cds_temp = defaultdict(list)
        gene_temp = defaultdict(list)
        with open(gtf_file, 'r') as f:
            for line in f:
                if line.startswith('#'): continue
                parts = line.strip().split('\t')
                if len(parts) < 9: continue

                feature = parts[2]
                attributes = self._parse_attributes(parts[8])
                chrom = parts[0]
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
    
    # For multi-exon ORFs, calculate exon-level overlap
    overlap_bp = 0
    for s1, e1 in cand1.blocks:
        for s2, e2 in cand2.blocks:
            os = max(s1, s2)
            oe = min(e1, e2)
            if os < oe:
                overlap_bp += oe - os
    
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
                continue
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
        # Sort by start position for a narrow comparison window
        indices_sorted = sorted(indices, key=lambda i: candidates[i].start)
        for ii in range(n):
            idx1 = indices_sorted[ii]
            c1 = candidates[idx1]
            for jj in range(ii + 1, n):
                idx2 = indices_sorted[jj]
                c2 = candidates[idx2]

                # Early-exit: if start positions are too far apart to share terminus
                # (and the ORFs are not likely to share a stop), skip
                # (max ORF is ~3000 aa → 9000 nt; we just avoid extremely distant pairs)
                if abs(c2.start - c1.start) > 100000:
                    break

                # Condition 2: shared start OR stop within tolerance
                shared_start = abs(c1.start - c2.start) <= terminus_tolerance
                shared_stop  = abs(c1.end   - c2.end)   <= terminus_tolerance
                if not (shared_start or shared_stop):
                    continue

                # Condition 3a: length ratio
                len1, len2 = c1.length_aa, c2.length_aa
                if len1 == 0 or len2 == 0:
                    continue
                shorter_len = min(len1, len2)
                longer_len  = max(len1, len2)
                if shorter_len / longer_len < min_length_ratio:
                    continue

                # Condition 3b: sequence identity (only if sequences available)
                aa1 = c1.aa_sequence
                aa2 = c2.aa_sequence
                if aa1 and aa2:
                    identity = _seq_identity_from_shared_terminus(aa1, aa2)
                    if identity < seq_identity:
                        continue
                # If sequences not available, fall back to length-ratio-only

                uf.union(idx1, idx2)

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
        with open(file_path, 'r') as f:
            first_line = f.readline().strip()
            if first_line.startswith('#') and ('placeholder' in first_line.lower() or 'no orfs' in first_line.lower() or 'insufficient' in first_line.lower()):
                print(f"  -> Placeholder file detected, skipping", file=sys.stderr)
                return []
    except Exception:
        pass
    
    filtered_tis = 0
    filtered_atg = 0
    
    try:
        with open(file_path, 'r') as f:
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
                        chrom = match.group(1)
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
        with open(file_path, 'r') as f:
            first_line = f.readline().strip()
            if first_line.startswith('#') and ('placeholder' in first_line.lower() or 'no orfs' in first_line.lower() or 'insufficient' in first_line.lower()):
                print(f"  -> Placeholder file detected, skipping", file=sys.stderr)
                return []
    except Exception:
        pass
    
    filtered_tis = 0
    filtered_atg = 0
    
    try:
        with open(file_path, 'r') as f:
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
        with open(file_path, 'r') as f:
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
        
        with open(file_path, 'r') as f:
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
                    current_attrs['chrom'] = parts[0]
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
            with open(bedgraph_file, 'r') as f:
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


def load_bedgraph_indices(bedgraph_dir, sample_list):
    """
    Pre-load all bedgraph files into indexed structures.
    Returns: dict of sample -> strand -> type -> BedgraphIndex
    """
    if not bedgraph_dir or not os.path.exists(bedgraph_dir):
        return None
    
    print(f"Pre-loading bedgraph files for {len(sample_list)} samples...", file=sys.stderr)
    indices = {}
    
    for sample in sample_list:
        indices[sample] = {}
        for strand_suffix in ['plus', 'minus']:
            indices[sample][strand_suffix] = {}
            
            # P-site bedgraphs
            psite_file = os.path.join(bedgraph_dir, f"{sample}_P_sites_{strand_suffix}.bedgraph")
            psite_uniq_file = os.path.join(bedgraph_dir, f"{sample}_P_sites_uniq_{strand_suffix}.bedgraph")
            coverage_file = os.path.join(bedgraph_dir, f"{sample}_coverage_{strand_suffix}.bedgraph")
            coverage_uniq_file = os.path.join(bedgraph_dir, f"{sample}_coverage_uniq_{strand_suffix}.bedgraph")
            
            indices[sample][strand_suffix]['psite'] = BedgraphIndex(psite_file)
            indices[sample][strand_suffix]['psite_uniq'] = BedgraphIndex(psite_uniq_file)
            indices[sample][strand_suffix]['coverage'] = BedgraphIndex(coverage_file)
            indices[sample][strand_suffix]['coverage_uniq'] = BedgraphIndex(coverage_uniq_file)
    
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
        with open(bedgraph_file, 'r') as f:
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


def calculate_statistics_from_indices(cand, bedgraph_indices, sample_list):
    """
    Calculate statistics using pre-loaded bedgraph indices (fast version)
    """
    if not bedgraph_indices:
        return
    
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
            total_psites += idx['psite'].count_in_region(cand.chrom, block_start, block_end)
            unique_psites += idx['psite_uniq'].count_in_region(cand.chrom, block_start, block_end)
            total_reads += idx['coverage'].count_in_region(cand.chrom, block_start, block_end)
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
    candidates, bedgraph_indices, sample_list = batch_data
    
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
                total_psites += idx['psite'].count_in_region(chrom, block_start, block_end)
                unique_psites += idx['psite_uniq'].count_in_region(chrom, block_start, block_end)
                total_reads += idx['coverage'].count_in_region(chrom, block_start, block_end)
                unique_reads += idx['coverage_uniq'].count_in_region(chrom, block_start, block_end)
        
        results.append({
            'id_key': cand_dict['id_key'],
            'total_psites': total_psites,
            'unique_psites': unique_psites,
            'total_reads': total_reads,
            'unique_reads': unique_reads
        })
    
    return results


def calculate_statistics_parallel(final_list, bedgraph_indices, sample_list, num_workers=None):
    """
    Calculate P-site statistics for all candidates in parallel.
    """
    if not bedgraph_indices or not sample_list:
        return
    
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
            calculate_statistics_from_indices(cand, bedgraph_indices, sample_list)
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
            futures.append(executor.submit(process_candidate_batch, (batch, bedgraph_indices, sample_list)))
        
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
        'ORFquant':   parse_orfquant,
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


def main():
    global _shared_gtf_index, _shared_fasta_path

    parser = argparse.ArgumentParser(description="Unify ORF predictions from multiple tools")
    parser.add_argument("--ribotish", nargs='+', help="Ribo-TISH output files")
    parser.add_argument("--ribotricer", nargs='+', help="Ribotricer output files")
    parser.add_argument("--orfquant", nargs='+', help="ORFquant GTF output files")
    parser.add_argument("--gtf", required=True, help="Reference GTF file for coordinate mapping")
    parser.add_argument("--fasta", required=True, help="Genome FASTA file for validation")
    parser.add_argument("--output", required=True, help="Output prefix")
    parser.add_argument("--min_len", type=int, default=6, help="Minimum amino acid length")
    parser.add_argument("--bedgraph-dir", help="Directory containing RiboseQC bedgraph files (optional)")
    parser.add_argument("--sample-list", help="Comma-separated list of sample names for bedgraph stats (optional)")
    parser.add_argument("--threads", type=int, default=4, help="Number of threads for parallel processing (default: 4)")
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
    
    args = parser.parse_args()
    
    # Build exclude_tistypes set
    exclude_tistypes = set()
    if args.exclude_tistypes:
        exclude_tistypes = {t.strip() for t in args.exclude_tistypes.split(',') if t.strip()}
    print(f"Stage 0 TisType exclusions: {exclude_tistypes or '(none)'}", file=sys.stderr)
    if args.atg_only:
        print("Stage 0: ATG-only mode enabled", file=sys.stderr)
    
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
    print(f"ORFquant files: {len(args.orfquant) if args.orfquant else 0}", file=sys.stderr)
    if args.orfquant:
        for f in args.orfquant:
            print(f"  - {f}", file=sys.stderr)

    all_candidates = []

    # Track statistics by tool
    tool_stats   = defaultdict(lambda: {'count': 0, 'samples': set()})
    sample_stats = defaultdict(lambda: {'ribotish': 0, 'ribotricer': 0, 'orfquant': 0})

    # Build task list: (tool_name, file_path, sample_id, min_len, excl, atg_only)
    parse_tasks = []
    if args.ribotish:
        for f in args.ribotish:
            sid = os.path.basename(f).split('.')[0].replace('_pred', '')
            parse_tasks.append(('Ribo-TISH', f, sid, args.min_len, exclude_tistypes, args.atg_only))
    if args.ribotricer:
        for f in args.ribotricer:
            sid = os.path.basename(f).split('.')[0].replace('_translating_ORFs', '')
            parse_tasks.append(('Ribotricer', f, sid, args.min_len, exclude_tistypes, args.atg_only))
    if args.orfquant:
        for f in args.orfquant:
            sid = os.path.basename(f).split('.')[0].replace('_Detected_ORFs', '')
            parse_tasks.append(('ORFquant', f, sid, args.min_len, exclude_tistypes, args.atg_only))

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
            'ORFquant':   parse_orfquant,
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
    for tool in ['ribotish', 'ribotricer', 'orfquant']:
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
        total = stats['ribotish'] + stats['ribotricer'] + stats['orfquant']
        if total > 0:
            print(f"  {sample:20s}: ribotish={stats['ribotish']:6d}, ribotricer={stats['ribotricer']:6d}, orfquant={stats['orfquant']:6d}, total={total:6d}", file=sys.stderr)
    
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
    # Sort deterministically so ORF_IDs are consistent regardless of parse order
    final_list.sort(key=lambda c: (c.chrom, c.strand, c.start, c.end))

    skip_stage3 = not args.seq_cluster

    # Stage 2 preparation: extract sequences + annotate CDS overlap for all candidates.
    # (Already sorted deterministically above; no additional sort needed for I/O locality.)

    if use_parallel and len(final_list) > 200:
        # Parallel: combine sequence extraction + CDS annotation in one worker pass.
        # Workers inherit _shared_gtf_index via fork; each opens its own Fasta handle.
        chunk_size = max(200, len(final_list) // (num_workers * 4))
        chunk_data = [
            (c.id_key, c.chrom, c.strand, c.blocks, c.frame, c.length_nt)
            for c in final_list
        ]
        chunks     = [chunk_data[i:i + chunk_size]
                      for i in range(0, len(chunk_data), chunk_size)]
        cand_lookup = {c.id_key: c for c in final_list}

        print(f"Extracting sequences + annotating CDS overlap "
              f"({len(final_list)} candidates, {num_workers} workers, "
              f"{len(chunks)} chunks)...", file=sys.stderr)
        try:
            ctx = multiprocessing.get_context('fork')
            with ctx.Pool(num_workers) as pool:
                all_chunk_results = pool.map(_extract_annotate_chunk, chunks)
            for chunk_results in all_chunk_results:
                for id_key, nt_seq, aa_seq, start_codon, is_cds, gene_ids in chunk_results:
                    cand = cand_lookup.get(id_key)
                    if cand:
                        cand.sequence             = nt_seq
                        cand.aa_sequence          = aa_seq
                        cand.start_codon          = start_codon
                        cand.is_cds_overlap       = is_cds
                        cand.overlapping_gene_ids = gene_ids
        except Exception as exc:
            print(f"Parallel extract/annotate failed ({exc}), retrying sequentially...",
                  file=sys.stderr)
            for cand in final_list:
                extract_sequence(cand, genome_fasta)
            annotate_cds_overlap(final_list, gtf_index)
    else:
        # Sequential path (threads=1, small dataset, or non-fork platform)
        print(f"Extracting sequences ({len(final_list)} candidates)...", file=sys.stderr)
        for cand in final_list:
            extract_sequence(cand, genome_fasta)
        print(f"Annotating CDS overlap and overlapping genes...", file=sys.stderr)
        annotate_cds_overlap(final_list, gtf_index)

    cds_overlap_count = sum(1 for c in final_list if c.is_cds_overlap)
    print(f"  CDS in-frame overlap: {cds_overlap_count} ORFs "
          f"({100*cds_overlap_count/max(len(final_list),1):.1f}%)", file=sys.stderr)

    # Stage 2: Frame-aware merging (single-exon only, overlap fraction threshold)
    if not args.no_frame_merge:
        frac = args.frame_merge_min_overlap
        print(f"Performing frame-aware merging (single-exon only, min_overlap_fraction={frac})...", file=sys.stderr)
        final_list = merge_frame_compatible_orfs(final_list, min_overlap_fraction=frac)
        print(f"After frame-aware merging: {len(final_list)}", file=sys.stderr)

    # Stage 3: Sequence-similarity clustering (disabled by default, opt-in with --seq-cluster)
    if not skip_stage3:
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
    TOOL_KEYS = {'Ribo-TISH': 'ribotish', 'Ribotricer': 'ribotricer', 'ORFquant': 'orfquant'}

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
        orfquant_cnt   = tool_counts.get('ORFquant', 0)
        total = ribotish_cnt + ribotricer_cnt + orfquant_cnt
        if total > 0:
            print(f"  {sample:20s}: ribotish={ribotish_cnt:6d}, ribotricer={ribotricer_cnt:6d}, orfquant={orfquant_cnt:6d}, total={total:6d}", file=sys.stderr)
    
    # Calculate statistics from bedgraphs if provided (using optimized indexed version)
    if args.bedgraph_dir and args.sample_list:
        sample_list = args.sample_list.split(',')
        print(f"Calculating statistics from bedgraphs for {len(sample_list)} samples...", file=sys.stderr)
        
        # Load bedgraph files into indexed structures (one-time cost)
        bedgraph_indices = load_bedgraph_indices(args.bedgraph_dir, sample_list)
        
        # Use parallel processing for large datasets
        calculate_statistics_parallel(final_list, bedgraph_indices, sample_list, args.threads)
    
    # Write metadata.tsv with extended columns
    with open(f"{args.output}.metadata.tsv", 'w') as out:
        header = ["orf_id", "chrom", "strand", "start", "end", "length_aa", "exon_blocks", 
                  "gene_id", "transcript_id", "tools", "samples", 
                  "tool_scores", "tool_pvalues", 
                  "total_reads", "unique_reads", "total_psites", "unique_psites", "pN", "unique_pN",
                  "num_subset_orfs", "subset_orfs",
                  "sequence", "start_codon", "aa_sequence", "is_cds_overlap", "overlapping_genes"]
        out.write('\t'.join(header) + '\n')
        
        for i, cand in enumerate(final_list):
            orf_id = f"ORF_{i+1}_{cand.gid}"
            tools = ",".join(sorted(list(set(t for t, s in cand.sources))))
            samples = ",".join(sorted(list(set(s for t, s in cand.sources))))
            blocks_str = ",".join(f"{s}-{e}" for s, e in cand.blocks)
            
            # Format tool_scores as tool1:score1,tool2:score2
            tool_scores_str = ",".join(f"{t}:{s:.3f}" if isinstance(s, float) else f"{t}:{s}" 
                                      for t, s in sorted(cand.tool_scores.items()) if s is not None) or "NA"
            
            # Format tool_pvalues as tool1:pval1,tool2:pval2
            tool_pvalues_str = ",".join(f"{t}:{p:.2e}" for t, p in sorted(cand.tool_pvalues.items()) if p is not None) or "NA"
            
            # Format subset_orfs as blocks|tools|samples;blocks|tools|samples;...
            if cand.subset_orfs:
                subset_strs = []
                for subset in cand.subset_orfs:
                    subset_str = f"{subset['blocks']}|{subset['tools']}|{subset['samples']}"
                    subset_strs.append(subset_str)
                subset_orfs_str = ";".join(subset_strs)
                num_subset_orfs = len(cand.subset_orfs)
            else:
                subset_orfs_str = "NA"
                num_subset_orfs = 0
            
            row = [orf_id, cand.chrom, cand.strand, str(cand.start), str(cand.end), str(cand.length_aa), 
                   blocks_str, cand.gid, cand.tid, tools, samples, 
                   tool_scores_str, tool_pvalues_str,
                   str(cand.total_reads), str(cand.unique_reads), 
                   str(cand.total_psites), str(cand.unique_psites),
                   f"{cand.pN:.6f}", f"{cand.unique_pN:.6f}",
                   str(num_subset_orfs), subset_orfs_str,
                   cand.sequence,
                   cand.start_codon,
                   cand.aa_sequence,
                   "1" if cand.is_cds_overlap else "0",
                   ",".join(cand.overlapping_gene_ids) if cand.overlapping_gene_ids else "NA"]
            out.write('\t'.join(row) + '\n')
            
    with open(f"{args.output}.bed", 'w') as out:
        for i, cand in enumerate(final_list):
            orf_id = f"ORF_{i+1}_{cand.gid}"
            
            chrom = cand.chrom
            start0 = cand.start - 1
            end1 = cand.end
            name = orf_id
            score = "0"
            strand = cand.strand
            thickStart = start0
            thickEnd = end1
            rgb = "0,0,0"
            blockCount = len(cand.blocks)
            
            blockSizes = []
            blockStarts = []
            
            for s, e in cand.blocks:
                size = e - s + 1
                rel_start = (s - 1) - start0
                blockSizes.append(str(size))
                blockStarts.append(str(rel_start))
                
            out.write(f"{chrom}\t{start0}\t{end1}\t{name}\t{score}\t{strand}\t{thickStart}\t{thickEnd}\t{rgb}\t{blockCount}\t{','.join(blockSizes)}\t{','.join(blockStarts)}\n")
    
    with open(f"{args.output}.gtf", 'w') as out:
        out.write("##gff-version 2\n")
        source = "UnifiedRiboseq"
        for i, cand in enumerate(final_list):
            orf_id = f"ORF_{i+1}_{cand.gid}"
            gene_id = cand.gid
            tid = cand.tid
            tools = ",".join(sorted(list(set(t for t, s in cand.sources))))
            samples = ",".join(sorted(list(set(s for t, s in cand.sources))))
            num_tools = len(set(t for t, s in cand.sources))
            
            attr_base = f'gene_id "{gene_id}"; transcript_id "{tid}"; orf_id "{orf_id}"; sources "{tools}"; samples "{samples}"; num_tools "{num_tools}";'
            
            for s, e in cand.blocks:
                out.write(f"{cand.chrom}\t{source}\texon\t{s}\t{e}\t.\t{cand.strand}\t.\t{attr_base}\n")
                out.write(f"{cand.chrom}\t{source}\tCDS\t{s}\t{e}\t.\t{cand.strand}\t.\t{attr_base}\n")

    print(f"Done. Outputs written to {args.output}.*", file=sys.stderr)

if __name__ == "__main__":
    main()
