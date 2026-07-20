#!/usr/bin/env python3
"""
Compute Yuanliang P-site purity metrics per ORF from RiboseQC bedgraph files.

Reads P-site and coverage bedgraph files, computes per-ORF per-sample:
  - p_site_GSE:        total P-site reads overlapping the ORF
  - reads_GSE:         total coverage reads overlapping the ORF
  - p_site_percentage: p_site_GSE / reads_GSE
  - p_site_postion_GSE: weighted average P-site position relative to ORF start
  - not_p_site_GSE:    reads_GSE - p_site_GSE

Outputs: psite_purity.tsv with per-sample and global aggregate columns.

Usage:
  python3 compute_psite_purity.py \
    --bed unified_orfs.bed \
    --riboseqc-dir results/riboseqc \
    --output psite_purity.tsv \
    --workers 8
"""

import gzip
import os
import sys
import argparse
import bisect
import time
from pathlib import Path
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np


# ---------------------------------------------------------------------------
# ORF index builder
# ---------------------------------------------------------------------------

def build_orf_index(bed_path):
    """
    Parse BED12 file and build per-chromosome sorted ORF index.

    Returns:
        orf_index: {chrom: {'starts': [], 'ends': [], 'ids': [], 'strands': [],
                             'blocks': [(s,e), ...]}}
        orf_meta:  [{orf_id, chrom, start, end, strand, blocks, orf_start, orf_end}]
    """
    orf_index = defaultdict(lambda: {
        'starts': [], 'ends': [], 'ids': [], 'strands': [], 'orf_starts': [], 'orf_ends': []
    })
    orf_meta = []

    open_func = gzip.open if str(bed_path).endswith('.gz') else open
    with open_func(bed_path, 'rt') as f:
        for line in f:
            if not line.strip():
                continue
            parts = line.strip().split('\t')
            if len(parts) < 12:
                continue

            chrom = parts[0]
            start0 = int(parts[1])     # 0-based
            end1 = int(parts[2])       # 1-based exclusive
            orf_id = parts[3]
            strand = parts[5]
            block_sizes = [int(x) for x in parts[10].split(',') if x]
            block_starts0 = [int(x) for x in parts[11].split(',') if x]

            # Convert to 1-based closed blocks [(s, e), ...]
            blocks = []
            for bs, bz in zip(block_starts0, block_sizes):
                s = start0 + bs + 1     # 1-based
                e = start0 + bs + bz     # 1-based inclusive
                blocks.append((s, e))

            orf_start = min(s for s, e in blocks)
            orf_end = max(e for s, e in blocks)

            idx = orf_index[chrom]
            idx['starts'].append(orf_start)
            idx['ends'].append(orf_end)
            idx['ids'].append(orf_id)
            idx['strands'].append(strand)
            idx['orf_starts'].append(orf_start)
            idx['orf_ends'].append(orf_end)

            orf_meta.append({
                'orf_id': orf_id,
                'chrom': chrom,
                'start': orf_start,
                'end': orf_end,
                'strand': strand,
                'blocks': blocks,
                'orf_start': orf_start,
                'orf_end': orf_end,
            })

    # Convert to numpy arrays for bisect
    for chrom in orf_index:
        idx = orf_index[chrom]
        order = np.argsort(idx['starts'])
        for key in ['starts', 'ends', 'ids', 'strands', 'orf_starts', 'orf_ends']:
            idx[key] = np.array(idx[key])[order]

    return dict(orf_index), orf_meta


# ---------------------------------------------------------------------------
# Bedgraph streamer — find overlapping ORFs via interval bisect
# ---------------------------------------------------------------------------

def find_overlapping_orfs(idx, bg_start, bg_end):
    """
    Find indices of ORFs overlapping a bedgraph interval [bg_start, bg_end] (1-based).

    Uses the fact that ORFs are sorted by start. Finds the range of ORFs whose
    start <= bg_end, then checks which of those have end >= bg_start.
    """
    starts = idx['starts']
    ends = idx['ends']

    # All ORFs with start <= bg_end
    right = int(bisect.bisect_right(starts, bg_end))
    if right == 0:
        return []

    # Of those, keep ones with end >= bg_start
    matched = []
    for i in range(right):
        if ends[i] >= bg_start:
            matched.append(i)
    return matched


def compute_orf_position(orf_entry, bg_start, bg_end):
    """
    Compute the relative position of a P-site interval within the ORF.

    For each ORF block, compute the overlap and weighted position.
    Returns list of (overlap_length, position_from_orf_start) tuples.
    """
    positions = []
    orf_start = orf_entry['orf_start']
    strand = orf_entry['strand']

    for bs, be in orf_entry['blocks']:
        ov_start = max(bs, bg_start)
        ov_end = min(be, bg_end)
        if ov_start <= ov_end:
            if strand == '+':
                pos_from_start = ov_start - orf_start + 1
            else:
                pos_from_start = orf_entry['orf_end'] - ov_end + 1
            positions.append((ov_end - ov_start + 1, pos_from_start))

    return positions


# ---------------------------------------------------------------------------
# Per-worker function: process a batch of samples
# ---------------------------------------------------------------------------

def process_samples(worker_args):
    """
    Process a batch of samples against all ORFs.

    Args:
        worker_args: (worker_id, sample_list, orf_index_serializable, bedgraph_dir)

    Returns:
        dict: {orf_id: {sample}_p_site_GSE, {sample}_reads_GSE, ...}
    """
    worker_id, sample_list, orf_index_dict, bedgraph_dir = worker_args

    # Reconstruct orf_index with numpy arrays
    orf_index = {}
    for chrom, idx in orf_index_dict.items():
        orf_index[chrom] = {
            'starts': np.array(idx['starts']),
            'ends': np.array(idx['ends']),
            'ids': np.array(idx['ids']),
            'strands': np.array(idx['strands']),
            'orf_starts': np.array(idx['orf_starts']),
            'orf_ends': np.array(idx['orf_ends']),
        }

    results = defaultdict(dict)  # orf_id -> {metric: value}

    for si, sample in enumerate(sample_list):
        for strand_label, strand_char in [('plus', '+'), ('minus', '-')]:
            psite_file = os.path.join(
                bedgraph_dir, f"{sample}_P_sites_{strand_label}.bedgraph")
            cov_file = os.path.join(
                bedgraph_dir, f"{sample}_coverage_{strand_label}.bedgraph")

            # ── P-site bedgraph ──
            if os.path.exists(psite_file):
                with open(psite_file, 'r') as f:
                    for line in f:
                        if line.startswith('#') or not line.strip():
                            continue
                        parts = line.strip().split('\t')
                        if len(parts) < 3:
                            continue
                        chrom, start0, end1 = parts[0], parts[1], parts[2]
                        value = float(parts[3]) if len(parts) >= 4 else 1.0
                        try:
                            bg_start = int(start0) + 1  # 1-based
                            bg_end = int(end1)
                        except ValueError:
                            continue

                        if chrom not in orf_index:
                            continue

                        idx = orf_index[chrom]
                        matched = find_overlapping_orfs(idx, bg_start, bg_end)
                        for midx in matched:
                            orf_id = idx['ids'][midx]
                            cur = results[orf_id]

                            # P-site reads
                            key_psite = f"{sample}_p_site_GSE"
                            cur[key_psite] = cur.get(key_psite, 0.0) + value

                            # P-site position (weighted by value)
                            # Compute position relative to ORF start
                            ov_start = max(idx['orf_starts'][midx], bg_start)
                            ov_end = min(idx['orf_ends'][midx], bg_end)
                            strand = idx['strands'][midx]
                            if strand == '+':
                                pos = ov_start - idx['orf_starts'][midx] + 1
                            else:
                                pos = idx['orf_ends'][midx] - ov_end + 1

                            key_pos_sum = f"{sample}_p_site_pos_sum"
                            key_pos_wt = f"{sample}_p_site_pos_wt"
                            cur[key_pos_sum] = cur.get(key_pos_sum, 0.0) + value * pos
                            cur[key_pos_wt] = cur.get(key_pos_wt, 0.0) + value

            # ── Coverage bedgraph ──
            if os.path.exists(cov_file):
                with open(cov_file, 'r') as f:
                    for line in f:
                        if line.startswith('#') or not line.strip():
                            continue
                        parts = line.strip().split('\t')
                        if len(parts) < 3:
                            continue
                        chrom, start0, end1 = parts[0], parts[1], parts[2]
                        value = float(parts[3]) if len(parts) >= 4 else 1.0
                        try:
                            bg_start = int(start0) + 1
                            bg_end = int(end1)
                        except ValueError:
                            continue

                        if chrom not in orf_index:
                            continue

                        idx = orf_index[chrom]
                        matched = find_overlapping_orfs(idx, bg_start, bg_end)
                        for midx in matched:
                            orf_id = idx['ids'][midx]
                            cur = results[orf_id]
                            key_reads = f"{sample}_reads_GSE"
                            cur[key_reads] = cur.get(key_reads, 0.0) + value

        if (si + 1) % 5 == 0 or (si + 1) == len(sample_list):
            print(f"  [worker {worker_id}] {si + 1}/{len(sample_list)} samples done",
                  flush=True)

    # Convert defaultdict to regular dict for pickling
    return dict(results)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Compute Yuanliang P-site purity metrics from RiboseQC bedgraphs')
    parser.add_argument('--bed', required=True,
                        help='Path to unified_orfs.bed (BED12)')
    parser.add_argument('--riboseqc-dir', required=True,
                        help='Directory containing RiboseQC bedgraph files')
    parser.add_argument('--output', default='psite_purity.tsv',
                        help='Output TSV path')
    parser.add_argument('--workers', type=int, default=8,
                        help='Number of parallel workers')
    parser.add_argument('--limit', type=int, default=0,
                        help='Limit ORFs for testing (0 = all)')
    args = parser.parse_args()

    t0 = time.time()

    # 1. Build ORF index
    print(f"Loading ORFs from {args.bed} ...", flush=True)
    orf_index, orf_meta = build_orf_index(args.bed)
    if args.limit > 0:
        orf_meta = orf_meta[:args.limit]
        # Filter orf_index to only keep limited ORFs
        limited_ids = {m['orf_id'] for m in orf_meta}
        for chrom in list(orf_index.keys()):
            idx = orf_index[chrom]
            keep = np.array([i for i, oid in enumerate(idx['ids']) if oid in limited_ids])
            if len(keep) == 0:
                del orf_index[chrom]
            else:
                for key in idx:
                    idx[key] = idx[key][keep]

    orf_lookup = {m['orf_id']: m for m in orf_meta}
    print(f"  {len(orf_meta)} ORFs indexed across {len(orf_index)} chromosomes",
          flush=True)

    # 2. Discover samples from bedgraph directory
    samples = set()
    for f in Path(args.riboseqc_dir).iterdir():
        if f.is_file() and '_P_sites_plus.bedgraph' in f.name:
            sample = f.name.replace('_P_sites_plus.bedgraph', '')
            samples.add(sample)
    sample_list = sorted(samples)
    print(f"Found {len(sample_list)} samples", flush=True)

    if not sample_list:
        print("ERROR: No P-site bedgraph files found", file=sys.stderr)
        sys.exit(1)

    # 3. Prepare serializable ORF index for workers
    orf_index_serializable = {}
    for chrom, idx in orf_index.items():
        orf_index_serializable[chrom] = {
            'starts': idx['starts'].tolist(),
            'ends': idx['ends'].tolist(),
            'ids': idx['ids'].tolist(),
            'strands': idx['strands'].tolist(),
            'orf_starts': idx['orf_starts'].tolist(),
            'orf_ends': idx['orf_ends'].tolist(),
        }

    # 4. Distribute samples across workers
    n_workers = min(args.workers, len(sample_list))
    worker_samples = [[] for _ in range(n_workers)]
    for i, s in enumerate(sample_list):
        worker_samples[i % n_workers].append(s)

    print(f"Processing with {n_workers} workers ...", flush=True)

    # 5. Parallel processing
    all_results = {}
    with ProcessPoolExecutor(max_workers=n_workers) as executor:
        futures = []
        for wid in range(n_workers):
            wargs = (wid, worker_samples[wid], orf_index_serializable, args.riboseqc_dir)
            futures.append(executor.submit(process_samples, wargs))

        for fut in as_completed(futures):
            worker_results = fut.result()
            for orf_id, metrics in worker_results.items():
                if orf_id not in all_results:
                    all_results[orf_id] = {}
                all_results[orf_id].update(metrics)

    print(f"  Collected metrics for {len(all_results)} ORFs", flush=True)

    # 6. Write output
    print(f"Writing {args.output} ...", flush=True)

    # Compute derived metrics and build header
    # Per-sample columns: {sample}_p_site_GSE, {sample}_reads_GSE,
    #   {sample}_p_site_pct, {sample}_p_site_pos
    header = ['orf_id', 'chrom', 'start', 'end', 'strand']
    for s in sample_list:
        header.append(f"{s}_p_site_GSE")
        header.append(f"{s}_reads_GSE")
        header.append(f"{s}_p_site_pct")
        header.append(f"{s}_p_site_pos")
        header.append(f"{s}_not_p_site_GSE")
    header.extend([
        'total_psites', 'total_reads',
        'global_p_site_pct', 'global_p_site_pos', 'global_p_site_pos_sd',
        'n_samples_with_psites'
    ])

    with open(args.output, 'w') as out:
        out.write('\t'.join(header) + '\n')

        n_written = 0
        for meta in orf_meta:
            orf_id = meta['orf_id']
            row = [orf_id, meta['chrom'], str(meta['start']),
                   str(meta['end']), meta['strand']]

            metrics = all_results.get(orf_id, {})

            # Per-sample values
            total_psites = 0.0
            total_reads = 0.0
            global_pos_sum = 0.0
            global_pos_wt = 0.0
            all_positions = []
            n_samples_psites = 0

            for s in sample_list:
                psite = metrics.get(f"{s}_p_site_GSE", 0.0)
                reads = metrics.get(f"{s}_reads_GSE", 0.0)
                pct = round(psite / reads, 4) if reads > 0 else 0.0
                pos_wt = metrics.get(f"{s}_p_site_pos_wt", 0.0)
                pos_avg = round(pos_wt / psite, 2) if psite > 0 else 0.0
                not_psite = max(0.0, reads - psite)

                row.extend([
                    str(int(psite)), str(int(reads)),
                    f"{pct:.4f}", f"{pos_avg:.2f}", str(int(not_psite))
                ])

                total_psites += psite
                total_reads += reads
                if psite > 0:
                    global_pos_sum += metrics.get(f"{s}_p_site_pos_sum", 0.0)
                    global_pos_wt += pos_wt
                    n_samples_psites += 1

            # Global aggregates
            global_pct = round(total_psites / total_reads, 4) if total_reads > 0 else 0.0
            global_pos = round(global_pos_sum / max(global_pos_wt, 1.0), 2)
            # Stddev of per-sample p_site_pos values
            sample_positions = []
            for s in sample_list:
                pos_wt = metrics.get(f"{s}_p_site_pos_wt", 0.0)
                psite = metrics.get(f"{s}_p_site_GSE", 0.0)
                if psite > 0:
                    sample_positions.append(round(pos_wt / psite, 2))
            if len(sample_positions) >= 2:
                pos_sd = round(float(np.std(sample_positions)), 2)
            else:
                pos_sd = 0.0

            row.extend([
                str(int(total_psites)), str(int(total_reads)),
                f"{global_pct:.4f}", f"{global_pos:.2f}", f"{pos_sd:.2f}",
                str(n_samples_psites)
            ])

            out.write('\t'.join(row) + '\n')
            n_written += 1

    elapsed = time.time() - t0
    print(f"Done. {n_written} ORFs written in {elapsed:.1f}s", flush=True)


if __name__ == '__main__':
    main()
