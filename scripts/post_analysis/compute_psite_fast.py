#!/usr/bin/env python3
"""
Fast P-site purity computation using bedtools map.

Uses bedtools map to sum bedgraph values over ORF intervals — orders of
magnitude faster than Python line-by-line scanning (seconds vs hours).

Usage:
  python3 compute_psite_fast.py \
    --bed prelim_orfs_for_psite.bed \
    --riboseqc-dir results/riboseqc \
    --output psite_purity.tsv \
    --bedtools /usr/bin/bedtools \
    --workers 8
"""

import argparse
import os
import subprocess
import sys
import tempfile
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path


def run_bedtools_map(bed_file, bedgraph_file, bedtools_path):
    """Run bedtools map to sum column 4 of bedgraph over BED intervals.

    Uses zero-padding (-null 0) so ORFs with no overlap get 0 instead of '.'.
    """
    cmd = [
        bedtools_path, 'map',
        '-a', bed_file,
        '-b', bedgraph_file,
        '-c', '4',          # column 4 (value)
        '-o', 'sum',         # sum over all overlapping intervals
        '-null', '0'         # ORFs with no overlap → 0
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0 and result.stderr:
        print(f"  WARNING: bedtools map failed for {bedgraph_file}: {result.stderr.strip()[:200]}",
              file=sys.stderr, flush=True)
    return result.stdout


def parse_bedtools_output(output_text):
    """Parse bedtools map output: BED columns + sum value.

    Returns dict: {orf_id: sum_value}
    """
    result = {}
    for line in output_text.strip().split('\n'):
        if not line.strip():
            continue
        parts = line.strip().split('\t')
        if len(parts) >= 5:
            orf_id = parts[3]  # BED col4 is the ORF ID
            try:
                val = float(parts[-1])
            except (ValueError, IndexError):
                val = 0.0
            result[orf_id] = val
    return result


def process_sample(sample, bed_file, riboseqc_dir, bedtools_path):
    """Process one sample: compute P-site and coverage sums per ORF.

    Returns:
        sample_name, {orf_id: {'p_site_GSE': float, 'reads_GSE': float}}
    """
    psite_plus_file = os.path.join(riboseqc_dir, f"{sample}_P_sites_plus.bedgraph")
    psite_minus_file = os.path.join(riboseqc_dir, f"{sample}_P_sites_minus.bedgraph")
    cov_plus_file = os.path.join(riboseqc_dir, f"{sample}_coverage_plus.bedgraph")
    cov_minus_file = os.path.join(riboseqc_dir, f"{sample}_coverage_minus.bedgraph")

    orf_data = {}

    # P-site: plus + minus
    for label, bg_file in [('plus', psite_plus_file), ('minus', psite_minus_file)]:
        if os.path.exists(bg_file):
            output = run_bedtools_map(bed_file, bg_file, bedtools_path)
            vals = parse_bedtools_output(output)
            for orf_id, val in vals.items():
                if orf_id not in orf_data:
                    orf_data[orf_id] = {'p_site_GSE': 0.0, 'reads_GSE': 0.0}
                orf_data[orf_id]['p_site_GSE'] += val

    # Coverage: plus + minus
    for label, bg_file in [('plus', cov_plus_file), ('minus', cov_minus_file)]:
        if os.path.exists(bg_file):
            output = run_bedtools_map(bed_file, bg_file, bedtools_path)
            vals = parse_bedtools_output(output)
            for orf_id, val in vals.items():
                if orf_id not in orf_data:
                    orf_data[orf_id] = {'p_site_GSE': 0.0, 'reads_GSE': 0.0}
                orf_data[orf_id]['reads_GSE'] += val

    # Compute p_site_pct per ORF
    for orf_id in orf_data:
        reads = orf_data[orf_id]['reads_GSE']
        psite = orf_data[orf_id]['p_site_GSE']
        orf_data[orf_id]['p_site_pct'] = round(psite / reads, 4) if reads > 0 else 0.0
        orf_data[orf_id]['not_p_site_GSE'] = max(0.0, reads - psite)

    return sample, orf_data


def main():
    parser = argparse.ArgumentParser(
        description='Fast P-site purity computation using bedtools map')
    parser.add_argument('--bed', required=True, help='BED file of ORFs')
    parser.add_argument('--riboseqc-dir', required=True, help='Directory with RiboseQC bedgraph files')
    parser.add_argument('--output', default='psite_purity.tsv', help='Output TSV')
    parser.add_argument('--bedtools', default='bedtools', help='Path to bedtools binary')
    parser.add_argument('--workers', type=int, default=8, help='Parallel workers')
    args = parser.parse_args()

    t0 = time.time()

    # Discover samples
    samples = set()
    for f in Path(args.riboseqc_dir).iterdir():
        if f.is_file() and '_P_sites_plus.bedgraph' in f.name:
            sample = f.name.replace('_P_sites_plus.bedgraph', '')
            samples.add(sample)
    sample_list = sorted(samples)
    print(f"Found {len(sample_list)} samples", flush=True)

    # Read ORF metadata from BED for output
    print("Reading ORF metadata ...", flush=True)
    orf_meta = {}
    orf_order = []
    with open(args.bed, 'r') as f:
        for line in f:
            if not line.strip():
                continue
            parts = line.strip().split('\t')
            if len(parts) >= 6:
                orf_id = parts[3]
                orf_meta[orf_id] = {
                    'chrom': parts[0],
                    'start': int(parts[1]) + 1,  # 0-based → 1-based
                    'end': int(parts[2]),
                    'strand': parts[5]
                }
                orf_order.append(orf_id)
    print(f"  {len(orf_meta)} ORFs", flush=True)

    # Process samples in parallel
    print(f"Processing {len(sample_list)} samples with {args.workers} workers ...", flush=True)
    all_results = {}
    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        futures = {}
        for sample in sample_list:
            fut = executor.submit(process_sample, sample, args.bed,
                                  args.riboseqc_dir, args.bedtools)
            futures[fut] = sample

        for i, fut in enumerate(as_completed(futures)):
            sample, orf_data = fut.result()
            all_results[sample] = orf_data
            if (i + 1) % 5 == 0 or (i + 1) == len(sample_list):
                print(f"  {i + 1}/{len(sample_list)} samples done ({time.time() - t0:.0f}s)",
                      flush=True)

    # Write output
    print(f"Writing {args.output} ...", flush=True)
    header = ['orf_id', 'chrom', 'start', 'end', 'strand']
    for s in sample_list:
        header.append(f"{s}_p_site_GSE")
        header.append(f"{s}_reads_GSE")
        header.append(f"{s}_p_site_pct")
        header.append(f"{s}_p_site_pos")       # not computed (bedtools doesn't track positions)
        header.append(f"{s}_not_p_site_GSE")

    header.extend([
        'total_psites', 'total_reads',
        'global_p_site_pct', 'global_p_site_pos', 'global_p_site_pos_sd',
        'n_samples_with_psites'
    ])

    with open(args.output, 'w') as out:
        out.write('\t'.join(header) + '\n')

        n_written = 0
        for orf_id in orf_order:
            meta = orf_meta[orf_id]
            row = [orf_id, meta['chrom'], str(meta['start']),
                   str(meta['end']), meta['strand']]

            total_psites = 0.0
            total_reads = 0.0
            n_samples_psites = 0

            for s in sample_list:
                data = all_results.get(s, {}).get(orf_id, {})
                psite = data.get('p_site_GSE', 0.0)
                reads = data.get('reads_GSE', 0.0)
                pct = data.get('p_site_pct', 0.0)
                not_psite = data.get('not_p_site_GSE', 0.0)

                row.extend([
                    str(int(psite)), str(int(reads)),
                    f"{pct:.4f}", "0.00",   # p_site_pos not computed
                    str(int(not_psite))
                ])

                total_psites += psite
                total_reads += reads
                if psite > 0:
                    n_samples_psites += 1

            global_pct = round(total_psites / max(total_reads, 1.0), 4)
            row.extend([
                str(int(total_psites)), str(int(total_reads)),
                f"{global_pct:.4f}", "0.00", "0.00",
                str(n_samples_psites)
            ])

            out.write('\t'.join(row) + '\n')
            n_written += 1

    elapsed = time.time() - t0
    print(f"Done. {n_written} ORFs, {len(sample_list)} samples in {elapsed:.1f}s",
          flush=True)


if __name__ == '__main__':
    main()
