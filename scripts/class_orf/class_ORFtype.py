#!/usr/bin/env python3
import sys
import os
import argparse
import csv
csv.field_size_limit(sys.maxsize)
from typing import List, Dict, Tuple, Any
from multiprocessing import Pool, cpu_count
from functools import partial

# --- 1. Parsing & Data Loading ---

def parse_coords(coord_str: str) -> List[Tuple[int, int]]:
    """Parses 'start-end,start-end' string into list of tuples."""
    exons = []
    if not coord_str or coord_str.lower() == 'na': return exons
    for part in coord_str.split(','):
        try:
            start, end = map(int, part.split('-'))
            exons.append((start, end))
        except ValueError:
            continue
    return exons

def load_gene_level_cds_from_gtf(gtf_file: str) -> Dict[str, Dict[str, Any]]:
    """
    Builds gene-level CDS reference from GTF.
    Merges all CDS intervals for a gene into a canonical set.
    Returns: gene_id -> {'strand': str, 'exons': [(s,e), ...]}
    """
    print(f"Loading GTF to build Gene-level CDS: {gtf_file}...", file=sys.stderr)
    gene_cds = {} # gid -> {strand, intervals: set of (s,e)}
    
    with open(gtf_file, 'r') as f:
        for line in f:
            if line.startswith('#'): continue
            parts = line.strip().split('\t')
            if len(parts) < 9: continue
            
            if parts[2] != 'CDS': continue
            
            # Parse attributes
            attrs = {}
            for p in parts[8].split(';'):
                p = p.strip()
                if not p: continue
                kv = p.split(' ')
                attrs[kv[0]] = ' '.join(kv[1:]).strip('"')
            
            gid = attrs.get('gene_id')
            if not gid: continue
            
            strand = parts[6]
            start, end = int(parts[3]), int(parts[4])
            
            if gid not in gene_cds:
                gene_cds[gid] = {'strand': strand, 'intervals': set()}
            
            gene_cds[gid]['intervals'].add((start, end))
    
    # Convert sets to sorted lists
    final_cds = {}
    for gid, data in gene_cds.items():
        sorted_intervals = sorted(list(data['intervals']))
        # Merge overlapping intervals? 
        # For gene-level CDS, usually we want the union of all CDS exons.
        merged = []
        if sorted_intervals:
            curr_s, curr_e = sorted_intervals[0]
            for next_s, next_e in sorted_intervals[1:]:
                if next_s <= curr_e + 1: # Overlap or adjacent
                    curr_e = max(curr_e, next_e)
                else:
                    merged.append((curr_s, curr_e))
                    curr_s, curr_e = next_s, next_e
            merged.append((curr_s, curr_e))
            
        final_cds[gid] = {
            'strand': data['strand'],
            'exons': merged
        }
    
    print(f"Gene-level CDS built for {len(final_cds)} genes.", file=sys.stderr)
    return final_cds

# --- 2. Classification Logic ---

def check_overlap(exons1: List[Tuple[int, int]], exons2: List[Tuple[int, int]]) -> bool:
    for s1, e1 in exons1:
        for s2, e2 in exons2:
            if max(s1, s2) <= min(e1, e2): return True
    return False

def classify_orf(orf: Dict, merged_cds: Dict) -> str:
    """
    Classifies an ORF based on its overlap with Gene-level CDS.
    Types: annotated, truncated, extension, isoform, novel, etc.
    This logic mimics the original R source code logic mentioned in previous version.
    """
    if not merged_cds: return "novel"
    
    if orf['strand'] != merged_cds['strand']: return "novel_antisenese"
    
    orf_exons = orf['exons']
    cds_exons = merged_cds['exons']
    
    # Check overlap
    if check_overlap(orf_exons, cds_exons):
        strand = orf['strand']
        
        # Define genomic boundaries
        cds_min_coord = min(s for s, e in cds_exons)
        cds_max_coord = max(e for s, e in cds_exons)
        
        # Gene Start/Stop depend on strand
        gen_sta, gen_sto = (cds_min_coord, cds_max_coord) if strand == '+' else (cds_max_coord, cds_min_coord)
        
        orf_min_coord = min(s for s, e in orf_exons)
        orf_max_coord = max(e for s, e in orf_exons)
        orf_sta, orf_sto = (orf_min_coord, orf_max_coord) if strand == '+' else (orf_max_coord, orf_min_coord)
        
        # Calculate shifts (positive if ORF is upstream/longer)
        if strand == '+':
            shift_start = gen_sta - orf_sta
            shift_stop = orf_sto - gen_sto
        else:
            shift_start = orf_sta - gen_sta # For minus, higher coord is start. 
            # If ORF start > Gene start, ORF is upstream.
            shift_stop = gen_sto - orf_sto
        
        # Classification rules
        if shift_start == 0 and shift_stop == 0:
            return "annotated" # Matches gene boundaries exactly
        elif shift_start == 0 and shift_stop < 0:
            return "truncated" # Same start, earlier stop
        elif shift_start == 0 and shift_stop > 0:
            return "extension" # Same start, later stop
        elif shift_start > 0:
            # Starts upstream
            if shift_stop == 0:
                return "extension" # Earlier start, same stop
            else:
                # Upstream start, different stop. 
                # Could be isoform or novel_upstream.
                return "isoform" 
        elif shift_start < 0:
            # Starts downstream
            if shift_stop == 0:
                return "truncated" # Later start, same stop (N-terminal truncation)
            elif shift_stop < 0:
                # Later start, earlier stop (nested)
                return "isoform" # or internal
            else:
                # Later start, later stop
                return "isoform"
        
        return "isoform" # Catch all for overlaps
        
    else:
        # No overlap
        # Check relative position
        strand = orf['strand']
        cds_min = min(s for s, e in cds_exons)
        cds_max = max(e for s, e in cds_exons)
        orf_min = min(s for s, e in orf_exons)
        orf_max = max(e for s, e in orf_exons)
        
        if strand == '+':
            if orf_max < cds_min: return "novel_upstream"
            if orf_min > cds_max: return "novel_downstream"
        else:
            if orf_min > cds_max: return "novel_upstream" # 5' is high coord
            if orf_max < cds_min: return "novel_downstream"
            
        return "novel"

def process_chunk(chunk_orfs: List[Dict], cds_data: Dict[str, Dict]) -> List[Tuple[str, str]]:
    results = []
    for orf in chunk_orfs:
        gid = orf['gene_id']
        merged_cds = cds_data.get(gid)
        category = classify_orf(orf, merged_cds)
        results.append((orf['orf_id'], category))
    return results

# --- 3. Main ---

CODON_TABLE = {
    'TTT':'F','TTC':'F','TTA':'L','TTG':'L',
    'CTT':'L','CTC':'L','CTA':'L','CTG':'L',
    'ATT':'I','ATC':'I','ATA':'I','ATG':'M',
    'GTT':'V','GTC':'V','GTA':'V','GTG':'V',
    'TCT':'S','TCC':'S','TCA':'S','TCG':'S',
    'CCT':'P','CCC':'P','CCA':'P','CCG':'P',
    'ACT':'T','ACC':'T','ACA':'T','ACG':'T',
    'GCT':'A','GCC':'A','GCA':'A','GCG':'A',
    'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*',
    'CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
    'AAT':'N','AAC':'N','AAA':'K','AAG':'K',
    'GAT':'D','GAC':'D','GAA':'E','GAG':'E',
    'TGT':'C','TGC':'C','TGA':'*','TGG':'W',
    'CGT':'R','CGC':'R','CGA':'R','CGG':'R',
    'AGT':'S','AGC':'S','AGA':'R','AGG':'R',
    'GGT':'G','GGC':'G','GGA':'G','GGG':'G',
}

def translate_nt(nt_seq: str) -> str:
    nt_seq = nt_seq.upper().replace('U', 'T')
    aa = []
    for i in range(0, len(nt_seq) - 2, 3):
        codon = nt_seq[i:i+3]
        aa.append(CODON_TABLE.get(codon, 'X'))
    return ''.join(aa)


def write_extra_outputs(orfs: List[Dict], classification_map: Dict[str, str], prefix: str, args) -> None:
    """Write BED12, GTF, nucleotide FASTA, protein FASTA, logs, and per-sample stats."""

    # Collect all unique samples across all ORFs
    all_samples: List[str] = []
    seen: set = set()
    for orf in orfs:
        for s in orf.get('samples', '').split(','):
            s = s.strip()
            if s and s not in seen:
                seen.add(s)
                all_samples.append(s)
    all_samples.sort()

    with open(f"{prefix}.orfs.bed", 'w') as bed_f, \
         open(f"{prefix}.orfs.gtf", 'w') as gtf_f, \
         open(f"{prefix}.orfs.fa", 'w') as fa_nt_f, \
         open(f"{prefix}.orfs.pep.fa", 'w') as fa_pep_f, \
         open(f"{prefix}.orfs.out", 'w') as out_f:

        # Write out header
        out_header = ['orf_id', 'chrom', 'start', 'end', 'strand',
                      'gene_id', 'transcript_id', 'length_aa',
                      'tools', 'n_samples', 'orf_type_category'] + all_samples
        out_f.write('\t'.join(out_header) + '\n')

        for orf in orfs:
            orf_id = orf.get('orf_id', 'NA')
            chrom = orf.get('chrom', 'NA')
            strand = orf.get('strand', '+')
            gene_id = orf.get('gene_id', 'NA')
            transcript_id = orf.get('transcript_id', 'NA')
            length_aa = orf.get('length_aa', '0')
            tools = orf.get('tools', 'NA')
            sequence = orf.get('sequence', '')
            classification = classification_map.get(orf_id, 'unknown')

            # Parse exon blocks (1-based inclusive genomic coords)
            blocks = parse_coords(orf.get('exon_blocks', ''))
            if not blocks:
                try:
                    start_1 = int(orf.get('start', 0))
                    end_1 = int(orf.get('end', 0))
                    blocks = [(start_1, end_1)]
                except (ValueError, TypeError):
                    continue

            # BED12: chromStart is 0-based, chromEnd is 0-based exclusive
            chrom_start = min(s for s, e in blocks) - 1
            chrom_end = max(e for s, e in blocks)
            block_count = len(blocks)
            block_sizes = ','.join(str(e - s + 1) for s, e in blocks)
            block_starts = ','.join(str(s - 1 - chrom_start) for s, e in blocks)

            bed_f.write(f"{chrom}\t{chrom_start}\t{chrom_end}\t{orf_id}\t0\t{strand}\t"
                        f"{chrom_start}\t{chrom_end}\t0,0,0\t{block_count}\t{block_sizes}\t{block_starts}\n")

            # GTF: CDS feature per exon block (1-based inclusive)
            for s, e in blocks:
                attrs = (f'gene_id "{gene_id}"; transcript_id "{transcript_id}"; '
                         f'orf_id "{orf_id}"; orf_type_category "{classification}";')
                gtf_f.write(f"{chrom}\torf_type\tCDS\t{s}\t{e}\t.\t{strand}\t0\t{attrs}\n")

            # Nucleotide FASTA
            if sequence:
                fa_nt_f.write(f">{orf_id} {chrom}:{chrom_start+1}-{chrom_end}({strand})\n{sequence}\n")
                pep = translate_nt(sequence)
                fa_pep_f.write(f">{orf_id} {chrom}:{chrom_start+1}-{chrom_end}({strand})\n{pep}\n")

            # Per-sample stats row
            sample_set = set(s.strip() for s in orf.get('samples', '').split(',') if s.strip())
            n_samples = len(sample_set)
            sample_cols = [('1' if s in sample_set else '0') for s in all_samples]
            out_row = [orf_id, chrom, orf.get('start',''), orf.get('end',''), strand,
                       gene_id, transcript_id, length_aa,
                       tools, str(n_samples), classification] + sample_cols
            out_f.write('\t'.join(out_row) + '\n')

    # Logs file
    with open(f"{prefix}.logs", 'w') as log_f:
        log_f.write(f"#input metadata: {args.input}\n")
        log_f.write(f"#reference gtf: {args.gtf}\n")
        log_f.write(f"#cpus: {args.cpus}\n")
        log_f.write(f"#total_orfs: {len(orfs)}\n")
        log_f.write(f"#total_samples: {len(all_samples)}\n")
        log_f.write(f"#samples: {';'.join(all_samples)}\n")

    print(f"Extra outputs written with prefix: {prefix}", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(description="Classify ORFs based on Gene-level CDS overlap (ORFtype)")
    parser.add_argument("--input", required=True, help="Unified Metadata TSV file")
    parser.add_argument("--gtf", required=True, help="Reference GTF file (to build gene CDS)")
    parser.add_argument("--output", required=True, help="Output file path")
    parser.add_argument("--output_prefix", default="orftype_results",
                        help="Prefix for extra output files (BED/GTF/FA/logs/out). Default: orftype_results")
    parser.add_argument("--cpus", type=int, default=1, help="Number of CPUs")
    
    args = parser.parse_args()
    
    # 1. Load CDS Data
    cds_data = load_gene_level_cds_from_gtf(args.gtf)
    
    # 2. Load ORFs
    print(f"Loading ORFs from {args.input}...", file=sys.stderr)
    orfs = []
    original_header = None
    with open(args.input, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        original_header = reader.fieldnames
        for row in reader:
            row['exons'] = parse_coords(row['exon_blocks'])
            orfs.append(row)
    
    print(f"Loaded {len(orfs)} ORFs.", file=sys.stderr)
    
    # 3. Parallel Processing
    num_processes = min(args.cpus, cpu_count())
    if num_processes < 1: num_processes = 1
    chunk_size = (len(orfs) // num_processes) + 1
    chunks = [orfs[i:i + chunk_size] for i in range(0, len(orfs), chunk_size)]
    
    print(f"Classifying with {num_processes} processes...", file=sys.stderr)
    
    with Pool(processes=num_processes) as pool:
        func = partial(process_chunk, cds_data=cds_data)
        results_list = pool.map(func, chunks)
    
    # Flatten results and create a mapping
    classification_map = {}
    for sublist in results_list:
        for orf_id, cat in sublist:
            classification_map[orf_id] = cat
    
    # 4. Write Output - Preserve ALL original columns and add classification
    print(f"Writing results to {args.output}...", file=sys.stderr)
    
    # Build output header: original columns + orf_type_category
    if original_header:
        output_header = list(original_header) + ['orf_type_category']
    else:
        output_header = ['orf_id', 'orf_type_category']
    
    with open(args.output, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=output_header, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        
        for orf in orfs:
            orf_id = orf.get('orf_id', 'NA')
            orf['orf_type_category'] = classification_map.get(orf_id, 'unknown')
            # Remove the temporary 'exons' field before writing
            orf_copy = {k: v for k, v in orf.items() if k != 'exons'}
            writer.writerow(orf_copy)

    # 5. Write extra outputs (BED, GTF, FA, logs, per-sample out)
    write_extra_outputs(orfs, classification_map, args.output_prefix, args)

    print("Done.", file=sys.stderr)

if __name__ == "__main__":
    main()
