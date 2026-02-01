#!/usr/bin/env python3
import argparse
import sys
import os
import csv
import re
from typing import List, Dict, Tuple, Set, Optional
from collections import defaultdict

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
        self._load_gtf(gtf_file)

    def _load_gtf(self, gtf_file):
        print(f"Loading GTF: {gtf_file}...", file=sys.stderr)
        with open(gtf_file, 'r') as f:
            for line in f:
                if line.startswith('#'): continue
                parts = line.strip().split('\t')
                if len(parts) < 9: continue
                
                feature = parts[2]
                attributes = self._parse_attributes(parts[8])
                
                if 'transcript_id' not in attributes: continue
                tid = attributes['transcript_id']
                gid = attributes.get('gene_id', 'NA')
                gname = attributes.get('gene_name', gid)
                
                self.gene_map[tid] = gid
                self.gene_names[gid] = gname
                
                if tid not in self.transcripts:
                    self.transcripts[tid] = {
                        'chrom': parts[0],
                        'strand': parts[6],
                        'exons': [],
                        'cds': []
                    }
                
                start, end = int(parts[3]), int(parts[4])
                if feature == 'exon':
                    self.transcripts[tid]['exons'].append((start, end))
                elif feature == 'CDS':
                    self.transcripts[tid]['cds'].append((start, end))
        
        # Sort exons and CDS
        for tid in self.transcripts:
            self.transcripts[tid]['exons'].sort()
            self.transcripts[tid]['cds'].sort()
        print(f"Loaded {len(self.transcripts)} transcripts.", file=sys.stderr)

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

class ORFCandidate:
    def __init__(self, chrom, strand, blocks, tid, gid, tool, sample, score=None, pvalue=None, sequence=None):
        self.chrom = chrom
        self.strand = strand
        self.blocks = tuple(sorted(blocks)) # List of (start, end) tuples, 1-based
        self.tid = tid
        self.gid = gid
        self.sources = {(tool, sample)} # Set of (tool, sample) tuples
        self.score = score
        self.sequence = sequence # Extracted sequence
        self.tool_scores = {tool: score} if score is not None else {}  # Dict: tool -> score
        self.tool_pvalues = {tool: pvalue} if pvalue is not None else {}  # Dict: tool -> pvalue
        
        # Statistics from bedgraph (calculated later)
        self.total_psites = 0
        self.unique_psites = 0
        self.total_reads = 0
        self.unique_reads = 0
        
        # Calculated fields
        self.start = self.blocks[0][0]
        self.end = self.blocks[-1][1]
        self.length_nt = sum(e - s + 1 for s, e in self.blocks)
        self.length_aa = self.length_nt // 3
        
        # ID for grouping
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
    
    @property
    def pN(self):
        """P-sites per nucleotide"""
        return self.total_psites / self.length_nt if self.length_nt > 0 else 0
    
    @property
    def unique_pN(self):
        """Unique P-sites per nucleotide"""
        return self.unique_psites / self.length_nt if self.length_nt > 0 else 0

# Parsers
def parse_ribotish(file_path, gtf_index, sample_id, min_len=0):
    candidates = []
    print(f"Parsing Ribo-TISH: {file_path}", file=sys.stderr)
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
    
    return candidates

def parse_ribotricer(file_path, gtf_index, sample_id, min_len=0):
    candidates = []
    print(f"Parsing Ribotricer: {file_path}", file=sys.stderr)
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
                
                # Parsing logic for Ribotricer ORF_ID: tid_start_end
                match = re.search(r'_(\d+)_(\d+)$', orf_id)
                if match:
                    t_start = int(match.group(1))
                    t_end = int(match.group(2))
                    
                    # Convert 1-based to 0-based
                    t_start -= 1
                    
                    if t_end - t_start < min_len * 3: continue
                    
                    c, s, b = gtf_index.get_genomic_blocks(tid, t_start, t_end, feature_type='exon')
                    if b:
                        chrom, strand, blocks = c, s, b
                
                if blocks:
                    gid = gtf_index.gene_map.get(tid, 'NA')
                    cand = ORFCandidate(chrom, strand, blocks, tid, gid, 'Ribotricer', sample_id, score=score)
                    candidates.append(cand)
    except Exception as e:
        print(f"Error parsing Ribotricer file: {e}", file=sys.stderr)

    return candidates

def parse_orfquant(file_path, gtf_index, sample_id, min_len=0):
    candidates = []
    print(f"Parsing ORFquant: {file_path}", file=sys.stderr)
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
                        if 'P_sites' in current_attrs:
                            try:
                                score = float(current_attrs['P_sites'])
                            except ValueError:
                                pass
                        
                        length_nt = sum(e-s+1 for s,e in current_blocks)
                        if length_nt // 3 >= min_len:
                            cand = ORFCandidate(chrom, strand, current_blocks, tid, gid, 'ORFquant', sample_id, score=score)
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
            if 'P_sites' in current_attrs:
                try:
                    score = float(current_attrs['P_sites'])
                except ValueError:
                    pass
            
            length_nt = sum(e-s+1 for s,e in current_blocks)
            if length_nt // 3 >= min_len:
                cand = ORFCandidate(chrom, strand, current_blocks, tid, gid, 'ORFquant', sample_id, score=score)
                candidates.append(cand)

    except Exception as e:
        print(f"Error parsing ORFquant file: {e}", file=sys.stderr)
    
    return candidates

def validate_sequence(cand, genome_fasta):
    try:
        seq_parts = []
        for s, e in cand.blocks:
            # pyfaidx 0-based indexing
            seq = genome_fasta[cand.chrom][s-1:e]
            seq_parts.append(str(seq))
        
        full_seq = "".join(seq_parts)
        seq_obj = Seq(full_seq)
        
        if cand.strand == '-':
            seq_obj = seq_obj.reverse_complement()
            
        cand.sequence = str(seq_obj)
        
    except Exception as e:
        cand.sequence = "N" * cand.length_nt

def count_psites_in_region(bedgraph_file, chrom, start, end):
    """
    Count P-sites in a genomic region from bedgraph file
    bedgraph format: chrom start end value
    Returns: total count
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
    Calculate statistics from RiboseQC bedgraph files
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

def main():
    parser = argparse.ArgumentParser(description="Unify ORF predictions from multiple tools")
    parser.add_argument("--ribotish", nargs='+', help="Ribo-TISH output files")
    parser.add_argument("--ribotricer", nargs='+', help="Ribotricer output files")
    parser.add_argument("--orfquant", nargs='+', help="ORFquant GTF output files")
    parser.add_argument("--gtf", required=True, help="Reference GTF file for coordinate mapping")
    parser.add_argument("--fasta", required=True, help="Genome FASTA file for validation")
    parser.add_argument("--output", required=True, help="Output prefix")
    parser.add_argument("--min_len", type=int, default=10, help="Minimum amino acid length")
    parser.add_argument("--bedgraph-dir", help="Directory containing RiboseQC bedgraph files (optional)")
    parser.add_argument("--sample-list", help="Comma-separated list of sample names for bedgraph stats (optional)")
    
    args = parser.parse_args()
    
    gtf_index = GTFIndex(args.gtf)
    
    print(f"Loading Genome FASTA: {args.fasta}...", file=sys.stderr)
    genome_fasta = Fasta(args.fasta)
    
    all_candidates = []
    
    if args.ribotish:
        for f in args.ribotish:
            sid = os.path.basename(f).split('.')[0].replace('_pred', '')
            all_candidates.extend(parse_ribotish(f, gtf_index, sid, args.min_len))
            
    if args.ribotricer:
        for f in args.ribotricer:
            sid = os.path.basename(f).split('.')[0].replace('_translating_ORFs', '')
            all_candidates.extend(parse_ribotricer(f, gtf_index, sid, args.min_len))
            
    if args.orfquant:
        for f in args.orfquant:
            sid = os.path.basename(f).split('.')[0].replace('_Detected_ORFs', '')
            all_candidates.extend(parse_orfquant(f, gtf_index, sid, args.min_len))

    print(f"Total raw candidates: {len(all_candidates)}", file=sys.stderr)
    
    merged_candidates = {} 
    for cand in all_candidates:
        if cand.id_key in merged_candidates:
            merged_candidates[cand.id_key].merge(cand)
        else:
            merged_candidates[cand.id_key] = cand
    
    print(f"Unique candidates after merging: {len(merged_candidates)}", file=sys.stderr)
    
    final_list = list(merged_candidates.values())
    
    # Calculate statistics from bedgraphs if provided
    if args.bedgraph_dir and args.sample_list:
        sample_list = args.sample_list.split(',')
        print(f"Calculating statistics from bedgraphs for {len(sample_list)} samples...", file=sys.stderr)
        for cand in final_list:
            calculate_statistics_from_bedgraphs(cand, args.bedgraph_dir, sample_list)
    
    # Validate sequences
    for cand in final_list:
        validate_sequence(cand, genome_fasta)
    
    # Write metadata.tsv with extended columns
    with open(f"{args.output}.metadata.tsv", 'w') as out:
        header = ["orf_id", "chrom", "strand", "start", "end", "length_aa", "exon_blocks", 
                  "gene_id", "transcript_id", "tools", "samples", 
                  "tool_scores", "tool_pvalues", 
                  "total_reads", "unique_reads", "total_psites", "unique_psites", "pN", "unique_pN",
                  "sequence"]
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
            
            row = [orf_id, cand.chrom, cand.strand, str(cand.start), str(cand.end), str(cand.length_aa), 
                   blocks_str, cand.gid, cand.tid, tools, samples, 
                   tool_scores_str, tool_pvalues_str,
                   str(cand.total_reads), str(cand.unique_reads), 
                   str(cand.total_psites), str(cand.unique_psites),
                   f"{cand.pN:.6f}", f"{cand.unique_pN:.6f}",
                   cand.sequence]
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
