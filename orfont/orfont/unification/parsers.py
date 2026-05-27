"""ORF prediction parsers for Ribo-TISH, Ribotricer, ORFquant, and RiboCode.

Ported from the nf-core/riboseq pipeline's unify_orf_predictions.py.
Each parser reads a tool-specific output format and returns a list of
ORFCandidate objects using the shared GTFIndex for coordinate resolution.
"""

import gzip
import math
import os
import re
import sys

from orfont.core.models import ORFCandidate


def _open(path, mode='r'):
    if path.endswith('.gz'):
        return gzip.open(path, mode + 't')
    return open(path, mode)


def infer_sample_id_from_prediction_path(file_path, tool_suffix):
    """Recover sample ID from a prediction filename."""
    name = os.path.basename(file_path)
    if name.endswith(tool_suffix):
        return name[:-len(tool_suffix)]
    return os.path.splitext(name)[0]


# ---------------------------------------------------------------------------
# Ribo-TISH parser
# ---------------------------------------------------------------------------

def parse_ribotish(file_path, gtf_index, sample_id, min_len=0,
                   exclude_tistypes=None, atg_only=False):
    candidates = []
    print(f"Parsing Ribo-TISH: {file_path}", file=sys.stderr)

    try:
        with _open(file_path, 'r') as f:
            first_line = f.readline().strip()
            if first_line.startswith('#') and (
                'placeholder' in first_line.lower()
                or 'no orfs' in first_line.lower()
                or 'insufficient' in first_line.lower()
            ):
                print("  -> Placeholder file detected, skipping", file=sys.stderr)
                return []
    except Exception:
        pass

    filtered_tis = 0
    filtered_atg = 0

    try:
        with _open(file_path, 'r') as f:
            header = f.readline().strip().split('\t')
            col_map = {name: i for i, name in enumerate(header)}

            has_start_stop = 'Tid' in col_map and 'Start' in col_map and 'Stop' in col_map
            has_genomepos = 'Tid' in col_map and 'GenomePos' in col_map

            if not has_start_stop and not has_genomepos:
                print("Warning: Ribo-TISH file missing required columns. Skipping.",
                      file=sys.stderr)
                return []

            for line in f:
                parts = line.strip().split('\t')
                tid = parts[col_map['Tid']]

                if exclude_tistypes and 'TisType' in col_map:
                    tis_type = parts[col_map['TisType']] if col_map['TisType'] < len(parts) else ''
                    if tis_type in exclude_tistypes:
                        filtered_tis += 1
                        continue

                if atg_only and 'StartCodon' in col_map:
                    start_codon = parts[col_map['StartCodon']] if col_map['StartCodon'] < len(parts) else ''
                    if start_codon.upper() != 'ATG':
                        filtered_atg += 1
                        continue

                score = None
                pvalue = None
                if 'TisPvalue' in col_map:
                    try:
                        pval = float(parts[col_map['TisPvalue']])
                        if pval > 0:
                            score = -math.log10(pval)
                            pvalue = pval
                    except (ValueError, OverflowError):
                        pass

                t_start = None
                t_stop = None
                chrom = None
                strand = None

                if has_start_stop:
                    try:
                        t_start = int(parts[col_map['Start']])
                        t_stop = int(parts[col_map['Stop']])
                    except ValueError:
                        continue

                    length_nt = t_stop - t_start
                    if length_nt // 3 < min_len:
                        continue

                    chrom, strand, blocks = gtf_index.get_genomic_blocks(
                        tid, t_start, t_stop, feature_type='exon')

                elif has_genomepos:
                    genome_pos = parts[col_map['GenomePos']]
                    match = re.search(r'(.+):(\d+)-(\d+):([+-])', genome_pos)
                    if match:
                        chrom = gtf_index.resolve_chrom(match.group(1))
                        t_start = int(match.group(2))
                        t_stop = int(match.group(3))
                        strand = match.group(4)

                        length_nt = t_stop - t_start + 1
                        if length_nt // 3 < min_len:
                            continue

                        blocks = [(t_start, t_stop)]
                    else:
                        continue
                else:
                    continue

                if blocks:
                    gid = gtf_index.gene_map.get(tid, 'NA')
                    cand = ORFCandidate(chrom, strand, blocks, tid, gid,
                                        'Ribo-TISH', sample_id,
                                        score=score, pvalue=pvalue)
                    candidates.append(cand)
    except Exception as e:
        print(f"Error parsing Ribo-TISH file: {e}", file=sys.stderr)

    if filtered_tis or filtered_atg:
        print(f"  -> Filtered: {filtered_tis} by TisType, {filtered_atg} by start codon",
              file=sys.stderr)
    return candidates


# ---------------------------------------------------------------------------
# Ribotricer parser
# ---------------------------------------------------------------------------

def parse_ribotricer(file_path, gtf_index, sample_id, min_len=0,
                     exclude_tistypes=None, atg_only=False):
    candidates = []
    print(f"Parsing Ribotricer: {file_path}", file=sys.stderr)

    try:
        with _open(file_path, 'r') as f:
            first_line = f.readline().strip()
            if first_line.startswith('#') and (
                'placeholder' in first_line.lower()
                or 'no orfs' in first_line.lower()
                or 'insufficient' in first_line.lower()
            ):
                print("  -> Placeholder file detected, skipping", file=sys.stderr)
                return []
    except Exception:
        pass

    filtered_tis = 0
    filtered_atg = 0

    try:
        with _open(file_path, 'r') as f:
            header_line = f.readline()
            while header_line.startswith('#'):
                header_line = f.readline()
            header = header_line.strip().split('\t')
            col_map = {name: i for i, name in enumerate(header)}

            if 'transcript_id' not in col_map or 'ORF_ID' not in col_map:
                print("Warning: Missing required columns in Ribotricer file.", file=sys.stderr)
                return []

            for line in f:
                parts = line.strip().split('\t')
                tid = parts[col_map['transcript_id']]
                orf_id = parts[col_map['ORF_ID']]

                if exclude_tistypes and 'ORF_type' in col_map:
                    orf_type = parts[col_map['ORF_type']] if col_map['ORF_type'] < len(parts) else ''
                    if orf_type in exclude_tistypes:
                        filtered_tis += 1
                        continue

                if atg_only and 'start_codon' in col_map:
                    start_codon = parts[col_map['start_codon']] if col_map['start_codon'] < len(parts) else ''
                    if start_codon.upper() != 'ATG':
                        filtered_atg += 1
                        continue

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

                match = re.search(r'_(\d+)_(\d+)(?:_\d+)?$', orf_id)
                if match:
                    t_start_raw = int(match.group(1))
                    t_end_raw = int(match.group(2))
                    t_start = t_start_raw
                    t_end = t_end_raw

                    t_start -= 1  # 1-based to 0-based

                    if t_end - t_start < min_len * 3:
                        continue

                    c, s, b = gtf_index.get_genomic_blocks(
                        tid, t_start, t_end, feature_type='exon')
                    if b:
                        chrom, strand, blocks = c, s, b
                    elif chrom and strand and t_end_raw >= t_start_raw:
                        length_nt = t_end_raw - t_start_raw + 1
                        if length_nt // 3 >= min_len:
                            blocks = [(t_start_raw, t_end_raw)]

                if blocks:
                    gid = gtf_index.gene_map.get(tid, 'NA')
                    cand = ORFCandidate(chrom, strand, blocks, tid, gid,
                                        'Ribotricer', sample_id, score=score)
                    candidates.append(cand)
    except Exception as e:
        print(f"Error parsing Ribotricer file: {e}", file=sys.stderr)

    if filtered_tis or filtered_atg:
        print(f"  -> Filtered: {filtered_tis} by ORF_type, {filtered_atg} by start codon",
              file=sys.stderr)
    return candidates


# ---------------------------------------------------------------------------
# ORFquant parser
# ---------------------------------------------------------------------------

def parse_orfquant(file_path, gtf_index, sample_id, min_len=0,
                   exclude_tistypes=None, atg_only=False):
    candidates = []
    print(f"Parsing ORFquant: {file_path}", file=sys.stderr)

    try:
        with _open(file_path, 'r') as f:
            first_line = f.readline().strip()
            if first_line.startswith('#') and (
                'placeholder' in first_line.lower()
                or 'no orfs' in first_line.lower()
                or 'insufficient' in first_line.lower()
            ):
                print("  -> Placeholder file detected, skipping", file=sys.stderr)
                return []
    except Exception:
        pass

    try:
        current_orf = None
        current_blocks = []
        current_attrs = {}

        with _open(file_path, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) < 9:
                    continue

                feature = parts[2]
                if feature != 'CDS':
                    continue

                attrs = {}
                for p in parts[8].split(';'):
                    p = p.strip()
                    if not p:
                        continue
                    kv = p.split(' ')
                    attrs[kv[0]] = ' '.join(kv[1:]).strip('"')

                orf_id = attrs.get('ORF_id') or attrs.get('transcript_id')
                if not orf_id:
                    continue

                if current_orf != orf_id:
                    if current_orf and current_blocks:
                        _finalize_orfquant_candidate(
                            candidates, current_attrs, current_blocks,
                            gtf_index, sample_id, min_len)

                    current_orf = orf_id
                    current_blocks = []
                    current_attrs = attrs
                    current_attrs['chrom'] = gtf_index.resolve_chrom(parts[0])
                    current_attrs['strand'] = parts[6]

                current_blocks.append((int(parts[3]), int(parts[4])))

        if current_orf and current_blocks:
            _finalize_orfquant_candidate(
                candidates, current_attrs, current_blocks,
                gtf_index, sample_id, min_len)

    except Exception as e:
        print(f"Error parsing ORFquant file: {e}", file=sys.stderr)

    return candidates


def _finalize_orfquant_candidate(candidates, attrs, blocks, gtf_index, sample_id, min_len):
    """Create an ORFCandidate from accumulated ORFquant CDS blocks."""
    chrom = attrs.get('chrom')
    strand = attrs.get('strand')
    tid = attrs.get('transcript_id', 'NA')
    gid = attrs.get('gene_id', 'NA')

    score = None
    psites_count = 0
    if 'P_sites' in attrs:
        try:
            psites_count = int(float(attrs['P_sites']))
            score = float(attrs['P_sites'])
        except ValueError:
            pass

    unique_psites_count = 0
    if 'P_sites_uniq' in attrs:
        try:
            unique_psites_count = int(float(attrs['P_sites_uniq']))
        except ValueError:
            pass

    length_nt = sum(e - s + 1 for s, e in blocks)
    if length_nt // 3 >= min_len:
        cand = ORFCandidate(chrom, strand, blocks, tid, gid,
                            'ORFquant', sample_id, score=score)
        cand.total_psites = psites_count
        cand.unique_psites = unique_psites_count if unique_psites_count > 0 else psites_count
        candidates.append(cand)


# ---------------------------------------------------------------------------
# RiboCode parser
# ---------------------------------------------------------------------------

def parse_ribocode(file_path, gtf_index, sample_id, min_len=0,
                   exclude_tistypes=None, atg_only=False):
    candidates = []
    print(f"Parsing RiboCode: {file_path}", file=sys.stderr)

    try:
        with _open(file_path, 'r') as f:
            first_line = f.readline().strip()
            if first_line.startswith('#') and (
                'placeholder' in first_line.lower()
                or 'no orfs' in first_line.lower()
                or 'insufficient' in first_line.lower()
            ):
                print("  -> Placeholder file detected, skipping", file=sys.stderr)
                return []
            if first_line.startswith('FAILED'):
                print("  -> Failure marker detected, skipping", file=sys.stderr)
                return []
    except Exception:
        pass

    def _parse_attrs(attr_str):
        attrs = {}
        for item in attr_str.split(';'):
            item = item.strip()
            if not item:
                continue
            parts = item.split(' ', 1)
            if len(parts) == 2:
                attrs[parts[0]] = parts[1].strip().strip('"')
        return attrs

    def _score_from_pvalue(pvalue):
        if pvalue is None:
            return None
        try:
            pval = float(pvalue)
            if pval > 0:
                return -math.log10(pval)
        except (TypeError, ValueError, OverflowError):
            return None
        return None

    def _load_txt_metrics(path):
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
                        'score': _score_from_pvalue(pvalue),
                    }
        except Exception as e:
            print(f"Warning: could not read RiboCode sidecar metrics {path}: {e}",
                  file=sys.stderr)
        return metrics

    def _sidecar_txt_path(path):
        root, ext = os.path.splitext(path)
        if ext in ('.gtf', '.bed'):
            return root + '.txt'
        return path if ext == '.txt' else None

    metrics = _load_txt_metrics(_sidecar_txt_path(file_path))

    # GTF format
    if file_path.endswith('.gtf'):
        grouped = {}
        try:
            with _open(file_path, 'r') as f:
                for line in f:
                    if line.startswith('#') or not line.strip():
                        continue
                    parts = line.rstrip('\n').split('\t')
                    if len(parts) < 9:
                        continue
                    attrs = _parse_attrs(parts[8])
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
                    rec['chrom'], rec['strand'], blocks, tid, gid,
                    'RiboCode', sample_id,
                    score=metric.get('score'),
                    pvalue=metric.get('pvalue'),
                )
                candidates.append(cand)
        except Exception as e:
            print(f"Error parsing RiboCode GTF file: {e}", file=sys.stderr)
        return candidates

    # Tab-delimited format
    try:
        with _open(file_path, 'r') as f:
            header = f.readline().strip().split('\t')
            col_map = {name: i for i, name in enumerate(header)}
            required = {'ORF_ID', 'chrom', 'strand', 'ORF_gstart', 'ORF_gstop'}
            if not required.issubset(col_map):
                print("Warning: RiboCode file missing required columns. Skipping.",
                      file=sys.stderr)
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
                tid = (parts[col_map['transcript_id']]
                       if 'transcript_id' in col_map and col_map['transcript_id'] < len(parts)
                       else 'NA')
                gid = (parts[col_map['gene_id']]
                       if 'gene_id' in col_map and col_map['gene_id'] < len(parts)
                       else gtf_index.gene_map.get(tid, 'NA'))
                pvalue = None
                for col in ('adjusted_pval', 'pval_combined'):
                    if col in col_map and col_map[col] < len(parts):
                        try:
                            pvalue = float(parts[col_map[col]])
                            break
                        except ValueError:
                            pass
                cand = ORFCandidate(
                    chrom, strand, [(start, end)], tid, gid,
                    'RiboCode', sample_id,
                    score=_score_from_pvalue(pvalue),
                    pvalue=pvalue,
                )
                candidates.append(cand)
    except Exception as e:
        print(f"Error parsing RiboCode file: {e}", file=sys.stderr)

    return candidates
