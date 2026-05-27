"""GTF loader — streams GTF annotation directly into DuckDB tables.

Replaces the GTFIndex in-memory dict with DuckDB-backed storage,
enabling SQL-based attribute queries and reducing memory footprint.

Provides two loading modes:
  - load_gtf(): full GTF load (genes + transcripts with exon/CDS JSON)
  - load_gene_names(): lightweight load of gene_id → gene_name only
"""

import json
import logging

from orfont.core.db import get_connection, append_rows

logger = logging.getLogger(__name__)


def _parse_attributes(attr_str):
    """Parse GTF attribute string into a dict."""
    attrs = {}
    for p in attr_str.split(';'):
        p = p.strip()
        if not p:
            continue
        parts = p.split(' ')
        if len(parts) >= 2:
            key = parts[0]
            val = ' '.join(parts[1:]).strip('"')
            attrs[key] = val
    return attrs


def load_gtf(gtf_file, con=None):
    """Load a GTF file into DuckDB gtf_genes and gtf_transcripts tables.

    Parses the GTF line-by-line, collecting genes and transcript structures,
    then bulk-inserts into DuckDB via append_rows.

    Args:
        gtf_file: path to GTF file
        con: DuckDB connection (uses thread-local default if None)

    Returns:
        dict with counts: {'genes': int, 'transcripts': int}
    """
    if con is None:
        con = get_connection()

    gene_rows = []
    tx_data = {}  # tid -> {chrom, strand, exons: [(s,e)], cds: [(s,e)], gid}

    logger.info("Loading GTF: %s", gtf_file)
    with open(gtf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue

            feature = parts[2]
            chrom = parts[0]
            strand = parts[6]
            start, end = int(parts[3]), int(parts[4])
            attrs = _parse_attributes(parts[8])

            if feature == 'gene':
                gid = attrs.get('gene_id', 'NA')
                gname = attrs.get('gene_name', gid)
                biotype = attrs.get('gene_biotype', attrs.get('gene_type', ''))
                gene_rows.append((gid, gname, chrom, strand, start, end, biotype))

            elif feature in ('exon', 'CDS'):
                tid = attrs.get('transcript_id')
                if not tid:
                    continue
                gid = attrs.get('gene_id', 'NA')

                if tid not in tx_data:
                    tx_data[tid] = {
                        'chrom': chrom,
                        'strand': strand,
                        'gene_id': gid,
                        'exons': [],
                        'cds': [],
                    }
                if feature == 'exon':
                    tx_data[tid]['exons'].append((start, end))
                elif feature == 'CDS':
                    tx_data[tid]['cds'].append((start, end))

    # Build transcript rows: sort exons/CDS, serialize as JSON
    tx_rows = []
    for tid, data in tx_data.items():
        data['exons'].sort()
        data['cds'].sort()
        exons_json = json.dumps(data['exons'])
        cds_json = json.dumps(data['cds'])
        tx_rows.append((
            tid, data['gene_id'], data['chrom'], data['strand'],
            exons_json, cds_json,
        ))

    # Bulk insert
    if gene_rows:
        append_rows(con, 'gtf_genes',
                    ['gene_id', 'gene_name', 'chrom', 'strand', '"start"', '"end"', 'biotype'],
                    gene_rows)
    if tx_rows:
        append_rows(con, 'gtf_transcripts',
                    ['transcript_id', 'gene_id', 'chrom', 'strand', 'exons', 'cds'],
                    tx_rows)

    logger.info("GTF loaded: %d genes, %d transcripts", len(gene_rows), len(tx_rows))
    return {'genes': len(gene_rows), 'transcripts': len(tx_rows)}


def load_gene_names(gtf_file, con=None):
    """Load only gene_id → gene_name mapping into DuckDB gtf_genes table.

    Skips transcript structures (exons/CDS) entirely. Suitable when
    coordinate queries are handled by GTFIndex + IntervalIndex, and
    only gene-name backfill needs DuckDB.

    Args:
        gtf_file: path to GTF file
        con: DuckDB connection (uses thread-local default if None)

    Returns:
        int: number of genes loaded
    """
    if con is None:
        con = get_connection()

    gene_rows = []
    logger.info("Loading gene names from GTF: %s", gtf_file)
    with open(gtf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
            if parts[2] != 'gene':
                continue
            attrs = _parse_attributes(parts[8])
            gid = attrs.get('gene_id', 'NA')
            gname = attrs.get('gene_name', gid)
            biotype = attrs.get('gene_biotype', attrs.get('gene_type', ''))
            gene_rows.append((
                gid, gname, parts[0], parts[6],
                int(parts[3]), int(parts[4]), biotype,
            ))

    if gene_rows:
        append_rows(con, 'gtf_genes',
                    ['gene_id', 'gene_name', 'chrom', 'strand', '"start"', '"end"', 'biotype'],
                    gene_rows)

    logger.info("Gene names loaded: %d genes", len(gene_rows))
    return len(gene_rows)
