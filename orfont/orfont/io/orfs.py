"""ORF I/O — streaming raw ORFs into DuckDB and Parquet intermediate storage.

Provides the bridge between parser output (lists of tuples) and the DuckDB
raw_orfs table, enabling memory-efficient bulk loading at scale.
"""

import json
import logging
import os

from orfont.core.db import get_connection, append_rows

logger = logging.getLogger(__name__)

RAW_ORF_COLUMNS = [
    'chrom', 'strand', 'genomic_start', 'genomic_end', 'blocks',
    'transcript_id', 'gene_id', 'tool', 'sample',
    'score', 'pvalue', 'sequence', 'start_codon',
    'length_nt', 'length_aa', 'frame', 'batch_id',
]


def orf_to_row(cand, tool, sample, batch_id=0):
    """Convert an ORFCandidate to a tuple for raw_orfs insertion.

    Args:
        cand: ORFCandidate object (from unify_orf_predictions)
        tool: tool name string (e.g. 'Ribo-TISH')
        sample: sample identifier
        batch_id: batch identifier for tracking

    Returns:
        tuple of 17 values matching RAW_ORF_COLUMNS
    """
    return (
        cand.chrom,
        cand.strand,
        cand.start,
        cand.end,
        json.dumps(cand.blocks),
        cand.tid,
        cand.gid,
        tool,
        sample,
        cand.score,
        getattr(cand, 'pvalue', None),
        getattr(cand, 'sequence', None) or '',
        getattr(cand, 'start_codon', None) or '',
        cand.length_nt,
        cand.length_aa,
        cand.frame,
        batch_id,
    )


def insert_raw_orfs(con, rows, columns=None):
    """Bulk-insert raw ORF rows into the raw_orfs table.

    Args:
        con: DuckDB connection
        rows: list of tuples (each must have 17 values)
        columns: column names (defaults to RAW_ORF_COLUMNS)
    """
    if not rows:
        return 0
    if columns is None:
        columns = RAW_ORF_COLUMNS
    append_rows(con, 'raw_orfs', columns, rows)
    return len(rows)


def candidates_to_rows(candidates, tool, sample, batch_id=0):
    """Convert a list of ORFCandidate objects to raw_orfs tuples.

    Returns:
        list of tuples ready for insert_raw_orfs
    """
    return [orf_to_row(c, tool, sample, batch_id) for c in candidates]


def export_raw_to_parquet(con, output_path):
    """Export the raw_orfs table to a Parquet file."""
    con.execute(f"COPY (SELECT * FROM raw_orfs) TO '{output_path}' (FORMAT parquet)")
    logger.info("Exported raw_orfs to %s", output_path)


def load_parquet_to_raw(con, glob_pattern):
    """Load Parquet files matching a glob into the raw_orfs table."""
    con.execute(f"""
        INSERT INTO raw_orfs
        SELECT * FROM read_parquet('{glob_pattern}')
    """)
    logger.info("Loaded Parquet files matching %s into raw_orfs", glob_pattern)


def export_unified_to_parquet(con, output_path):
    """Export the unified_orfs table to a Parquet file."""
    con.execute(
        f"COPY (SELECT * FROM unified_orfs) TO '{output_path}' (FORMAT parquet)"
    )
    logger.info("Exported unified_orfs to %s", output_path)
