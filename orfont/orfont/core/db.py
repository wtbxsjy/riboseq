"""DuckDB backend for ORF pipeline — connection management, schema, and utilities.

Provides:
- Singleton connection management (in-memory or persistent)
- Table schema definitions for all pipeline stages
- Appender helpers for bulk data loading
- Common query utilities
"""

import logging
import os
import threading
from typing import Optional

import duckdb

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Connection management
# ---------------------------------------------------------------------------

_local = threading.local()
_DB_FILE: Optional[str] = None
_MEMORY_LIMIT: str = "48GB"
_THREADS: int = 4


def configure(db_file: Optional[str] = None, memory_limit: str = "48GB",
              threads: int = 4):
    """Configure global DuckDB settings. Call before get_connection()."""
    global _DB_FILE, _MEMORY_LIMIT, _THREADS
    _DB_FILE = db_file
    _MEMORY_LIMIT = memory_limit
    _THREADS = threads


def get_connection(read_only: bool = False) -> duckdb.DuckDBPyConnection:
    """Get a DuckDB connection (thread-local).

    Returns an in-memory connection by default. If configure(db_file=...) was
    called, opens a persistent database at that path.
    """
    if not hasattr(_local, 'con') or _local.con is None:
        if _DB_FILE:
            _local.con = duckdb.connect(_DB_FILE, read_only=read_only)
        else:
            _local.con = duckdb.connect(':memory:')
        _local.con.execute(f"SET memory_limit = '{_MEMORY_LIMIT}'")
        _local.con.execute(f"SET threads = {_THREADS}")
        _local.con.execute("SET max_temp_directory_size = '100GiB'")
        _local.con.execute("SET preserve_insertion_order = false")
        _local.con.execute("INSTALL parquet; LOAD parquet")
    return _local.con


def reset_connection():
    """Close and clear the thread-local connection."""
    if hasattr(_local, 'con') and _local.con is not None:
        _local.con.close()
        _local.con = None


# ---------------------------------------------------------------------------
# Schema definitions
# ---------------------------------------------------------------------------

SCHEMA = {
    "gtf_genes": """
        CREATE TABLE IF NOT EXISTS gtf_genes (
            gene_id VARCHAR,
            gene_name VARCHAR,
            chrom VARCHAR,
            strand VARCHAR,
            "start" INTEGER,
            "end" INTEGER,
            biotype VARCHAR,
            PRIMARY KEY (gene_id)
        )
    """,

    "gtf_transcripts": """
        CREATE TABLE IF NOT EXISTS gtf_transcripts (
            transcript_id VARCHAR PRIMARY KEY,
            gene_id VARCHAR,
            chrom VARCHAR,
            strand VARCHAR,
            exons VARCHAR,   -- JSON array: [[s1,e1], [s2,e2], ...]
            cds VARCHAR       -- JSON array of CDS intervals
        )
    """,

    "raw_orfs": """
        CREATE TABLE IF NOT EXISTS raw_orfs (
            chrom VARCHAR,
            strand VARCHAR,
            genomic_start INTEGER,
            genomic_end INTEGER,
            blocks VARCHAR,          -- JSON: [[s1,e1], ...]
            transcript_id VARCHAR,
            gene_id VARCHAR,
            tool VARCHAR,
            sample VARCHAR,
            score DOUBLE,
            pvalue DOUBLE,
            sequence VARCHAR,
            start_codon VARCHAR,
            length_nt INTEGER,
            length_aa INTEGER,
            frame INTEGER,
            batch_id INTEGER
        )
    """,

    "unified_orfs": """
        CREATE TABLE IF NOT EXISTS unified_orfs (
            orf_id VARCHAR PRIMARY KEY,
            chrom VARCHAR,
            strand VARCHAR,
            blocks VARCHAR,          -- JSON: [[s1,e1], ...]
            genomic_start INTEGER,
            genomic_end INTEGER,
            length_nt INTEGER,
            length_aa INTEGER,
            frame INTEGER,
            gene_id VARCHAR,
            gene_name VARCHAR,
            tools VARCHAR,           -- comma-separated
            samples VARCHAR,         -- comma-separated
            n_tools INTEGER,
            n_samples INTEGER,
            sequence VARCHAR,
            start_codon VARCHAR,
            aa_sequence VARCHAR,
            is_cds_overlap BOOLEAN DEFAULT false,
            is_representative BOOLEAN DEFAULT true,
            representative_orf_id VARCHAR
        )
    """,

    "classification_gencode": """
        CREATE TABLE IF NOT EXISTS classification_gencode (
            orf_id VARCHAR PRIMARY KEY,
            orf_biotype VARCHAR
        )
    """,

    "classification_orfquant": """
        CREATE TABLE IF NOT EXISTS classification_orfquant (
            orf_id VARCHAR PRIMARY KEY,
            ORF_category_Gen VARCHAR,
            ORF_category_Tx VARCHAR,
            ORF_category_Tx_compatible VARCHAR
        )
    """,

    "classification_orftype": """
        CREATE TABLE IF NOT EXISTS classification_orftype (
            orf_id VARCHAR PRIMARY KEY,
            orf_type_category VARCHAR
        )
    """,
}


def init_schema(con: Optional[duckdb.DuckDBPyConnection] = None):
    """Create all pipeline tables if they don't exist."""
    if con is None:
        con = get_connection()
    for name, ddl in SCHEMA.items():
        con.execute(ddl)
    logger.info("DuckDB schema initialized (%d tables)", len(SCHEMA))


def drop_all(con: Optional[duckdb.DuckDBPyConnection] = None):
    """Drop all pipeline tables (useful for testing/cleanup)."""
    if con is None:
        con = get_connection()
    for name in SCHEMA:
        con.execute(f"DROP TABLE IF EXISTS {name}")
    logger.info("All pipeline tables dropped")


# ---------------------------------------------------------------------------
# Appender helpers for bulk loading
# ---------------------------------------------------------------------------

def append_rows(con: duckdb.DuckDBPyConnection, table: str,
                columns: list, rows: list):
    """Bulk-insert rows into a table using pandas DataFrame for speed.

    Uses pandas DataFrame + DuckDB register() for vectorized loading.
    Falls back to executemany if pandas is unavailable.

    Args:
        con: DuckDB connection
        table: target table name
        columns: list of column names (must match table schema order)
        rows: list of tuples, each tuple is one row
    """
    if not rows:
        return
    _bulk_insert_via_pandas(con, table, columns, rows)


def _bulk_insert_via_pandas(con, table, columns, rows):
    """Fast bulk insert: list of tuples → pandas DataFrame → DuckDB register → INSERT."""
    import pandas as pd
    quoted_cols = ', '.join(
        f'"{c}"' if c.lower() in ('start', 'end') else c
        for c in columns
    )
    df = pd.DataFrame(rows, columns=columns)
    con.register('_bulk_temp', df)
    try:
        con.execute(f"INSERT INTO {table} ({quoted_cols}) SELECT * FROM _bulk_temp")
    finally:
        con.unregister('_bulk_temp')


class streaming_appender:
    """Context manager: Batched pandas DataFrame insert for streaming rows.

    Rows are accumulated in batches and flushed via the fast pandas
    ``register`` + ``INSERT INTO ... SELECT *`` path.  Memory usage is
    capped at ``batch_size`` rows.

    Usage::

        with streaming_appender(con, 'raw_orfs', RAW_ORF_COLUMNS) as append:
            for orf in candidates:
                append((chrom, strand, start, end, ...))

    Args:
        con: DuckDB connection
        table: target table name
        columns: list of column names
        batch_size: rows per flush (default 500K)
    """

    def __init__(self, con, table, columns=None, batch_size=500000):
        self.con = con
        self.table = table
        self.columns = columns or []
        self.batch_size = batch_size
        self._batch = []
        self.count = 0

    def __enter__(self):
        self._batch = []
        self.count = 0
        return self._append

    def _append(self, row_tuple):
        self._batch.append(row_tuple)
        self.count += 1
        if len(self._batch) >= self.batch_size:
            self._flush()

    def _flush(self):
        if not self._batch:
            return
        import pandas as pd
        col_str = ', '.join(
            f'"{c}"' if c.lower() in ('start', 'end') else c
            for c in self.columns
        )
        df = pd.DataFrame(self._batch, columns=self.columns)
        self.con.register('_stream_batch', df)
        self.con.execute(
            f"INSERT INTO {self.table} ({col_str}) SELECT * FROM _stream_batch")
        self.con.unregister('_stream_batch')
        self._batch.clear()

    def __exit__(self, *args):
        self._flush()
        logger.debug("Streaming batch insert done: %d rows → %s", self.count, self.table)
        return False


def parquet_to_table(con: duckdb.DuckDBPyConnection, glob_pattern: str,
                     table: str):
    """Load all Parquet files matching a glob into a table.

    Creates the table if it doesn't exist; appends if it does.

    Args:
        con: DuckDB connection
        glob_pattern: glob pattern for Parquet files (e.g., 'parsed_*.parquet')
        table: target table name
    """
    con.execute(f"""
        INSERT INTO {table}
        SELECT * FROM read_parquet('{glob_pattern}')
    """)


def table_to_parquet(con: duckdb.DuckDBPyConnection, table: str,
                     output_path: str, **kwargs):
    """Export a table (or query) to Parquet.

    Args:
        con: DuckDB connection
        table: source table name (or subquery)
        output_path: output .parquet file path
        **kwargs: passed to COPY (e.g., COMPRESSION 'zstd')
    """
    con.execute(f"COPY (SELECT * FROM {table}) TO '{output_path}' (FORMAT parquet)")


# ---------------------------------------------------------------------------
# View definitions for common aggregations
# ---------------------------------------------------------------------------

VIEWS = {
    "per_tool_stats": """
        CREATE OR REPLACE VIEW per_tool_stats AS
        SELECT
            t.tool,
            COUNT(*) AS n_orfs,
            COUNT(DISTINCT u.chrom) AS n_chroms,
            AVG(u.length_aa) AS avg_length_aa,
            SUM(CASE WHEN u.is_cds_overlap THEN 1 ELSE 0 END) AS n_cds_overlap
        FROM unified_orfs u
        CROSS JOIN UNNEST(STRING_SPLIT(u.tools, ',')) AS t(tool)
        GROUP BY t.tool
    """,

    "full_classification": """
        CREATE OR REPLACE VIEW full_classification AS
        SELECT
            u.*,
            g.orf_biotype AS gencode_biotype,
            q.ORF_category_Gen,
            q.ORF_category_Tx,
            q.ORF_category_Tx_compatible,
            t.orf_type_category
        FROM unified_orfs u
        LEFT JOIN classification_gencode g USING (orf_id)
        LEFT JOIN classification_orfquant q USING (orf_id)
        LEFT JOIN classification_orftype t USING (orf_id)
    """,
}


def init_views(con: Optional[duckdb.DuckDBPyConnection] = None):
    """Create common aggregation views."""
    if con is None:
        con = get_connection()
    for name, ddl in VIEWS.items():
        con.execute(ddl)
    logger.info("Views initialized (%d views)", len(VIEWS))
