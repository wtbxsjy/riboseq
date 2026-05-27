"""I/O module for streaming data into and out of DuckDB."""
from orfont.io.gtf import load_gtf, load_gene_names
from orfont.io.orfs import (
    orf_to_row, insert_raw_orfs, candidates_to_rows,
    export_raw_to_parquet, load_parquet_to_raw, export_unified_to_parquet,
    RAW_ORF_COLUMNS,
)
