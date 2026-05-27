"""SQL-based ORF deduplication using DuckDB.

Replaces Python dict-key-based exact-match merging with DuckDB SQL GROUP BY.
Handles Phase 1 dedup (exact match) and provides statistics queries.

All functions work on the DuckDB raw_orfs / unified_orfs tables — no
in-memory ORFCandidate lists required.
"""

import logging

from orfont.core.db import get_connection

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Exact-match dedup (Phase 1) — SQL GROUP BY
# ---------------------------------------------------------------------------

EXACT_DEDUP_SQL = """
    INSERT INTO unified_orfs
    SELECT
        -- Deterministic ORF ID: orf_{chrom}_{start}_{end}_{strand}_{hash8}
        format('orf_{}_{}_{}_{}_{}',
            chrom,
            MIN(genomic_start),
            MAX(genomic_end),
            strand,
            left(md5(blocks), 8)
        ) AS orf_id,
        chrom,
        strand,
        blocks,
        MIN(genomic_start) AS genomic_start,
        MAX(genomic_end) AS genomic_end,
        MIN(length_nt) AS length_nt,
        MIN(length_aa) AS length_aa,
        MIN(frame) AS frame,
        MIN_BY(gene_id, length_nt) AS gene_id,
        '' AS gene_name,                -- populated later via GTF join
        STRING_AGG(DISTINCT tool, ',') AS tools,
        STRING_AGG(DISTINCT sample, ',') AS samples,
        COUNT(DISTINCT tool) AS n_tools,
        COUNT(DISTINCT sample) AS n_samples,
        FIRST(sequence) AS sequence,
        FIRST(start_codon) AS start_codon,
        '' AS aa_sequence,              -- populated later via translation
        false AS is_cds_overlap,
        true AS is_representative,
        '' AS representative_orf_id
    FROM raw_orfs
    GROUP BY chrom, strand, blocks
"""


def dedup_exact(con=None):
    """Run exact-match deduplication via SQL GROUP BY.

    Reads from raw_orfs, groups by (chrom, strand, blocks), and inserts
    into unified_orfs with deterministic ORF IDs.

    Args:
        con: DuckDB connection (uses default if None)

    Returns:
        int: number of unified ORFs inserted
    """
    if con is None:
        con = get_connection()
    con.execute(EXACT_DEDUP_SQL)
    count = con.execute("SELECT COUNT(*) FROM unified_orfs").fetchone()[0]
    logger.info("Exact dedup: %d unified ORFs", count)
    return count


# ---------------------------------------------------------------------------
# Gene name annotation — backfill from GTF table
# ---------------------------------------------------------------------------

ANNOTATE_GENE_NAMES_SQL = """
    UPDATE unified_orfs u
    SET gene_name = COALESCE(g.gene_name, u.gene_id)
    FROM gtf_genes g
    WHERE u.gene_id = g.gene_id
"""


def annotate_gene_names(con=None):
    """Backfill gene_name from the gtf_genes table."""
    if con is None:
        con = get_connection()
    con.execute(ANNOTATE_GENE_NAMES_SQL)
    logger.info("Gene names annotated")


# ---------------------------------------------------------------------------
# AA sequence translation — DuckDB Python UDF
# ---------------------------------------------------------------------------

def translate_aa_sequences(con=None):
    """Translate nucleotide sequences to AA using a Python UDF.

    Creates a temporary DuckDB Python UDF for translation and updates
    the aa_sequence column in unified_orfs.
    """
    if con is None:
        con = get_connection()

    from orfont.core.utils import translate_sequence

    con.create_function(
        "translate_nt",
        lambda nt: translate_sequence(nt) if nt else None,
        ['VARCHAR'], 'VARCHAR',
        type='native',
    )
    con.execute("""
        UPDATE unified_orfs
        SET aa_sequence = translate_nt(sequence)
        WHERE sequence IS NOT NULL AND sequence != ''
    """)
    logger.info("AA sequences translated")


# ---------------------------------------------------------------------------
# Statistics queries
# ---------------------------------------------------------------------------

def _to_dicts(rows, columns):
    """Convert fetchall() rows to list of dicts."""
    return [dict(zip(columns, row)) for row in rows]


def per_tool_summary(con=None):
    """Return per-tool statistics as a list of dicts."""
    if con is None:
        con = get_connection()
    rows = con.execute("""
        SELECT
            t.tool,
            COUNT(*) AS n_orfs,
            COUNT(DISTINCT u.chrom) AS n_chroms,
            ROUND(AVG(u.length_aa), 1) AS avg_length_aa,
            SUM(CASE WHEN u.is_cds_overlap THEN 1 ELSE 0 END) AS n_cds_overlap
        FROM unified_orfs u
        CROSS JOIN UNNEST(STRING_SPLIT(u.tools, ',')) AS t(tool)
        GROUP BY t.tool
        ORDER BY n_orfs DESC
    """).fetchall()
    return _to_dicts(rows, ['tool', 'n_orfs', 'n_chroms', 'avg_length_aa', 'n_cds_overlap'])


def per_sample_summary(con=None):
    """Return per-sample 0/1 detection matrix."""
    if con is None:
        con = get_connection()
    rows = con.execute("""
        SELECT
            s.sample,
            COUNT(*) AS n_orfs,
            COUNT(DISTINCT u.chrom) AS n_chroms,
            ROUND(AVG(u.length_aa), 1) AS avg_length_aa,
            SUM(CASE WHEN u.is_cds_overlap THEN 1 ELSE 0 END) AS n_cds_overlap
        FROM unified_orfs u
        CROSS JOIN UNNEST(STRING_SPLIT(u.samples, ',')) AS s(sample)
        GROUP BY s.sample
        ORDER BY n_orfs DESC
    """).fetchall()
    return _to_dicts(rows, ['sample', 'n_orfs', 'n_chroms', 'avg_length_aa', 'n_cds_overlap'])


def gene_level_summary(con=None):
    """Return per-gene ORF count summary."""
    if con is None:
        con = get_connection()
    rows = con.execute("""
        SELECT
            gene_id,
            gene_name,
            COUNT(*) AS n_orfs,
            COUNT(DISTINCT chrom) AS n_chroms,
            ROUND(AVG(length_aa), 1) AS avg_length_aa,
            SUM(CASE WHEN is_cds_overlap THEN 1 ELSE 0 END) AS n_cds_overlap
        FROM unified_orfs
        GROUP BY gene_id, gene_name
        ORDER BY n_orfs DESC
    """).fetchall()
    return _to_dicts(rows, ['gene_id', 'gene_name', 'n_orfs', 'n_chroms', 'avg_length_aa', 'n_cds_overlap'])


def full_classification_summary(con=None):
    """Return unified ORFs joined with all 3 classifiers.

    Requires classification tables to be populated first.
    """
    if con is None:
        con = get_connection()
    rows = con.execute("""
        SELECT
            u.orf_id,
            u.chrom, u.strand,
            u.genomic_start, u.genomic_end,
            u.length_aa, u.gene_id, u.gene_name,
            u.tools, u.samples, u.n_tools, u.n_samples,
            g.orf_biotype AS gencode_biotype,
            q.ORF_category_Gen,
            q.ORF_category_Tx,
            q.ORF_category_Tx_compatible,
            t.orf_type_category
        FROM unified_orfs u
        LEFT JOIN classification_gencode g USING (orf_id)
        LEFT JOIN classification_orfquant q USING (orf_id)
        LEFT JOIN classification_orftype t USING (orf_id)
    """).fetchall()
    return _to_dicts(rows, [
        'orf_id', 'chrom', 'strand', 'genomic_start', 'genomic_end',
        'length_aa', 'gene_id', 'gene_name', 'tools', 'samples',
        'n_tools', 'n_samples', 'gencode_biotype',
        'ORF_category_Gen', 'ORF_category_Tx', 'ORF_category_Tx_compatible',
        'orf_type_category',
    ])
