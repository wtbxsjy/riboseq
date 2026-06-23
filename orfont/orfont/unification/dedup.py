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
# Frame-aware merge (Phase 2) — sweep-line clustering on single-exon ORFs
# ---------------------------------------------------------------------------

def dedup_frame_aware(con, min_overlap_fraction=0.9):
    """Merge single-exon ORFs that are frame-compatible with >= threshold overlap.

    Implements the same logic as merge_frame_compatible_orfs() from the
    original script:
      - Only single-exon ORFs participate (multi-exon pass through unchanged)
      - ORFs must share chrom, strand, and frame (start % 3)
      - Overlap >= min_overlap_fraction of the shorter ORF
      - Transitive closure: A overlaps B, B overlaps C → A, B, C merged
      - Longest ORF in each merged cluster becomes representative

    Operates directly on the unified_orfs DuckDB table.  Multi-exon ORFs
    are untouched; merged single-exon ORFs have their tools/samples/n_tools/
    n_samples merged into the representative, and the shorter ORFs are deleted.
    """
    logger.info("Frame-aware merge: fetching single-exon ORFs from unified_orfs...")

    rows = con.execute("""
        SELECT orf_id, chrom, strand, genomic_start, genomic_end, length_nt,
               tools, samples, n_tools, n_samples
        FROM unified_orfs
        WHERE json_array_length(blocks::json) = 1
        ORDER BY chrom, strand, genomic_start % 3, genomic_start
    """).fetchall()

    if not rows:
        logger.info("Frame-aware merge: no single-exon ORFs to merge")
        total = con.execute("SELECT COUNT(*) FROM unified_orfs").fetchone()[0]
        return total

    # Sweep-line clustering within (chrom, strand, frame) groups
    used = set()
    merges = []       # list of (representative_orf_id, [merged_orf_ids])

    for i, row_i in enumerate(rows):
        if i in used:
            continue
        orf_id, chrom, strand, start_i, end_i, len_i, tools_i, samples_i, n_tools_i, n_samples_i = row_i
        frame_i = start_i % 3

        cluster = [(orf_id, len_i, start_i, end_i, tools_i, samples_i)]
        max_end = end_i
        used.add(i)

        for j in range(i + 1, len(rows)):
            if j in used:
                continue
            row_j = rows[j]
            chrom_j, strand_j, start_j, end_j = row_j[1], row_j[2], row_j[3], row_j[4]
            frame_j = start_j % 3

            # Different group → stop
            if chrom_j != chrom or strand_j != strand or frame_j != frame_i:
                break
            # Early exit: starts beyond furthest end of cluster
            if start_j > max_end:
                break

            # Compute overlap fraction relative to shorter ORF
            overlap_start = max(start_i, start_j)
            overlap_end = min(max_end, end_j)

            if overlap_start >= overlap_end:
                continue

            shorter_len = min(end_i - start_i, end_j - start_j)
            overlap_bp = overlap_end - overlap_start
            if overlap_bp / shorter_len >= min_overlap_fraction:
                orf_id_j, len_j, tools_j, samples_j = row_j[0], row_j[4], row_j[6], row_j[7]
                cluster.append((orf_id_j, len_j, start_j, end_j, tools_j, samples_j))
                max_end = max(max_end, end_j)
                used.add(j)

        if len(cluster) > 1:
            # Find representative (longest)
            rep = max(cluster, key=lambda x: x[1])
            merged_ids = [x[0] for x in cluster if x[0] != rep[0]]
            merges.append((rep[0], merged_ids, rep[4], rep[5], cluster))

    logger.info(
        "Frame-aware merge: %d single-exon ORFs → %d merge groups saving %d ORFs",
        len(rows), len(merges), sum(len(m[1]) for m in merges),
    )

    if not merges:
        total = con.execute("SELECT COUNT(*) FROM unified_orfs").fetchone()[0]
        return total

    # Build updates: merge tools/samples into representatives
    update_values = []
    for rep_id, merged_ids, _, _, cluster in merges:
        all_tools = set()
        all_samples = set()
        for _, _, _, _, tools_str, samples_str in cluster:
            for t in tools_str.split(','):
                if t.strip():
                    all_tools.add(t.strip())
            for s in samples_str.split(','):
                if s.strip():
                    all_samples.add(s.strip())
        merged_tools = ','.join(sorted(all_tools))
        merged_samples = ','.join(sorted(all_samples))
        update_values.append((rep_id, merged_tools, merged_samples,
                              len(all_tools), len(all_samples)))

    # Apply updates via temp table for efficiency
    import pandas as pd
    df_updates = pd.DataFrame(
        update_values,
        columns=['orf_id', 'tools', 'samples', 'n_tools', 'n_samples'],
    )
    con.register('_frame_merge_updates', df_updates)
    con.execute("""
        UPDATE unified_orfs u
        SET
            tools = fm.tools,
            samples = fm.samples,
            n_tools = fm.n_tools,
            n_samples = fm.n_samples
        FROM _frame_merge_updates fm
        WHERE u.orf_id = fm.orf_id
    """)
    con.unregister('_frame_merge_updates')

    # Delete merged-away ORFs
    all_merged_ids = []
    for _, merged_ids, _, _, _ in merges:
        all_merged_ids.extend(merged_ids)

    if all_merged_ids:
        placeholders = ','.join(['?'] * len(all_merged_ids))
        con.execute(
            f"DELETE FROM unified_orfs WHERE orf_id IN ({placeholders})",
            all_merged_ids,
        )

    total = con.execute("SELECT COUNT(*) FROM unified_orfs").fetchone()[0]
    logger.info(
        "Frame-aware merge complete: %d ORFs after merge (removed %d)",
        total, len(all_merged_ids),
    )
    return total


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
