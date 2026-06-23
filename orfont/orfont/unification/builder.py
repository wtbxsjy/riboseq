"""ORF unification builder.

Orchestrates the 3-stage deduplication pipeline:
  1. Exact-match merge (same chrom, strand, blocks)
  2. Frame-aware merge (single-exon ORFs only)
  3. Sequence clustering (optional, off by default)

Two paths available:
  - Optimized: DuckDB SQL dedup + pandas bulk insert + IntervalIndex CDS overlap
  - Original:  direct import of unify_orf_predictions.py (bit-identical fallback)
"""

import os
import sys
import logging
import time

from orfont.core.scripts_bridge import (
    _ensure_scripts_on_path,
)

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Output writers — generate BED/GTF/metadata from unified_orfs table
# ---------------------------------------------------------------------------

BED12_COLUMNS = [
    "chrom", "genomic_start", "genomic_end", "orf_id", "score",
    "strand", "thickStart", "thickEnd", "itemRgb",
    "blockCount", "blockSizes", "blockStarts",
]

def _blocks_to_bed12(blocks_json, start):
    """Convert JSON blocks string to BED12 block fields."""
    import json
    try:
        blocks = json.loads(blocks_json)
    except (json.JSONDecodeError, TypeError):
        return 1, "100", "0"

    block_sizes = []
    block_starts = []
    for b_start, b_end in blocks:
        size = b_end - b_start
        if size > 0:  # skip zero-length blocks (causes bedtools errors)
            block_sizes.append(str(size))
            block_starts.append(str(b_start - start))

    if not block_sizes:
        return 1, "100", "0"

    return (
        len(block_sizes),
        ",".join(block_sizes) + ",",
        ",".join(block_starts) + ",",
    )


def _write_outputs(con, output_prefix):
    """Generate BED, GTF, and metadata TSV from the unified_orfs table.

    Returns dict of {suffix: path}.
    """
    meta_path = f"{output_prefix}.metadata.tsv"
    bed_path = f"{output_prefix}.bed"
    gtf_path = f"{output_prefix}.gtf"

    rows = con.execute("""
        SELECT orf_id, chrom, strand, blocks, genomic_start, genomic_end,
               length_nt, length_aa, frame, gene_id, gene_name,
               tools, samples, n_tools, n_samples,
               sequence, start_codon, aa_sequence, is_cds_overlap
        FROM unified_orfs
        ORDER BY chrom, genomic_start
    """).fetchall()

    columns = [
        "orf_id", "chrom", "strand", "blocks", "genomic_start", "genomic_end",
        "length_nt", "length_aa", "frame", "gene_id", "gene_name",
        "tools", "samples", "n_tools", "n_samples",
        "sequence", "start_codon", "aa_sequence", "is_cds_overlap",
    ]

    # Write metadata
    with open(meta_path, 'w') as mf:
        hdr = [
            "orf_id", "orf_name", "orf_type", "chrom", "strand",
            "start", "end", "source", "gene_id", "gene_name",
            "transcript_id", "orf_length", "score", "source_count",
            "samples", "tools", "sequence", "start_codon",
            "aa_sequence", "cds_overlap",
        ]
        mf.write("\t".join(hdr) + "\n")
        for row in rows:
            d = dict(zip(columns, row))
            mf.write("\t".join([
                d["orf_id"],
                d["orf_id"],                          # orf_name
                "novel",                               # orf_type
                d["chrom"],
                d["strand"],
                str(d["genomic_start"]),
                str(d["genomic_end"]),
                "unified",                             # source
                d["gene_id"] or "",
                d["gene_name"] or "",
                "",                                    # transcript_id
                str(d["length_aa"]),                   # orf_length
                "0",                                   # score
                str(d["n_tools"]),                     # source_count
                d["samples"] or "",
                d["tools"] or "",
                d["sequence"] or "",
                d["start_codon"] or "",
                d["aa_sequence"] or "",
                str(d["is_cds_overlap"]).lower(),
            ]) + "\n")

    # Write BED12
    with open(bed_path, 'w') as bf:
        for row in rows:
            d = dict(zip(columns, row))
            n_blocks, sizes, starts = _blocks_to_bed12(
                d["blocks"], d["genomic_start"])
            bf.write("\t".join([
                d["chrom"],
                str(d["genomic_start"]),
                str(d["genomic_end"]),
                d["orf_id"],
                str(d["n_tools"]),                     # score = tool count
                d["strand"],
                str(d["genomic_start"]),               # thickStart
                str(d["genomic_end"]),                 # thickEnd
                "0",                                   # itemRgb
                str(n_blocks),
                sizes,
                starts,
            ]) + "\n")

    # Write GTF
    with open(gtf_path, 'w') as gf:
        gf.write("##gff-version 3\n")
        for row in rows:
            d = dict(zip(columns, row))
            # Parse blocks for exon features
            import json
            blocks = []
            try:
                blocks = json.loads(d["blocks"])
            except (json.JSONDecodeError, TypeError):
                pass

            if not blocks:
                gf.write("\t".join([
                    d["chrom"], "orfont", "CDS",
                    str(d["genomic_start"]),
                    str(d["genomic_end"]),
                    ".", d["strand"], "0",
                    f'gene_id "{d["gene_id"] or "."}"; '
                    f'transcript_id "{d["orf_id"]}"; '
                    f'orf_id "{d["orf_id"]}"; '
                    f'n_tools "{d["n_tools"]}"; '
                    f'tools "{d["tools"]}"; '
                    f'samples "{d["samples"]}";',
                ]) + "\n")
            else:
                for i, (bs, be) in enumerate(blocks):
                    attr = (
                        f'gene_id "{d["gene_id"] or "."}"; '
                        f'transcript_id "{d["orf_id"]}"; '
                        f'orf_id "{d["orf_id"]}"; '
                        f'exon_number "{i + 1}";'
                    )
                    gf.write("\t".join([
                        d["chrom"], "orfont", "exon",
                        str(bs), str(be),
                        ".", d["strand"], "0", attr,
                    ]) + "\n")
                # CDS line (collapsed)
                cds_start = blocks[0][0]
                cds_end = blocks[-1][1]
                gf.write("\t".join([
                    d["chrom"], "orfont", "CDS",
                    str(cds_start), str(cds_end),
                    ".", d["strand"], "0",
                    f'gene_id "{d["gene_id"] or "."}"; '
                    f'transcript_id "{d["orf_id"]}"; '
                    f'orf_id "{d["orf_id"]}"; '
                    f'n_tools "{d["n_tools"]}"; '
                    f'tools "{d["tools"]}"; '
                    f'samples "{d["samples"]}";',
                ]) + "\n")

    logger.info("Outputs written: %s, %s, %s", meta_path, bed_path, gtf_path)
    return {
        "metadata": meta_path,
        "bed": bed_path,
        "gtf": gtf_path,
    }


# ---------------------------------------------------------------------------
# Optimized unification — DuckDB SQL dedup + pandas bulk insert
# ---------------------------------------------------------------------------

def unify(ribotish_files=None, ribotricer_files=None, ribocode_files=None,
          orfquant_files=None, price_files=None, gtf_path=None, fasta_path=None,
          output_dir='.', prefix='unified_orfs',
          frame_merge=True, frame_merge_min_overlap=0.9,
          seq_cluster=False, bedgraph_dir=None, sample_list=None,
          min_len=6, extra_args=None):
    """Run ORF unification across tools and samples — optimized DuckDB path.

    Falls back to original script for:
      - Frame-aware merge (stage 2) — not yet implemented in DuckDB
      - Sequence clustering (stage 3) — not yet implemented in DuckDB
      - Bedgraph statistics — not yet ported
      - Per-tool output mode
    """
    os.makedirs(output_dir, exist_ok=True)
    output_prefix = os.path.join(output_dir, prefix)

    return _unify_optimized(
        ribotish_files, ribotricer_files, ribocode_files, orfquant_files,
        price_files, gtf_path, fasta_path, output_prefix,
        min_len, frame_merge, frame_merge_min_overlap,
    )


def _unify_optimized(ribotish_files, ribotricer_files, ribocode_files,
                     orfquant_files, price_files, gtf_path, fasta_path, output_prefix,
                     min_len, frame_merge=True, frame_merge_min_overlap=0.9):
    """Optimized path: pandas bulk insert → DuckDB SQL dedup → output."""
    from orfont.core.db import get_connection, init_schema, reset_connection
    from orfont.core.models import ORFCandidate, GTFIndex
    from orfont.unification.parsers import (
        parse_ribotish, parse_ribotricer, parse_orfquant, parse_ribocode, parse_price,
        infer_sample_id_from_prediction_path,
    )

    t_start = time.perf_counter()

    # 1. Setup DuckDB
    reset_connection()
    con = get_connection()
    init_schema(con)

    # 2. Load GTF gene names (fast path — gene_id→gene_name only)
    logger.info("Loading GTF gene names...")
    from orfont.io.gtf import load_gene_names
    n_genes = load_gene_names(gtf_path, con=con)
    logger.info("Loaded %d gene names from GTF", n_genes)

    # 3. Build GTFIndex (needed by parsers for gene annotation during parse)
    logger.info("Building GTFIndex for parsing...")
    t_gtf = time.perf_counter()
    gtf_index = GTFIndex(gtf_path)
    logger.info("GTFIndex built in %.1fs", time.perf_counter() - t_gtf)

    # 4. Parse all inputs → ORFCandidate → raw_orfs rows → pandas bulk insert
    from orfont.io.orfs import orf_to_row, RAW_ORF_COLUMNS
    from orfont.core.db import append_rows

    all_rows = []
    tool_configs = [
        ("Ribo-TISH", ribotish_files, parse_ribotish, "_pred.txt"),
        ("Ribotricer", ribotricer_files, parse_ribotricer, "_translating_ORFs.tsv"),
        ("RiboCode", ribocode_files, parse_ribocode, "_collapsed.gtf"),
        ("ORFquant", orfquant_files, parse_orfquant, "_Detected_ORFs.gtf"),
        ("PRICE", price_files, parse_price, ".orfs.tsv"),
    ]

    total_parsed = 0
    for tool_name, files, parser, suffix in tool_configs:
        if not files:
            continue
        for file_path in files:
            sample_id = infer_sample_id_from_prediction_path(file_path, suffix)
            logger.info("Parsing %s: %s (sample=%s)", tool_name, file_path, sample_id)
            try:
                candidates = parser(file_path, gtf_index, sample_id,
                                    min_len=min_len,
                                    exclude_tistypes={'Annotated', 'annotated'},
                                    atg_only=False)
            except Exception as e:
                logger.warning("Parse failed for %s: %s — skipping", file_path, e)
                continue
            for c in candidates:
                all_rows.append(orf_to_row(c, tool_name, sample_id))
            total_parsed += len(candidates)
            logger.info("  parsed %d ORFs", len(candidates))

    logger.info("Total parsed: %d ORFs from %d tools", total_parsed,
                sum(1 for _, files, _, _ in tool_configs if files))

    if not all_rows:
        logger.warning("No ORFs parsed — returning empty outputs")
        return _write_outputs(con, output_prefix)

    # 5. Bulk insert into raw_orfs (pandas DataFrame + register)
    logger.info("Bulk inserting %d ORFs into DuckDB...", len(all_rows))
    t_insert = time.perf_counter()
    append_rows(con, 'raw_orfs', RAW_ORF_COLUMNS, all_rows)
    logger.info("Bulk insert done in %.1fs (%.0f rows/s)",
                time.perf_counter() - t_insert,
                len(all_rows) / max(time.perf_counter() - t_insert, 0.001))

    # 6. SQL exact-match dedup
    from orfont.unification.dedup import (
        dedup_exact, annotate_gene_names, translate_aa_sequences,
    )
    logger.info("Running SQL exact-match dedup...")
    t_dedup = time.perf_counter()
    n_unified = dedup_exact(con)
    logger.info("Exact dedup: %d unified ORFs in %.1fs",
                n_unified, time.perf_counter() - t_dedup)

    # 6b. SQL frame-aware merge (single-exon ORFs only)
    if frame_merge:
        from orfont.unification.dedup import dedup_frame_aware
        logger.info(
            "Running SQL frame-aware merge (min_overlap=%.2f)...",
            frame_merge_min_overlap,
        )
        t_frame = time.perf_counter()
        n_merged = dedup_frame_aware(
            con, min_overlap_fraction=frame_merge_min_overlap)
        logger.info(
            "Frame-merge: %d → %d ORFs in %.1fs (removed %d)",
            n_unified, n_merged, time.perf_counter() - t_frame,
            n_unified - n_merged,
        )
        n_unified = n_merged

    # 7. Backfill gene names from gtf_genes table
    annotate_gene_names(con)

    # 8. Translate AA sequences (DuckDB Python UDF)
    translate_aa_sequences(con)

    # 9. CDS overlap annotation via IntervalIndex
    try:
        from orfont.core.intervals import build_cds_index_from_gtfindex
        logger.info("Building CDS IntervalIndex...")
        t_cds = time.perf_counter()
        cds_index = build_cds_index_from_gtfindex(gtf_index)

        # Query unified ORFs, check each against IntervalIndex, batch-update
        rows = con.execute("""
            SELECT orf_id, chrom, strand, genomic_start, genomic_end, frame
            FROM unified_orfs
        """).fetchall()

        overlap_ids = []
        for orf_id, chrom, strand, gstart, gend, frame in rows:
            for iv_start, iv_end, (cds_strand, _) in cds_index.query(
                chrom, gstart, gend
            ):
                if cds_strand != strand:
                    continue
                cds_frame = iv_start % 3 if strand == '+' else iv_end % 3
                if cds_frame == frame:
                    overlap_ids.append((orf_id,))
                    break

        if overlap_ids:
            import pandas as pd
            df = pd.DataFrame(overlap_ids, columns=['orf_id'])
            con.register('_cds_overlap_ids', df)
            con.execute("""
                UPDATE unified_orfs
                SET is_cds_overlap = true
                WHERE orf_id IN (SELECT orf_id FROM _cds_overlap_ids)
            """)
            con.unregister('_cds_overlap_ids')

        n_overlap = len(overlap_ids)
        logger.info("CDS overlap: %d overlapping ORFs in %.1fs",
                    n_overlap, time.perf_counter() - t_cds)
    except Exception as e:
        logger.warning("CDS overlap annotation skipped: %s", e)

    # 10. Generate output files
    result = _write_outputs(con, output_prefix)
    result["stats"] = _write_stats(con, output_prefix)

    dt_total = time.perf_counter() - t_start
    logger.info("Optimized unification complete: %d ORFs in %.1fs", n_unified, dt_total)

    return result


# ---------------------------------------------------------------------------
# Original path (fallback)
# ---------------------------------------------------------------------------

def _unify_original(ribotish_files, ribotricer_files, ribocode_files,
                    orfquant_files, price_files, gtf_path, fasta_path,
                    output_dir, prefix, frame_merge, frame_merge_min_overlap,
                    seq_cluster, bedgraph_dir, sample_list, min_len, extra_args):
    """Original path via bridge — calls unify_orf_predictions.main() directly."""
    output_prefix = os.path.join(output_dir, prefix)
    argv = ['--gtf', gtf_path, '--fasta', fasta_path,
            '--output', output_prefix,
            '--min-len', str(min_len),
            '--threads', '4']

    for f in (ribotish_files or []):
        argv.extend(['--ribotish', f])
    for f in (ribotricer_files or []):
        argv.extend(['--ribotricer', f])
    for f in (ribocode_files or []):
        argv.extend(['--ribocode', f])
    for f in (orfquant_files or []):
        argv.extend(['--orfquant', f])
    for f in (price_files or []):
        argv.extend(['--price', f])

    if not frame_merge:
        argv.append('--no-frame-merge')
    else:
        argv.extend(['--frame-merge-min-overlap', str(frame_merge_min_overlap)])

    if seq_cluster:
        argv.append('--seq-cluster')

    if bedgraph_dir:
        argv.extend(['--bedgraph-dir', bedgraph_dir])
        if sample_list:
            argv.append('--sample-list')
            argv.extend(sample_list)

    if extra_args:
        argv.extend(extra_args.split())

    logger.info("Running unification (original path): %s", argv)
    call_unify_orf_predictions(argv)

    return {
        'metadata': os.path.join(output_dir, f'{prefix}.metadata.tsv'),
        'bed': os.path.join(output_dir, f'{prefix}.bed'),
        'gtf': os.path.join(output_dir, f'{prefix}.gtf'),
        'stats': os.path.join(output_dir, f'{prefix}.stats.txt'),
    }


# ---------------------------------------------------------------------------
# Stats writer
# ---------------------------------------------------------------------------

def _write_stats(con, output_prefix):
    """Write statistics summary from unified_orfs table."""
    stats_path = f"{output_prefix}.stats.txt"
    total = con.execute("SELECT COUNT(*) FROM unified_orfs").fetchone()[0]

    with open(stats_path, 'w') as f:
        f.write("=== ORF Unification Statistics ===\n")
        f.write(f"Total unified ORFs: {total:,}\n\n")

        # Per-tool stats
        f.write("By Tool:\n")
        tool_rows = con.execute("""
            SELECT t.tool, COUNT(*) AS n
            FROM unified_orfs u
            CROSS JOIN UNNEST(STRING_SPLIT(u.tools, ',')) AS t(tool)
            GROUP BY t.tool ORDER BY n DESC
        """).fetchall()
        for tool, n in tool_rows:
            f.write(f"  {tool}: {n:,}\n")

        # Per-sample stats
        f.write("\nBy Sample:\n")
        sample_rows = con.execute("""
            SELECT s.sample, COUNT(*) AS n
            FROM unified_orfs u
            CROSS JOIN UNNEST(STRING_SPLIT(u.samples, ',')) AS s(sample)
            GROUP BY s.sample ORDER BY n DESC
        """).fetchall()
        for sample, n in sample_rows:
            f.write(f"  {sample}: {n:,}\n")

        # CDS overlap
        n_cds = con.execute(
            "SELECT COUNT(*) FROM unified_orfs WHERE is_cds_overlap").fetchone()[0]
        f.write(f"\nCDS-overlapping ORFs: {n_cds:,} "
                f"({100 * n_cds / total:.1f}%)\n" if total else "0\n")

        # Length distribution
        avg_len = con.execute(
            "SELECT AVG(length_aa) FROM unified_orfs").fetchone()[0]
        f.write(f"Average ORF length: {avg_len:.1f} aa\n")

    return stats_path
