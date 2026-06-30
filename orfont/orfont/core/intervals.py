"""Fast interval overlap queries using binning-based indexing.

Provides a generic IntervalIndex for O(1)-average overlap queries against
large sets of genomic intervals (genes, CDS, exons). Uses the same binning
strategy as BedgraphIndex but extended with strand-aware and frame-aware
filtering for CDS overlap detection.

Memory: ~2 × number_of_intervals (one pointer per bin crossing).
Build:  O(n) single pass through sorted intervals.
Query:  O(1) amortised for sparse intervals, O(bins_per_query + k)
        where k is the number of true overlaps.
"""

from collections import defaultdict


class IntervalIndex:
    """Binning-based genomic interval index for fast overlap queries.

    Stores intervals in fixed-size bins (default 10 kb). Query by (chrom,
    start, end) returns all intervals overlapping the query region.

    Usage:
        idx = IntervalIndex()
        for chrom, start, end, data in intervals:
            idx.add(chrom, start, end, data)
        idx.build()   # sort and dedup within bins
        hits = idx.query('chr1', 100, 500)
    """

    def __init__(self, bin_size=10000):
        self.bin_size = bin_size
        self._bins = defaultdict(lambda: defaultdict(list))
        self._count = 0

    def add(self, chrom, start, end, data=None):
        """Add an interval (1-based, inclusive)."""
        start0 = start - 1
        end0 = end
        first_bin = start0 // self.bin_size
        last_bin = (end0 - 1) // self.bin_size
        entry = (start, end, data)
        for b in range(first_bin, last_bin + 1):
            self._bins[chrom][b].append(entry)
        self._count += 1

    def build(self):
        """Sort and deduplicate entries within each bin (call after all adds)."""
        for chrom_bins in self._bins.values():
            for b in chrom_bins:
                chrom_bins[b] = list(set(chrom_bins[b]))
                chrom_bins[b].sort(key=lambda x: x[0])

    def query(self, chrom, start, end):
        """Yield (interval_start, interval_end, data) overlapping (start, end).

        Coordinates are 1-based inclusive.
        """
        chrom_bins = self._bins.get(chrom)
        if not chrom_bins:
            return

        start0 = start - 1
        end0 = end
        first_bin = start0 // self.bin_size
        last_bin = (end0 - 1) // self.bin_size

        seen = set()
        for b in range(first_bin, last_bin + 1):
            for iv_start, iv_end, data in chrom_bins.get(b, []):
                iv_key = (iv_start, iv_end, id(data))
                if iv_key in seen:
                    continue
                seen.add(iv_key)
                if iv_start <= end and iv_end >= start:
                    yield iv_start, iv_end, data

    def has_overlap(self, chrom, start, end):
        """Return True if any interval overlaps the query region."""
        for _ in self.query(chrom, start, end):
            return True
        return False

    def count_overlaps(self, chrom, start, end):
        """Return count of overlapping intervals."""
        return sum(1 for _ in self.query(chrom, start, end))

    def __len__(self):
        return self._count


def build_cds_index(con, bin_size=10000):
    """Build an IntervalIndex of CDS regions from the gtf_transcripts table.

    Args:
        con: DuckDB connection
        bin_size: bin width in bp

    Returns:
        IntervalIndex keyed by chrom, with data = (strand, gene_id)
    """
    import json
    idx = IntervalIndex(bin_size=bin_size)

    rows = con.execute("""
        SELECT chrom, strand, cds, gene_id
        FROM gtf_transcripts
        WHERE cds IS NOT NULL AND cds != '[]'
    """).fetchall()

    for chrom, strand, cds_json, gene_id in rows:
        for cds_start, cds_end in json.loads(cds_json):
            idx.add(chrom, cds_start, cds_end, (strand, gene_id))

    idx.build()
    return idx


def query_cds_overlap_inframe(idx, chrom, strand, orf_start, orf_end, orf_frame):
    """Check if an ORF overlaps in-frame with any annotated CDS.

    Two intervals are in-frame when the CDS frame (start % 3 for + strand,
    end % 3 for - strand) matches the ORF frame.

    Args:
        idx: IntervalIndex from build_cds_index()
        chrom, strand: ORF chromosome and strand
        orf_start, orf_end: ORF genomic coordinates (1-based)
        orf_frame: ORF reading frame (0, 1, 2)

    Returns:
        bool: True if in-frame CDS overlap exists
    """
    for iv_start, iv_end, (cds_strand, gene_id) in idx.query(chrom, orf_start, orf_end):
        if cds_strand != strand:
            continue
        cds_frame = iv_start % 3 if strand == '+' else iv_end % 3
        if cds_frame == orf_frame:
            return True
    return False


def query_overlapping_genes(idx, chrom, strand, orf_start, orf_end):
    """Return gene_ids whose CDS intervals overlap the ORF on the same strand.

    Args:
        idx: IntervalIndex from build_cds_index()
        chrom, strand: ORF chromosome and strand
        orf_start, orf_end: ORF genomic coordinates (1-based)

    Returns:
        list of gene_id strings
    """
    genes = set()
    for iv_start, iv_end, (cds_strand, gene_id) in idx.query(chrom, orf_start, orf_end):
        if cds_strand == strand:
            genes.add(gene_id)
    return list(genes)


# ---------------------------------------------------------------------------
# Adapter: build IntervalIndex from GTFIndex (in-memory dict)
# ---------------------------------------------------------------------------

def build_cds_index_from_gtfindex(gtf_index, bin_size=10000):
    """Build an IntervalIndex from a GTFIndex instance's CDS data.

    Enables drop-in replacement of gtf_index.find_cds_overlap_inframe()
    with IntervalIndex-based O(1) queries.

    Args:
        gtf_index: GTFIndex instance with populated cds_by_chrom
        bin_size: bin width in bp

    Returns:
        IntervalIndex with data = (strand, gene_id)
    """
    idx = IntervalIndex(bin_size=bin_size)
    for chrom, intervals in gtf_index.cds_by_chrom.items():
        for iv_start, iv_end, iv_strand, iv_gid in intervals:
            idx.add(chrom, iv_start, iv_end, (iv_strand, iv_gid))
    idx.build()
    return idx


def build_gene_index_from_gtfindex(gtf_index, bin_size=10000):
    """Build an IntervalIndex from a GTFIndex instance's gene data.

    Enables O(1) gene-overlap queries replacing gtf_index.find_overlapping_genes().

    Args:
        gtf_index: GTFIndex instance with populated gene_by_chrom
        bin_size: bin width in bp

    Returns:
        IntervalIndex with data = (strand, gene_id)
    """
    idx = IntervalIndex(bin_size=bin_size)
    for chrom, intervals in gtf_index.gene_by_chrom.items():
        for iv_start, iv_end, iv_strand, iv_gid in intervals:
            idx.add(chrom, iv_start, iv_end, (iv_strand, iv_gid))
    idx.build()
    return idx


def annotate_cds_overlap_fast(candidates, cds_index, gene_index=None):
    """Annotate CDS overlap and overlapping genes using IntervalIndex.

    Drop-in replacement for annotate_cds_overlap() in unify_orf_predictions.py.
    Uses binning-based IntervalIndex instead of bisect scanning.

    Args:
        candidates: list of ORFCandidate objects
        cds_index: IntervalIndex from build_cds_index_from_gtfindex()
        gene_index: optional IntervalIndex from build_gene_index_from_gtfindex()
    """
    for cand in candidates:
        if gene_index:
            cand.overlapping_gene_ids = []
            for _, _, (gs, gid) in gene_index.query(
                cand.chrom, cand.start, cand.end
            ):
                if gs == cand.strand:
                    cand.overlapping_gene_ids.append(gid)

        cand.is_cds_overlap = False
        for iv_start, iv_end, (cds_strand, _) in cds_index.query(
            cand.chrom, cand.start, cand.end
        ):
            if cds_strand != cand.strand:
                continue
            cds_frame = iv_start % 3 if cand.strand == '+' else iv_end % 3
            if cds_frame == cand.frame:
                cand.is_cds_overlap = True
                break


# ---------------------------------------------------------------------------
# Chromosome sharding — split ORF tables by chrom for parallel processing
# ---------------------------------------------------------------------------

def shard_by_chromosome(con, source_table='raw_orfs', output_dir='.',
                         file_prefix='orfs', format='parquet'):
    """Export each chromosome's ORFs to a separate Parquet/CSV file.

    Enables parallel per-chromosome dedup and classification, reducing
    peak memory by processing one chromosome at a time.

    Args:
        con: DuckDB connection
        source_table: table to shard (raw_orfs or unified_orfs)
        output_dir: directory for output files
        file_prefix: prefix for output filenames
        format: 'parquet' or 'csv'

    Returns:
        list of (chrom, filepath) tuples
    """
    import os
    os.makedirs(output_dir, exist_ok=True)

    chroms = con.execute(
        f"SELECT DISTINCT chrom FROM {source_table} ORDER BY chrom"
    ).fetchall()

    results = []
    for (chrom,) in chroms:
        path = os.path.join(output_dir, f'{file_prefix}_{chrom}.{format}')
        if format == 'parquet':
            con.execute(f"""
                COPY (SELECT * FROM {source_table} WHERE chrom = '{chrom}')
                TO '{path}' (FORMAT parquet)
            """)
        else:
            con.execute(f"""
                COPY (SELECT * FROM {source_table} WHERE chrom = '{chrom}')
                TO '{path}' (FORMAT CSV, HEADER, DELIMITER '\t')
            """)
        results.append((chrom, path))

    return results


def merge_shards(con, glob_pattern, target_table):
    """Load all shard files matching a glob into a table.

    Args:
        con: DuckDB connection
        glob_pattern: glob pattern for shard files
        target_table: table to insert into (must exist)
    """
    con.execute(f"""
        INSERT INTO {target_table}
        SELECT * FROM read_parquet('{glob_pattern}')
    """)
