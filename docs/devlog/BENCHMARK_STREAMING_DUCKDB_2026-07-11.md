# Streaming DuckDB Unification Benchmark

> **Date**: 2026-07-11 | **Branch**: `dev` | **Dataset**: GSE157490 8 RPF samples

## Test Configuration

```
8 RPF samples × 5 tools (Ribo-TISH, Ribotricer, RiboCode, ORFquant, PRICE)
--duckdb-db unify_bench.db
--duckdb-memory-limit 32GB
--threads 2
--frame-merge-min-overlap 0.9
```

## Results

### Timing Breakdown

| Phase | Time (s) | % |
|-------|----------|---|
| GTF loading + gene names | 60 | 4.9 |
| Parse 12.4M raw ORFs (40 files) | 690.5 | 56.2 |
| SQL exact dedup (12.4M → 943K) | 29.0 | 2.4 |
| SQL frame-merge (943K → 794K) | 11.0 | 0.9 |
| Gene name backfill | 0.5 | 0.0 |
| Sequence extraction (794K ORFs) | 239 | 19.5 |
| AA translation (Python UDF) | 157 | 12.8 |
| CDS overlap (IntervalIndex) | 10.3 | 0.8 |
| Write BED/GTF/metadata | 89 | 7.2 |
| **Total** | **1,228.6 (20.5 min)** | 100 |

### Resource Usage

| Metric | Value |
|--------|-------|
| Python peak RSS | ~6 GB |
| DuckDB file (disk) | 5.0 GB |
| BED output | 143 MB |
| GTF output | 1.8 GB |
| Metadata output | 2.4 GB |

### Tool Distribution

| Tool | ORFs |
|------|------|
| Ribo-TISH | 684,667 |
| ORFquant | 72,490 |
| Ribotricer | 29,124 |
| PRICE | 23,810 |
| RiboCode | 616 |
| **Total** | **793,896** |

### Comparison: Old vs New

| Metric | Old (batch, original script) | New (streaming DuckDB) | Improvement |
|--------|---------------------------|----------------------|-------------|
| Total time | ~5h | 20.5 min | **15x** |
| Parse phase | ~4h | 11.5 min | 21x |
| Dedup + merge | ~30 min | 40s | 45x |
| Python peak RSS | ~50 GB | ~6 GB | **8x** |
| All 5 tools | ✅ | ✅ | — |
| Single run (no batches) | ❌ (OOM) | ✅ | — |

### Extrapolation: 32 RPF Samples

| Metric | Estimate |
|--------|----------|
| Parse 160 files | ~46 min |
| Total time | ~80 min (1.3h) |
| Python peak RSS | ~10 GB |
| DB file | ~20 GB |
| Unified ORFs | ~3.1M |

## Key Changes

1. `orfont/core/db.py`: `streaming_appender` — batch pandas register+INSERT (500K rows/batch)
2. `orfont/unification/builder.py`: Replace `all_rows` list accumulation; add `db_file`/`memory_limit` params
3. `scripts/run_orf.py`: `--duckdb-db` / `--duckdb-memory-limit` CLI args
4. `orfont/unification/parsers.py`: Fix RiboCode `.gtf.gz` detection + PRICE `chr+strand` format
5. `scripts/unify_orf_predictions.py`: Skip CDS annotation in sequential path (OOM workaround)
