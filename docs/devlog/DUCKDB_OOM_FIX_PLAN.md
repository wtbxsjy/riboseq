# DuckDB Streaming + Persistent DB — OOM 消除方案

> **Date**: 2026-07-11 | **Status**: Implemented & tested (1 sample) | **Benchmark**: below

## 1. Root Cause

orfont 路径 OOM 根因是**三重内存累积**:

```
all_rows = []           ← ① Python list 累积 42M dicts (~12.6 GB)
append_rows(con, ...)   ← ② Pandas DataFrame 再复制 (~10.5 GB)
duckdb(':memory:')      ← ③ 全内存 DB + 48GB 上限 → 超限 SIGKILL
                        ────────────────
                        ≈ 82 GB → 冲到 174 GB → OOM
```

## 2. Solution Implemented

### 2.1 Persistent DuckDB (`orfont/core/db.py`)

```python
configure(db_file='unify.db', memory_limit='32GB', threads=2)
```

- `db_file`: 数据写磁盘，memory_limit 仅控制缓存
- `memory_limit`: 通过 `--duckdb-memory-limit` 暴露为 CLI 参数

### 2.2 Streaming Batch Insert (`orfont/core/db.py`)

```python
class streaming_appender:
    """Batched pandas register + INSERT — 500K rows per batch."""
    def _append(self, row):
        self._batch.append(row)      # accumulate in Python (500K max)
        if len(batch) >= 500000:
            self._flush()            # register → INSERT → unregister → clear
```

替换了原来的 `all_rows = []` 全量累积 + 单次 bulk insert。

### 2.3 解析循环 (`orfont/unification/builder.py`)

```python
with streaming_appender(con, 'raw_orfs', columns=RAW_ORF_COLUMNS, batch_size=500000) as append:
    for file_path in files:
        candidates = parser(file_path, ...)
        for c in candidates:
            append(orf_to_row(c, tool, sample))
        candidates.clear()     # 显式释放
        del candidates
```

### 2.4 CLI 参数 (`run_orf.py`)

```
--duckdb-db PATH           # 持久化数据库文件 (default: in-memory)
--duckdb-memory-limit 32GB # DuckDB 缓存上限
```

## 3. Memory Budget

| 组件 | 修改前 | 修改后 |
|------|--------|--------|
| Python all_rows list | 12.6 GB | ~150 MB (500K batch) |
| Pandas DataFrame | 10.5 GB | ~150 MB (500K batch) |
| DuckDB in-memory | 8.4 GB | 0 (on disk) |
| DuckDB buffer | 48 GB crash | 32 GB (safe with 188 GB RAM) |
| **Python peak** | **~174 GB → OOM** | **~5 GB** |
| **Disk** | 0 | ~3 GB per sample (.db file) |

## 4. Benchmark Results

### 4.1 Test: 1 sample × 5 tools (SRR12588864)

```bash
python3 run_orf.py unify \
    --duckdb-db unify_test.db --duckdb-memory-limit 32GB \
    --threads 2 --frame-merge-min-overlap 0.9 \
    --ribotish ... --ribotricer ... --ribocode ... --orfquant ... --price ...
```

| Metric | Value |
|--------|-------|
| **Total time** | 546.7s (9.1 min) |
| **GTF loading** | ~55s |
| **Parsing + insert** | ~240s (5 tools) |
| **Dedup** | ~60s |
| **Frame-merge** | ~30s |
| **Sequence extraction** | ~120s |
| **AA translation + CDS overlap** | ~40s |
| **DB file size** | 3.1 GB |
| **Python peak RSS** | ~5 GB |
| **Unified ORFs** | 597,301 |

### 4.2 Tool Contribution

| Tool | ORFs |
|------|------|
| Ribo-TISH | 558,972 |
| Ribotricer | 19,482 |
| ORFquant | 19,394 |
| PRICE | 5,290 |
| RiboCode | 401 |

### 4.3 Extrapolation: 32 RPF samples

| Metric | Estimate |
|--------|----------|
| Total time | ~4.8h (32 × 9.1 min, linear scaling) |
| DB file size | ~100 GB |
| Python peak RSS | ~8 GB |
| Unified ORFs | ~3.1M (before cross-batch dedup) |

## 5. Comparison with Batch Approach

| Dimension | Batch (4 × 8) | Streaming DuckDB |
|-----------|--------------|-----------------|
| Runs needed | 4 | **1** |
| Manual steps | Launch 4 + merge | **1 command** |
| Memory peak | 50 GB/batch | **8 GB** |
| Total time | 4 × 5h = 20h | **~4.8h** |
| Cross-batch dedup | Extra step | **Built-in** |
| Disk output | 4 × 2.5 GB = 10 GB | ~100 GB (.db) + 10 GB outputs |

## 6. Files Changed

| File | Changes |
|------|---------|
| `orfont/core/db.py` | `streaming_appender` class (batch pandas register+INSERT) |
| `orfont/unification/builder.py` | Replace `all_rows` list with streaming appender; add `db_file`/`memory_limit` params; explicit cleanup |
| `scripts/run_orf.py` | Add `--duckdb-db` and `--duckdb-memory-limit` CLI args |

## 7. Usage

```bash
# Full 32 RPF samples, single run, no OOM
python3 run_orf.py unify \
    --gtf human_SARS2.genome.filtered.gtf \
    --fasta human_SARS2.genome.fa \
    --output unified_orfs \
    --min_len 6 --threads 2 \
    --frame-merge-min-overlap 0.9 \
    --duckdb-db unify.db \
    --duckdb-memory-limit 32GB \
    --ribotish *_pred.txt \
    --ribotricer *_translating_ORFs.tsv \
    --ribocode *.gtf.gz \
    --orfquant *_Detected_ORFs.gtf.gz \
    --price *.orfs.tsv
```

## 8. Pending

- [ ] Full 32 RPF sample benchmark (estimated 4.8h)
- [ ] Verify no cross-sample ORF duplication issue
- [ ] GENCODE + ORF-type classification on final output
- [ ] Compare results with batch_1-4 merged output
