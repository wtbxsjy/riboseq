# RiboCode \_predict\_psite Parameter Sweep — Lishuqi Ribo\_13

## Date: 2026-07-13

## Background
RiboCode failed on Ribo\_12 and Ribo\_13 (0 ORFs) despite 100% read retention in transcriptome sORF filter (`sorf_unique_mode_transcriptome='off'`). The error: "No obviously periodicity are detected from bam file". Root cause: hardcoded `f0.sum() >= 10` threshold in `\_predict\_psite()` (`metaplots.py:166`).

## Parameter Sweep Results (Ribo\_13, 2.11M reads)

| f0\_percent | -m/-M | Result | Notes |
|-----------|-------|--------|-------|
| 0.2 | 20-40 | FAIL | Widest possible range, lowest f0 |
| 0.3 | 20-40 | FAIL | |
| 0.3 | 24-40 | FAIL | |
| 0.4 | 20-40 | FAIL | |
| 0.4 | 24-40 | FAIL | |
| 0.5 | 20-40 | FAIL | |
| 0.5 | 24-40 | FAIL | |
| 0.65 | 24-36 | FAIL | Default (pipeline standard) |

## CLI-Adjustable Thresholds (all tested, none effective)

| Parameter | Default | Tested Range | Effect |
|-----------|---------|-------------|--------|
| `-f0_percent` | 0.65 | 0.2–0.5 | None — bottleneck is the 10-read minimum, not the proportion |
| `-m/-M` (read length) | 24/36 | 20–40 | None — no read length has ≥10 frame0 reads at start codon |

## Hardcoded Thresholds (require source patch)

| Threshold | Location | Value |
|-----------|----------|-------|
| `f0.sum() >= 10` | `metaplots.py:166` | 10 |
| `d.sum() >= 10` | `metaplots.py:200` | 10 |
| `pvalue1/2_cutoff` | `RiboCode_onestep.py:48` | 0.001 |
| `PSITE_SUM_CUTOFF` | `detectORF.py:306` | 5 |

## Bypass Option: P-site Config File

`RiboCode_onestep` → `meta_analysis()` → `\_predict\_psite()` → `_pre_config.txt` → `detectORF`. A manually-created config (using riboWaltz P-site offsets: all lengths = 12 nt) could bypass `\_predict\_psite()` entirely via the `-c` flag to `RiboCode` (non-onestep mode).

## Final RiboCode Results (Lishuqi, mode='off')

| Sample | Reads | ORFs |
|--------|-------|------|
| Ribo_11 | 1.09M | 4,529 |
| Ribo_12 | 2.18M | 0 (periodicity fail) |
| Ribo_13 | 2.11M | 0 (periodicity fail) |
| Ribo_14 | 2.40M | 11,419 |
| Ribo_15 | 2.37M | 16,907 |
| Ribo_16 | 2.25M | 6,306 |

Success rate: 4/6 (67%). Total ORFs: 39,161.

## Conclusion
Ribo_12 and Ribo_13 genuinely lack sufficient periodic ribosome footprint signal at start codons. CLI parameters cannot resolve this — the hardcoded `f0.sum() >= 10` threshold must be lowered or bypassed with external P-site config. Neither the read-length range nor the f0_percent ratio is the limiting factor.
