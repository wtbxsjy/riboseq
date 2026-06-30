#!/usr/bin/env python3
"""
Unified ORF pipeline entry point — auto-selects optimized or fallback path.

Detects whether orfont is installed and routes to the fastest available
implementation.  Falls back to the original scripts when orfont is unavailable.

Usage:
    python run_orf.py unify --ribotish a.txt --gtf ref.gtf --fasta ref.fa -o .
    python run_orf.py classify-gencode --bed u.bed --metadata u.tsv --ensembl-dir ./Ens110
"""

import os
import sys
import subprocess
import logging

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
RIBO_ROOT = os.path.dirname(SCRIPT_DIR)             # riboseq/
ORFONT_ROOT = os.path.join(RIBO_ROOT, 'orfont')      # riboseq/orfont/ (project root)

logging.basicConfig(
    level=logging.INFO,
    format='[%(asctime)s] %(levelname)s: %(message)s',
    datefmt='%H:%M:%S')
logger = logging.getLogger('run_orf')


def _orfont_available():
    """Return True if orfont package can be imported."""
    try:
        if ORFONT_ROOT not in sys.path:
            sys.path.insert(0, ORFONT_ROOT)
        from orfont.core.scripts_bridge import _scripts_dir
        _scripts_dir()  # verify script path resolves
        return True
    except Exception as e:
        logger.debug("orfont unavailable: %s", e)
        return False


def _run_original(argv):
    """Fall back to running the original scripts as subprocess."""
    cmd = argv[0]
    args = argv[1:]

    # Commands that go through classify_orfs_wrapper need --mode
    mode_map = {
        'classify-gencode': 'gencode',
        'classify-orftype': 'orf_type',
        'classify-orfquant': 'orfquant',
    }

    script_map = {
        'unify': os.path.join(SCRIPT_DIR, 'unify_orf_predictions.py'),
        'classify-gencode': os.path.join(SCRIPT_DIR, 'classify_orfs_wrapper.py'),
        'classify-orftype': os.path.join(SCRIPT_DIR, 'classify_orfs_wrapper.py'),
        'classify-orfquant': os.path.join(SCRIPT_DIR, 'classify_orfs_wrapper.py'),
    }
    script = script_map.get(cmd)
    if not script:
        logger.error("Unknown command: %s", cmd)
        return 1

    mode = mode_map.get(cmd)
    if mode:
        full_cmd = [sys.executable, script, '--mode', mode] + args
    else:
        full_cmd = [sys.executable, script] + args

    logger.info("Running original: %s ...", script)
    return subprocess.run(full_cmd).returncode


def _run_orfont(argv):
    """Route through orfont — uses optimized DuckDB path for unify, bridge for others."""
    if ORFONT_ROOT not in sys.path:
        sys.path.insert(0, ORFONT_ROOT)

    cmd = argv[0]
    args = argv[1:]

    logger.info("Running optimized (orfont): %s", cmd)

    if cmd == 'unify':
        return _run_orfont_unify(args)
    elif cmd == 'classify-gencode':
        from orfont.core.scripts_bridge import call_classify_gencode
        call_classify_gencode(['--mode', 'gencode'] + args)
    elif cmd == 'classify-orftype':
        from orfont.core.scripts_bridge import call_classify_orftype
        call_classify_orftype(args)
    else:
        logger.error("Unknown command: %s", cmd)
        return 1
    return 0


def _run_orfont_unify(args):
    """Parse unify args and call builder.unify() directly (optimized DuckDB path).

    Falls back to the original script when features not yet ported to orfont
    are requested: frame-merge, seq-cluster, bedgraph stats, or per-tool output.
    """
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('--ribotish', nargs='+', default=[])
    ap.add_argument('--ribotricer', nargs='+', default=[])
    ap.add_argument('--ribocode', nargs='+', default=[])
    ap.add_argument('--orfquant', nargs='+', default=[])
    ap.add_argument('--price', nargs='+', default=[])
    ap.add_argument('--gtf', required=True)
    ap.add_argument('--fasta', required=True)
    ap.add_argument('--output', required=True)
    # Accept both --min_len (Nextflow) and --min-len (CLI)
    ap.add_argument('--min-len', type=int, default=None)
    ap.add_argument('--min_len', type=int, default=None, dest='min_len_alt')
    ap.add_argument('--no-frame-merge', action='store_true')
    ap.add_argument('--frame-merge-min-overlap', type=float, default=0.9)
    ap.add_argument('--seq-cluster', action='store_true')
    ap.add_argument('--bedgraph-dir', default=None)
    ap.add_argument('--sample-list', nargs='+', default=None)
    ap.add_argument('--threads', type=int, default=4)
    ap.add_argument('--per-tool-output', default=None)
    parsed, extra = ap.parse_known_args(args)

    # Resolve --min_len / --min-len (NF uses --min_len)
    min_len = parsed.min_len if parsed.min_len is not None else (parsed.min_len_alt or 6)

    from orfont.unification.builder import unify

    # Determine output_dir and prefix from --output
    output_path = parsed.output
    output_dir = os.path.dirname(output_path) or '.'
    prefix = os.path.basename(output_path)

    unify(
        ribotish_files=parsed.ribotish or None,
        ribotricer_files=parsed.ribotricer or None,
        ribocode_files=parsed.ribocode or None,
        orfquant_files=parsed.orfquant or None,
        price_files=parsed.price or None,
        gtf_path=parsed.gtf,
        fasta_path=parsed.fasta,
        output_dir=output_dir,
        prefix=prefix,
        frame_merge=not parsed.no_frame_merge,
        frame_merge_min_overlap=parsed.frame_merge_min_overlap,
        seq_cluster=parsed.seq_cluster,
        bedgraph_dir=parsed.bedgraph_dir,
        sample_list=parsed.sample_list,
        min_len=min_len,
    )
    return 0


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        return 1

    cmd = sys.argv[1]
    rest = sys.argv[2:]

    use_orfont = _orfont_available()
    logger.info("orfont: %s",
                "available (optimized)" if use_orfont else "not available (fallback)")

    if use_orfont:
        rc = _run_orfont([cmd] + rest)
    else:
        rc = _run_original([cmd] + rest)

    sys.exit(rc)


if __name__ == '__main__':
    main()
