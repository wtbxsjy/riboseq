"""Bridge to external pipeline scripts — imports modules directly instead of subprocess.

All external scripts read sys.argv directly (argparse pattern). We temporarily
override sys.argv before calling main() to pass arguments without subprocess overhead.
"""

import os
import sys


def _scripts_dir():
    """Return the path to riboseq/scripts directory.

    Resolves relative to this file. Supports two layouts:
      - riboseq/orfont/orfont/core/  →  ../../../scripts  (integrated)
      - orf_ontology/orfont/core/    →  ../../../riboseq/scripts  (standalone)

    Tries integrated layout first, then standalone.
    """
    import os as _os
    this_dir = _os.path.dirname(_os.path.abspath(__file__))

    # Integrated: riboseq/orfont/orfont/core/ → ../../../scripts
    candidate = _os.path.abspath(_os.path.join(this_dir, '..', '..', '..', 'scripts'))
    if _os.path.isdir(candidate):
        return candidate

    # Standalone: orf_ontology/orfont/core/ → ../../../riboseq/scripts
    candidate = _os.path.abspath(_os.path.join(this_dir, '..', '..', '..', 'riboseq', 'scripts'))
    if _os.path.isdir(candidate):
        return candidate

    raise FileNotFoundError(
        f"Cannot find riboseq/scripts dir from {this_dir}. "
        f"Tried ../../../scripts and ../../../riboseq/scripts")


def _ensure_scripts_on_path():
    """Add scripts directory to sys.path if not already present."""
    scripts = _scripts_dir()
    if scripts not in sys.path:
        sys.path.insert(0, scripts)


def _call_main_with_argv(module_name, prog_name, argv):
    """Import module_name, set sys.argv, call module.main(), then restore sys.argv.

    All pipeline scripts use argparse with main() taking no arguments,
    reading from sys.argv directly.
    """
    import importlib
    _saved = sys.argv
    sys.argv = [prog_name] + argv
    try:
        mod = importlib.import_module(module_name)
        mod.main()
    finally:
        sys.argv = _saved


def call_unify_orf_predictions(argv):
    """Call unify_orf_predictions.py main() directly. Returns the output paths dict."""
    _ensure_scripts_on_path()
    _call_main_with_argv('unify_orf_predictions', 'unify_orf_predictions', argv)
    output_dir = None
    for i, arg in enumerate(argv):
        if arg == '--output' and i + 1 < len(argv):
            output_dir = argv[i + 1]
            break
    if output_dir is None:
        raise ValueError("--output not found in argv")
    return output_dir


def call_classify_gencode(argv):
    """Call classify_orfs_wrapper.py main() with gencode mode."""
    _ensure_scripts_on_path()
    _call_main_with_argv('classify_orfs_wrapper', 'classify_orfs_wrapper', argv)


def call_classify_orfquant(gtf_path, annotation_path, output_path, metadata_path=None,
                         output_prefix=None, cpus=1, parallel=True):
    """Run ORFquant classification via Rscript with mirai-based parallelization.

    Passes --parallel and --threads to the R script, which dispatches ORF
    chunks to mirai daemons.  Annotation is shared via temp RDS files so each
    daemon loads its own copy (no serialization overhead).  Falls back to
    serial orfquant_classify_orfs() when mirai is unavailable or n_cores <= 1.
    """
    import subprocess
    script = os.path.join(
        _scripts_dir(), 'class_orf', 'run_orfquant_classify.R')
    cmd = ['Rscript', script,
           '--input', gtf_path,
           '--annotation', annotation_path,
           '--output', output_path]
    if metadata_path:
        cmd.extend(['--metadata', metadata_path])
    if output_prefix:
        cmd.extend(['--output_prefix', output_prefix])
    if parallel and cpus > 1:
        cmd.extend(['--parallel', '--threads', str(cpus)])
    subprocess.run(cmd, check=True)


def call_classify_orftype(argv):
    """Call class_ORFtype.py main() directly."""
    _ensure_scripts_on_path()
    import importlib.util
    script_path = os.path.join(_scripts_dir(), 'class_orf', 'class_ORFtype.py')
    spec = importlib.util.spec_from_file_location('class_ORFtype', script_path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    _saved = sys.argv
    sys.argv = ['class_ORFtype'] + argv
    try:
        mod.main()
    finally:
        sys.argv = _saved


def call_convert_ribotish(argv):
    """Call ribotish_to_gencode.py main() directly."""
    _ensure_scripts_on_path()
    sys.path.insert(0, os.path.join(_scripts_dir(), 'gencode_converters'))
    _call_main_with_argv('ribotish_to_gencode', 'ribotish_to_gencode', argv)


def call_convert_ribotricer(argv):
    """Call ribotricer_to_gencode.py main() directly."""
    _ensure_scripts_on_path()
    sys.path.insert(0, os.path.join(_scripts_dir(), 'gencode_converters'))
    _call_main_with_argv('ribotricer_to_gencode', 'ribotricer_to_gencode', argv)
