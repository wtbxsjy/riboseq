#!/bin/bash
# ─── 在 unify_orf.sif 容器内执行 run_orf.py（自动处理 duckdb/pyuser 细节）───
# Usage:
#   run_orf_in_container.sh <sif> [options] -- <run_orf.py 子命令> [args...]
#
# Options (在 -- 之前):
#   --runtime-dir DIR   可写运行时目录（容器 HOME/PYTHONUSERBASE，默认 ./orf_runtime）
#   --bind PATH         bind 进容器的路径（默认仓库根）
#   --no-duckdb         跳过 pip install duckdb（orfont 不可用/原脚本路径）
#   --dry-run           只打印将要执行的 singularity 命令
#   --log FILE          同时 tee 输出到日志文件
#
# Examples:
#   run_orf_in_container.sh run/rice/containers/unify_orf.sif \
#       -- unify --gtf ref.gtf --fasta ref.fa --output unified_orfs --ribotish '*.txt'
#   run_orf_in_container.sh run/rice/containers/unify_orf.sif --dry-run \
#       -- classify-gencode --input unified_orfs --output_dir out --ensembl_dir Ens58

set -euo pipefail

[ $# -ge 2 ] || { echo "Usage: $0 <sif> [options] -- <subcommand> [args...]"; exit 2; }

SIF="$1"; shift
[ -f "$SIF" ] || { echo "ERROR: sif not found: $SIF"; exit 2; }

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../../.." && pwd)"
RUN_ORF="$REPO_ROOT/scripts/run_orf.py"

RUNTIME_DIR="$(pwd)/orf_runtime"
BIND="$REPO_ROOT"
INSTALL_DUCKDB=true
DRY_RUN=false
LOG=""

while [ $# -gt 0 ]; do
  case "$1" in
    --runtime-dir) RUNTIME_DIR="$2"; shift 2 ;;
    --bind)        BIND="$2"; shift 2 ;;
    --no-duckdb)   INSTALL_DUCKDB=false; shift ;;
    --dry-run)     DRY_RUN=true; shift ;;
    --log)         LOG="$2"; shift 2 ;;
    --)            shift; break ;;
    *) echo "ERROR: unknown option: $1 (use -- before run_orf args)"; exit 2 ;;
  esac
done

[ $# -ge 1 ] || { echo "ERROR: missing run_orf subcommand after --"; exit 2; }

mkdir -p "$RUNTIME_DIR"
case "$RUNTIME_DIR" in
  "$BIND"/*) ;;
  *) echo "WARNING: runtime dir $RUNTIME_DIR not inside bind path $BIND (container may not see it)" ;;
esac

INNER="export HOME='$RUNTIME_DIR'; "
INNER+="export PYTHONUSERBASE='$RUNTIME_DIR/.pylibs'; "
INNER+="export PATH=\"\$PYTHONUSERBASE/bin:\$PATH\"; "
INNER+="export PIP_NO_CACHE_DIR=1; "
$INSTALL_DUCKDB && INNER+="pip install --user --no-cache-dir duckdb >/dev/null 2>&1 || true; "
INNER+="python3 -u '$RUN_ORF' $(printf '%q ' "$@")"

CMD=(singularity exec --no-home --pid -B "$BIND" "$SIF" bash -c "$INNER")

if $DRY_RUN; then
  echo "sif:       $SIF"
  echo "bind:      $BIND"
  echo "runtime:   $RUNTIME_DIR"
  echo "command:"
  printf '  %s\n' "${CMD[@]}" | head -c 2000
  echo ""
  exit 0
fi

echo "[run_orf_in_container] sif=$SIF runtime=$RUNTIME_DIR"
echo "[run_orf_in_container] $*"
if [ -n "$LOG" ]; then
  "${CMD[@]}" 2>&1 | tee "$LOG"
else
  "${CMD[@]}"
fi
