#!/bin/bash
# ─── Nextflow 运行残留清理（默认干跑预览，--apply 才真删）────────────────
# Usage:
#   cleanup_nextflow.sh <process_dir> [options]
#
# Options:
#   --result-dir DIR      # 同时清理 DIR 下 report 文件（timeline/flowchart/pipeline_report.html）
#   --clear-lock          # 清理会话锁（.nextflow.pid + .nextflow/cache/*/db/LOCK）
#   --failed-only         # 删除 .exitcode 非 0 的失败任务 work dir
#   --task-pattern REGEX  # 删除 .command.sh 匹配 REGEX 的任务 work dir（如 "UNIFY_ORF|ORF_QC"）
#   --apply               # 真正执行删除；缺省只打印预览
#   --help
#
# Examples:
#   cleanup_nextflow.sh run/rice/process --result-dir run/rice/result --clear-lock
#   cleanup_nextflow.sh run/rice/process --failed-only
#   cleanup_nextflow.sh run/rice/process --task-pattern "UNIFY_ORF|ORF_QC" --apply

set -euo pipefail

apply=false
result_dir=""
clear_lock=false
failed_only=false
task_pattern=""

usage() {
  sed -n '2,16p' "$0" | sed 's/^# \{0,1\}//'
  exit "${1:-0}"
}

[ $# -ge 1 ] || usage 1
case "$1" in --help|-h) usage 0 ;; esac

PROCESS_DIR="$1"; shift
[ -d "$PROCESS_DIR" ] || { echo "ERROR: process dir not found: $PROCESS_DIR"; exit 1; }

while [ $# -gt 0 ]; do
  case "$1" in
    --result-dir)    result_dir="$2"; shift 2 ;;
    --clear-lock)    clear_lock=true; shift ;;
    --failed-only)   failed_only=true; shift ;;
    --task-pattern)  task_pattern="$2"; shift 2 ;;
    --apply)         apply=true; shift ;;
    *) echo "ERROR: unknown option: $1"; usage 1 ;;
  esac
done

WORK_DIR="$PROCESS_DIR/work"
NEXTFLOW_DIR="$PROCESS_DIR/.nextflow"

# ─── Collect actions ──────────────────────────────────────────────────────
rm_files=()
rm_dirs=()

if [ -n "$result_dir" ]; then
  for f in timeline.html flowchart.html pipeline_report.html; do
    [ -f "$result_dir/$f" ] && rm_files+=("$result_dir/$f")
  done
fi

if $clear_lock; then
  [ -f "$PROCESS_DIR/.nextflow.pid" ] && rm_files+=("$PROCESS_DIR/.nextflow.pid")
  if [ -d "$NEXTFLOW_DIR/cache" ]; then
    for lock in "$NEXTFLOW_DIR"/cache/*/db/LOCK; do
      [ -f "$lock" ] && rm_files+=("$lock")
    done
  fi
fi

if $failed_only && [ -d "$WORK_DIR" ]; then
  while IFS= read -r f; do
    [ -n "$f" ] && rm_dirs+=("$(dirname "$f")")
  done < <(find "$WORK_DIR" -name ".exitcode" -exec grep -l "^[^0]" {} \; 2>/dev/null)
fi

if [ -n "$task_pattern" ] && [ -d "$WORK_DIR" ]; then
  while IFS= read -r f; do
    [ -n "$f" ] && rm_dirs+=("$(dirname "$f")")
  done < <(find "$WORK_DIR" -name ".command.sh" -exec grep -l "$task_pattern" {} \; 2>/dev/null)
fi

# 去重
rm_dirs=($(printf '%s\n' "${rm_dirs[@]}" | sort -u))

# ─── Execute or preview ───────────────────────────────────────────────────
echo "process dir : $PROCESS_DIR"
echo "files to rm : ${#rm_files[@]}"
for f in "${rm_files[@]}"; do
  if $apply; then rm -f "$f"; echo "  RM   $f"; else echo "  would rm $f"; fi
done
echo "dirs to rm  : ${#rm_dirs[@]}"
for d in "${rm_dirs[@]}"; do
  if $apply; then rm -rf "$d"; echo "  RM   $d"; else echo "  would rm $d"; fi
done
$apply || echo "DRY RUN — pass --apply to actually delete."
