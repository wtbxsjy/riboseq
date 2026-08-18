#!/bin/bash
# ─── Nextflow pipeline 后台监控 ──────────────────────────────────────────
# Usage:
#   monitor_pipeline.sh <logfile> [--interval 300] [--timeout-hours 0] [--once]
#
#   <logfile>       .nextflow.log 路径（run/{project}/process/.nextflow.log）
#   --interval N    轮询间隔秒数（默认 300）
#   --timeout-hours H  超过 H 小时仍未完成则以非零码退出（0 = 不设限，默认 0）
#   --once          只打印一次统计后退出（快速查询）
#
# 行为：持续展示新增的进度行（Submitted/Cached/ERROR 等），
#       出现 "Pipeline completed" 时以 0 退出；出错行数持续增长时给出提示。
#
# Examples:
#   monitor_pipeline.sh run/rice/process/.nextflow.log --interval 300 --timeout-hours 96
#   monitor_pipeline.sh run/rice/process/.nextflow.log --once

set -euo pipefail

[ $# -ge 1 ] || { echo "Usage: monitor_pipeline.sh <logfile> [--interval 300] [--timeout-hours H] [--once]"; exit 1; }

LOG="$1"; shift
INTERVAL=300
TIMEOUT_HOURS=0
ONCE=false

while [ $# -gt 0 ]; do
  case "$1" in
    --interval)       INTERVAL="$2"; shift 2 ;;
    --timeout-hours)  TIMEOUT_HOURS="$2"; shift 2 ;;
    --once)           ONCE=true; shift ;;
    *) echo "ERROR: unknown option: $1"; exit 1 ;;
  esac
done

[ -f "$LOG" ] || { echo "ERROR: logfile not found: $LOG"; exit 1; }

PATTERN="Submitted|Cached|ERROR|Failed|Pipeline completed|Workflow execution interrupted"
START_TIME=$(date +%s)

print_stats() {
  local submitted errors completed
  submitted=$(grep -c "Submitted" "$LOG" 2>/dev/null || true)
  errors=$(grep -c "ERROR" "$LOG" 2>/dev/null || true)
  completed=$(grep -c "Pipeline completed" "$LOG" 2>/dev/null || true)
  echo "[$(date '+%H:%M:%S')] submitted=$submitted errors=$errors completed=$completed"
}

# 基线行号：只展示监控期间的新增行
BASELINE=$(wc -l < "$LOG")
print_stats

while true; do
  $ONCE && break
  sleep "$INTERVAL"

  NEW_LINES=$(wc -l < "$LOG")
  if [ "$NEW_LINES" -gt "$BASELINE" ]; then
    tail -n "+$((BASELINE + 1))" "$LOG" | grep -E "$PATTERN" | grep -v DEBUG || true
    BASELINE=$NEW_LINES
  fi

  if grep -q "Pipeline completed" "$LOG"; then
    echo "PIPELINE COMPLETED"
    print_stats
    exit 0
  fi

  if [ "$TIMEOUT_HOURS" -gt 0 ]; then
    ELAPSED=$(( ($(date +%s) - START_TIME) / 3600 ))
    if [ "$ELAPSED" -ge "$TIMEOUT_HOURS" ]; then
      echo "TIMEOUT: ${TIMEOUT_HOURS}h elapsed without completion"
      print_stats
      exit 2
    fi
  fi
done
