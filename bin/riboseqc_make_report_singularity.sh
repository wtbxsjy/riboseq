#!/usr/bin/env bash
# If this script is invoked with /bin/sh (e.g. `sh script.sh`), Bash-only options
# like `set -o pipefail` will fail. Re-exec with bash in that case.
if [ -z "${BASH_VERSION:-}" ]; then
  exec /usr/bin/env bash "$0" "$@"
fi

set -euo pipefail

usage() {
  cat <<'USAGE'
Usage:
  riboseqc_make_report_singularity.sh \
    --results-dir /path/to/riboseq_qc/riboseqc \
    --output-html /path/to/report.html \
    [--image IMAGE_OR_URI] \
    [--copy-plots] \
    [--extended]

What it does:
  - Finds all "*_results_RiboseQC" files under --results-dir (excluding "*_results_RiboseQC_all")
  - Runs RiboseQC::create_html_report() inside a Singularity container
  - Writes the HTML report to --output-html, plus a "*_plots/" folder next to it

Defaults:
  --image docker://quay.io/biocontainers/riboseqc:1.1--r36_1

Notes:
  - This generates the HTML report FROM existing RiboseQC results files.
  - It does not rerun RiboseQC_analysis(), so it avoids annotation/genome-package incompatibilities.
  - Font rendering: the container image is minimal and may not ship fonts.
    This script will (when present) bind host fontconfig + fonts into the container:
      /etc/fonts, /usr/share/fonts, /usr/local/share/fonts
  - On Windows-mounted paths (e.g. /mnt/c/... OneDrive), RiboseQC can create very long filenames
    under "*_plots/" which may hit path-length or sync/locking issues.
    If staging is enabled, this script copies the HTML back by default; add --copy-plots to also
    copy the "*_plots/" folder back.
USAGE
}

IMAGE='docker://quay.io/biocontainers/riboseqc:1.1--r36_1'
RESULTS_DIR=''
OUTPUT_HTML=''
EXTENDED='FALSE'
STAGE_DIR=''
COPY_PLOTS='FALSE'

while [[ $# -gt 0 ]]; do
  case "$1" in
    --results-dir)
      RESULTS_DIR=${2:-}; shift 2 ;;
    --output-html)
      OUTPUT_HTML=${2:-}; shift 2 ;;
    --image)
      IMAGE=${2:-}; shift 2 ;;
    --extended)
      EXTENDED='TRUE'; shift 1 ;;
    --copy-plots)
      COPY_PLOTS='TRUE'; shift 1 ;;
    -h|--help)
      usage; exit 0 ;;
    *)
      echo "Unknown argument: $1" >&2
      usage >&2
      exit 2
      ;;
  esac
done

if [[ -z "$RESULTS_DIR" || -z "$OUTPUT_HTML" ]]; then
  usage >&2
  exit 2
fi

if [[ ! -d "$RESULTS_DIR" ]]; then
  echo "ERROR: --results-dir is not a directory: $RESULTS_DIR" >&2
  exit 1
fi

# Make paths absolute for safer bind mounting
RESULTS_DIR=$(readlink -f "$RESULTS_DIR")
OUTPUT_HTML=$(readlink -m "$OUTPUT_HTML")
OUTPUT_DIR=$(dirname "$OUTPUT_HTML")
OUTPUT_BASE=$(basename "$OUTPUT_HTML")

mkdir -p "$OUTPUT_DIR"

# Preflight: avoid the common failure mode where the intended HTML output path
# already exists as a directory (then rmarkdown cannot write the HTML file).
if [[ -d "$OUTPUT_HTML" ]]; then
  echo "ERROR: --output-html points to an existing DIRECTORY, not a file: $OUTPUT_HTML" >&2
  echo "       Please choose a different filename or remove/rename that directory." >&2
  exit 1
fi

# Preflight: ensure we can write to the output directory
if [[ ! -w "$OUTPUT_DIR" ]]; then
  echo "ERROR: Output directory is not writable: $OUTPUT_DIR" >&2
  exit 1
fi

# If an old HTML file exists, remove it to avoid partial/locked overwrites.
if [[ -f "$OUTPUT_HTML" ]]; then
  rm -f "$OUTPUT_HTML" || true
fi

# Collect results files (RData-like) produced by RiboseQC_analysis
mapfile -t RESULTS_FILES < <(
  find "$RESULTS_DIR" -maxdepth 2 -type f -name '*_results_RiboseQC' ! -name '*_results_RiboseQC_all' | sort
)

if [[ ${#RESULTS_FILES[@]} -eq 0 ]]; then
  echo "ERROR: No '*_results_RiboseQC' files found under: $RESULTS_DIR" >&2
  exit 1
fi

# Create an R script for report rendering.
# Use short sample names (<=5 chars) to satisfy RiboseQC's expectations.
R_SCRIPT=$(mktemp)

cleanup() {
  rm -f "$R_SCRIPT"
  if [[ -n "${STAGE_DIR:-}" && -d "$STAGE_DIR" ]]; then
    rm -rf "$STAGE_DIR"
  fi
}
trap cleanup EXIT

# Apptainer/Singularity parses --bind as a comma-separated list.
# If host paths contain commas (common on Windows/OneDrive, e.g. "Co., Limited"),
# binds will break even when quoted. Work around by staging into a comma-free temp dir.
NEEDS_STAGING='FALSE'
if [[ "$RESULTS_DIR" == *","* || "$OUTPUT_DIR" == *","* ]]; then
  NEEDS_STAGING='TRUE'
fi

# Even without commas, writing many long plot filenames to Windows-mounted paths
# can fail (path-length / sync / locking). If output is on /mnt/*, prefer staging.
if [[ "$OUTPUT_DIR" == /mnt/* ]]; then
  NEEDS_STAGING='TRUE'
fi

IN_DIR_FOR_CONTAINER="$RESULTS_DIR"
OUT_DIR_FOR_CONTAINER="$OUTPUT_DIR"

if [[ "$NEEDS_STAGING" == 'TRUE' ]]; then
  STAGE_DIR=$(mktemp -d -p /tmp riboseqc_report_stage.XXXXXX)
  STAGE_IN="$STAGE_DIR/in"
  STAGE_OUT="$STAGE_DIR/out"
  mkdir -p "$STAGE_IN" "$STAGE_OUT"

  echo "[riboseqc] Using staging under: $STAGE_DIR" >&2
  if [[ "$RESULTS_DIR" == *","* || "$OUTPUT_DIR" == *","* ]]; then
    echo "[riboseqc] Reason: comma in path breaks Apptainer --bind parsing." >&2
  elif [[ "$OUTPUT_DIR" == /mnt/* ]]; then
    echo "[riboseqc] Reason: output on Windows mount (/mnt/*); avoiding long-filename/sync issues in *_plots." >&2
  fi

  # Copy input results to staging (avoid commas in bind path)
  idx=1
  STAGED_BASENAMES=()
  for f in "${RESULTS_FILES[@]}"; do
    b=$(basename "$f")
    staged_name="$b"
    # Ensure unique basenames in staging
    if [[ -e "$STAGE_IN/$staged_name" ]]; then
      staged_name="S${idx}_$b"
    fi
    cp -f "$f" "$STAGE_IN/$staged_name"
    STAGED_BASENAMES+=("$staged_name")
    idx=$((idx+1))
  done

  IN_DIR_FOR_CONTAINER="$STAGE_IN"
  OUT_DIR_FOR_CONTAINER="$STAGE_OUT"
fi

{
  echo 'suppressPackageStartupMessages(library(RiboseQC))'
  echo 'suppressPackageStartupMessages(library(rmarkdown))'
  echo ''
  echo '# input_files: named character vector, names are sample labels in the report'
  echo 'input_files <- c()'

  idx=1
  if [[ "${NEEDS_STAGING}" == 'TRUE' ]]; then
    for b in "${STAGED_BASENAMES[@]}"; do
      echo "input_files <- c(input_files, S${idx}=file.path('/in', '${b}'))"
      idx=$((idx+1))
    done
  else
    for f in "${RESULTS_FILES[@]}"; do
      b=$(basename "$f")
      # Inside container we bind results dir to /in
      echo "input_files <- c(input_files, S${idx}=file.path('/in', '${b}'))"
      idx=$((idx+1))
    done
  fi

  echo ''
  echo 'cat("RiboseQC version:", as.character(packageVersion("RiboseQC")), "\n")'
  echo 'cat("Samples:\n")'
  echo 'for (nm in names(input_files)) cat("  ", nm, " -> ", input_files[[nm]], "\n", sep="")'
  echo ''

  echo '# NOTE: RiboseQC 1.1 create_html_report() renders the Rmd template in-place from the'
  echo '# package library directory and tries to write "riboseqc_template.knit.md" there.'
  echo '# In containerized environments the library directory is typically read-only, causing:'
  echo '#   Error in file(con, "w") : cannot open the connection'
  echo '# Workaround: copy the template to a writable tempdir and render that copy.'
  echo ''
  echo "output_file <- file.path('/out', '${OUTPUT_BASE}')"
  echo 'output_fig_path <- paste0(output_file, "_plots/")'
  echo 'dir.create(paste0(output_fig_path, "rds/"), recursive=TRUE, showWarnings=FALSE)'
  echo 'dir.create(paste0(output_fig_path, "pdf/"), recursive=TRUE, showWarnings=FALSE)'
  echo 'sink(file = paste0(output_file, "_report_text_output.txt"))'
  echo 'on.exit({ while (sink.number() > 0) sink(NULL) }, add = TRUE)'

  echo 'rmd_path <- file.path(system.file(package="RiboseQC", mustWork = TRUE), "rmd", "riboseqc_template.Rmd")'
  echo "if (${EXTENDED}) rmd_path <- file.path(system.file(package=\"RiboseQC\", mustWork = TRUE), \"rmd\", \"riboseqc_template_full.Rmd\")"
  echo 'rmd_tmp <- file.path(tempdir(), basename(rmd_path))'
  echo 'file.copy(rmd_path, rmd_tmp, overwrite=TRUE)'
  echo 'setwd(tempdir())'

  echo 'render(rmd_tmp,'
  echo '  params = list('
  echo '    input_files = input_files,'
  echo '    input_sample_names = names(input_files),'
  echo '    output_fig_path = output_fig_path'
  echo '  ),'
  echo '  output_file = output_file,'
  echo '  intermediates_dir = tempdir(),' 
  echo '  quiet = FALSE'
  echo ')'
} > "$R_SCRIPT"

# Bind-mount the exact results/output directories.
# (Windows/OneDrive paths: spaces are ok as long as quoted.)
BIND_ARGS=(
  --bind "$IN_DIR_FOR_CONTAINER:/in"
  --bind "$OUT_DIR_FOR_CONTAINER:/out"
)

# If available on the host, bind fonts + fontconfig into the container.
# This fixes missing/garbled text in plots when the container lacks system fonts.
if [[ -d /etc/fonts ]]; then
  BIND_ARGS+=(--bind "/etc/fonts:/etc/fonts")
fi
if [[ -d /usr/share/fonts ]]; then
  BIND_ARGS+=(--bind "/usr/share/fonts:/usr/share/fonts")
fi
if [[ -d /usr/local/share/fonts ]]; then
  BIND_ARGS+=(--bind "/usr/local/share/fonts:/usr/local/share/fonts")
fi

singularity exec \
  --cleanenv \
  --home /tmp \
  "${BIND_ARGS[@]}" \
  "$IMAGE" \
  env \
    R_PROFILE_USER=/dev/null \
    R_ENVIRON_USER=/dev/null \
    FONTCONFIG_PATH=/etc/fonts \
    FONTCONFIG_FILE=/etc/fonts/fonts.conf \
    Rscript "$R_SCRIPT"

# If we staged, copy results back to the requested output directory.
if [[ "$NEEDS_STAGING" == 'TRUE' ]]; then
  staged_html="$OUT_DIR_FOR_CONTAINER/$OUTPUT_BASE"
  staged_plots_dir="$OUT_DIR_FOR_CONTAINER/${OUTPUT_BASE}_plots"

  if [[ ! -f "$staged_html" ]]; then
    echo "ERROR: Expected report not found in staging output: $staged_html" >&2
    exit 1
  fi

  cp -f "$staged_html" "$OUTPUT_DIR/$OUTPUT_BASE"

  if [[ -d "$staged_plots_dir" ]]; then
    if [[ "$COPY_PLOTS" == 'TRUE' ]]; then
      rm -rf "$OUTPUT_DIR/${OUTPUT_BASE}_plots" || true
      cp -a "$staged_plots_dir" "$OUTPUT_DIR/${OUTPUT_BASE}_plots"
    else
      echo "[riboseqc] NOTE: Not copying plots folder back (use --copy-plots to enable):" >&2
      echo "[riboseqc]   $staged_plots_dir" >&2
    fi
  fi
fi
