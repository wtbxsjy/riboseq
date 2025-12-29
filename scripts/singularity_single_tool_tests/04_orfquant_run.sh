#!/usr/bin/env bash
set -euo pipefail

# Avoid noisy locale warnings inside the container (host locales may not exist there)
export SINGULARITYENV_LANG=${SINGULARITYENV_LANG:-C}
export SINGULARITYENV_LC_ALL=${SINGULARITYENV_LC_ALL:-C}
export APPTAINERENV_LANG=${APPTAINERENV_LANG:-C}
export APPTAINERENV_LC_ALL=${APPTAINERENV_LC_ALL:-C}

# Reference-aligned default:
# The provided Nextflow example uses a dedicated ORFquant biocontainer (R 4.0).
# This generally avoids runtime installation and reduces namespace/version issues.
IMG_URL="https://depot.galaxyproject.org/singularity/orfquant:1.1.0--r40_1"

usage() {
  cat <<'EOF'
Usage:
  04_orfquant_run.sh --sample ID --for-orfquant X_for_ORFquant --annotation X_Rannot --fasta genome.fa [options]

Required:
  --sample        sample ID (prefix)
  --for-orfquant  RiboseQC output file (*_for_ORFquant)
  --annotation    RiboseQC annotation file (*_Rannot)
  --fasta         genome FASTA (can be .gz)

Options:
  --orfquant-pkg  local ORFquant source tar.gz (optional; avoids GitHub download)
  --container     SIF file path to run ORFquant in (optional). If not set, will
                  use ./orfquant.sif when present, otherwise pull riboseqc:1.1.
  --cpus          threads (default: 4)
  --outdir        output directory (default: ./out_orfquant)

Env:
  BIND_EXTRA extra singularity binds, comma-separated (e.g. /mnt:/mnt)
EOF
}

SAMPLE=""
FOR_ORFQUANT=""
ANNOT=""
FASTA=""
ORFQUANT_PKG=""
CONTAINER_SIF=""
CPUS=4
OUTDIR="./out_orfquant"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --sample) SAMPLE="$2"; shift 2;;
    --for-orfquant) FOR_ORFQUANT="$2"; shift 2;;
    --annotation) ANNOT="$2"; shift 2;;
    --fasta) FASTA="$2"; shift 2;;
    --orfquant-pkg) ORFQUANT_PKG="$2"; shift 2;;
    --container) CONTAINER_SIF="$2"; shift 2;;
    --cpus) CPUS="$2"; shift 2;;
    --outdir) OUTDIR="$2"; shift 2;;
    -h|--help) usage; exit 0;;
    *) echo "Unknown arg: $1"; usage; exit 2;;
  esac
done

if [[ -z "$SAMPLE" || -z "$FOR_ORFQUANT" || -z "$ANNOT" || -z "$FASTA" ]]; then
  usage
  exit 2
fi

if [[ -n "$ORFQUANT_PKG" ]]; then
  if [[ ! -f "$ORFQUANT_PKG" ]]; then
    echo "[ERROR] --orfquant-pkg not found: $ORFQUANT_PKG" >&2
    exit 2
  fi
  ORFQUANT_PKG="$(python3 -c 'import os,sys; print(os.path.abspath(sys.argv[1]))' "$ORFQUANT_PKG")"
fi

mkdir -p "$OUTDIR" "./containers"
OUTDIR="$(cd "$OUTDIR" && pwd)"
WORKDIR="$(pwd)"

detect_runtime() {
  if command -v apptainer >/dev/null 2>&1; then
    echo "apptainer"
    return 0
  fi
  if command -v singularity >/dev/null 2>&1; then
    echo "singularity"
    return 0
  fi
  echo "[ERROR] Neither 'apptainer' nor 'singularity' is available on PATH." >&2
  echo "        Install Apptainer/Singularity or run on a node that provides it." >&2
  exit 2
}

RUNTIME="$(detect_runtime)"

pull_img() {
  local url="$1"
  local base
  base="$(basename "$url")"
  base="${base//:/_}"
  local sif="$(pwd)/containers/${base}.sif"
  if [[ ! -f "$sif" ]]; then
    "$RUNTIME" pull --disable-cache --force "$sif" "$url"
  fi
  echo "$sif"
}

if [[ -n "$CONTAINER_SIF" ]]; then
  if [[ ! -f "$CONTAINER_SIF" ]]; then
    echo "[ERROR] --container not found: $CONTAINER_SIF" >&2
    exit 2
  fi
  IMG="$(python3 -c 'import os,sys; print(os.path.abspath(sys.argv[1]))' "$CONTAINER_SIF")"
elif [[ -f "$WORKDIR/orfquant.sif" ]]; then
  IMG="$WORKDIR/orfquant.sif"
else
  IMG="$(pull_img "$IMG_URL")"
fi

echo "[INFO] Container runtime: $RUNTIME"
echo "[INFO] Container image:   $IMG"

collect_binds() {
  local binds=()
  local add
  add() {
    local p="$1"
    [[ -z "$p" ]] && return 0
    if [[ -e "$p" ]]; then
      local d
      if [[ -d "$p" ]]; then
        d="$p"
      else
        d="$(cd "$(dirname "$p")" && pwd)"
      fi
      binds+=("$d:$d")
    fi
  }

  add "$WORKDIR"
  add "$OUTDIR"
  add "$FOR_ORFQUANT"
  add "$ANNOT"
  add "$FASTA"
  add "$ORFQUANT_PKG"

  # Add any extra user binds.
  if [[ -n "${BIND_EXTRA:-}" ]]; then
    IFS=',' read -r -a extra <<<"$BIND_EXTRA"
    for b in "${extra[@]}"; do
      [[ -n "$b" ]] && binds+=("$b")
    done
  fi

  # Deduplicate binds while preserving order.
  python3 - "$@" <<'PY'
import sys
seen=set()
out=[]
for b in sys.argv[1:]:
    if b not in seen:
        seen.add(b)
        out.append(b)
print(",".join(out))
PY
}

BIND_SPEC="$(collect_binds "${WORKDIR}:${WORKDIR}" "${OUTDIR}:${OUTDIR}" )"

# Write a small wrapper to execute inside the container (avoids fragile nested quoting)
quote_sh() {
  python3 -c 'import sys,shlex; s=sys.stdin.read().rstrip("\n"); print(shlex.quote(s))' <<<"$1"
}

ENV_SH="$OUTDIR/orfquant_env.sh"
cat > "$ENV_SH" <<EOF
export SAMPLE=$(quote_sh "$SAMPLE")
export FOR_ORFQUANT=$(quote_sh "$FOR_ORFQUANT")
export ANNOT=$(quote_sh "$ANNOT")
export FASTA=$(quote_sh "$FASTA")
export CPUS=$(quote_sh "$CPUS")
export ORFQUANT_PKG=$(quote_sh "$ORFQUANT_PKG")
EOF

INNER_SH="$OUTDIR/run_orfquant_container.sh"
cat > "$INNER_SH" <<'EOSH'
#!/usr/bin/env bash
set -euo pipefail

cd "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source ./orfquant_env.sh

# Handle gz fasta like module does
FASTA_IN="$FASTA"
if [[ "$FASTA_IN" == *.gz ]]; then
  gunzip -c "$FASTA_IN" > "$(basename "$FASTA_IN" .gz)"
  FASTA_IN="$(basename "$FASTA_IN" .gz)"
fi

# Append task-local Rlibs to R_LIBS_USER instead of replacing, so that packages
# pre-installed in the container (e.g. ORFquant) remain accessible.
_local_rlibs="$(pwd)/Rlibs"
mkdir -p "$_local_rlibs"
export R_LIBS_USER="${_local_rlibs}${R_LIBS_USER:+:${R_LIBS_USER}}"

# Prevent host/user R startup files from interfering inside the container.
export R_PROFILE_USER=/dev/null
export R_ENVIRON_USER=/dev/null

cat > run_orfquant.R <<'RSCRIPTEOF'
options(show.error.locations = TRUE)
options(error = function() {
  cat("\n=== R ERROR ===\n")
  cat(geterrmessage(), "\n")
  cat("\n=== sys.calls() ===\n")
  print(sys.calls())
  cat("\n=== traceback() ===\n")
  try(traceback(50), silent = TRUE)
  cat("\n=== sessionInfo() ===\n")
  print(sessionInfo())
  quit(save = "no", status = 1)
})

# Some environments can throw a fatal error from R6 finalizers during garbage
# collection, e.g. "x$.self$finalize() : attempt to apply non-function".
# This is usually unrelated to the analysis logic but aborts the whole run.
# We defensively wrap R6-registered finalizers so they cannot terminate the job.
patch_r6_finalizers <- function() {
  if (!requireNamespace("R6", quietly = TRUE)) return(invisible(FALSE))
  safe_reg_finalizer <- function(e, f, onexit = FALSE) {
    wrapped <- function(x) {
      tryCatch(f(x), error = function(err) {
        message("[WARN] Suppressed error in finalizer: ", conditionMessage(err))
        invisible(NULL)
      })
    }
    environment(wrapped) <- baseenv()
    base::reg.finalizer(e, wrapped, onexit = onexit)
  }
  ok <- tryCatch({
    assignInNamespace("reg.finalizer", safe_reg_finalizer, ns = "R6")
    TRUE
  }, error = function(e) FALSE)
  invisible(ok)
}

patch_r6_finalizers()

install_orfquant <- function(local_pkg_tgz = NULL, tag = "1.02") {
  work <- file.path(getwd(), "orfquant_src")
  dir.create(work, showWarnings = FALSE, recursive = TRUE)

  tgz <- local_pkg_tgz
  if (!is.null(tgz) && nzchar(tgz) && file.exists(tgz) && file.info(tgz)[1, "size"] > 0) {
    message("Installing ORFquant from local tar.gz: ", tgz)
  } else {
    url <- sprintf("https://github.com/lcalviell/ORFquant/archive/refs/tags/%s.tar.gz", tag)
    tgz <- file.path(work, sprintf("ORFquant-%s.tar.gz", tag))
    message("Downloading ORFquant from GitHub: ", url)
    utils::download.file(url, tgz, mode = "wb", quiet = FALSE)
  }

  utils::untar(tgz, exdir = work, tar = "internal")
  pkg_dir <- list.dirs(work, recursive = FALSE, full.names = TRUE)
  if (length(pkg_dir) != 1) stop("Unexpected ORFquant source layout in: ", work)

  cmd <- sprintf("R CMD INSTALL %s", shQuote(pkg_dir[[1]]))
  message(cmd)
  status <- system(cmd)
  if (status != 0) stop("R CMD INSTALL failed with status ", status)
}

if (!requireNamespace("ORFquant", quietly = TRUE)) {
  local_pkg <- Sys.getenv("ORFQUANT_PKG", unset = "")
  if (!nzchar(local_pkg)) local_pkg <- NULL
  tryCatch({
    install_orfquant(local_pkg_tgz = local_pkg, tag = "1.02")
  }, error = function(e) {
    stop(
      "ORFquant is not installed and automatic installation failed: ", conditionMessage(e), "\n",
      "Provide a pre-downloaded tarball with --orfquant-pkg, or use a container with ORFquant pre-installed."
    )
  })
}

library(ORFquant)

# Try to load RiboseQC if available (some containers like orfquant:1.1.0--r40_1 don't include it).
# RiboseQC is only needed for patching import collisions; ORFquant can run without it.
has_riboseqc <- requireNamespace("RiboseQC", quietly = TRUE)
if (has_riboseqc) {
  suppressPackageStartupMessages(library(RiboseQC))
} else {
  message("Note: RiboseQC package not available. Skipping RiboseQC-specific import patches.")
}

# Work around import-environment collisions observed in some environments:
# - ggplot2::Position is a ggproto object (not a function) and can override base/BiocGenerics::Position,
#   causing: "attempt to apply non-function".
#
# Note: these conflicting symbols live in the *imports environment* (parent of the namespace), not
# necessarily in the namespace environment itself.
repair_imported_symbol <- function(pkg, sym, replacement_env, replacement_sym = sym) {
  ns <- tryCatch(asNamespace(pkg), error = function(e) NULL)
  if (is.null(ns)) return(invisible(FALSE))
  imp <- parent.env(ns)
  if (!exists(sym, envir = imp, inherits = FALSE)) return(invisible(FALSE))
  obj <- get(sym, envir = imp, inherits = FALSE)
  if (is.function(obj)) return(invisible(FALSE))

  repl <- if (identical(replacement_env, "base")) {
    get(replacement_sym, envir = baseenv(), inherits = FALSE)
  } else {
    get(replacement_sym, envir = asNamespace(replacement_env), inherits = FALSE)
  }

  was_locked <- FALSE
  if (bindingIsLocked(sym, imp)) {
    was_locked <- TRUE
    unlockBinding(sym, imp)
  }
  assign(sym, repl, envir = imp)
  if (was_locked) lockBinding(sym, imp)
  invisible(TRUE)
}

# Patch both ORFquant and RiboseQC imports.
repair_imported_symbol("ORFquant", "Position", "base")
if (has_riboseqc) {
  repair_imported_symbol("RiboseQC", "Position", "base")
}

# combine is a BiocGenerics generic; some setups import gridExtra::combine instead.
repair_imported_symbol("ORFquant", "combine", "BiocGenerics")
if (has_riboseqc) {
  repair_imported_symbol("RiboseQC", "combine", "BiocGenerics")
}

# Fail fast with a clear message if the repair did not take effect.
check_import_is_function <- function(pkg, sym) {
  if (!requireNamespace(pkg, quietly = TRUE)) return(invisible(TRUE))
  ns <- asNamespace(pkg)
  imp <- parent.env(ns)
  if (!exists(sym, envir = imp, inherits = FALSE)) return(invisible(TRUE))
  obj <- get(sym, envir = imp, inherits = FALSE)
  if (!is.function(obj)) {
    stop(sprintf("After repair, %s import '%s' is still not a function (class: %s)",
                 pkg, sym, paste(class(obj), collapse = ",")))
  }
  invisible(TRUE)
}

check_import_is_function("ORFquant", "Position")
if (has_riboseqc) {
  check_import_is_function("RiboseQC", "Position")
}

# Force gc() to trigger any pending finalizers BEFORE entering parallel code.
# This can reduce the chance of GC-triggered errors in child processes.
invisible(gc(verbose = FALSE, full = TRUE))

# IMPORTANT: Temporarily disable the global error handler so tryCatch can catch errors.
# The global error handler (options(error = ...)) calls quit() which prevents tryCatch
# from handling errors for the fallback mechanism.
saved_error_handler <- getOption("error")
options(error = NULL)

# Wrap run_ORFquant in a tryCatch; if it fails with the "attempt to apply non-function"
# error (common in parallel mode), retry with n_cores = 1.
n_cores_requested <- as.integer(Sys.getenv("CPUS"))
run_orfquant_result <- tryCatch({
  run_ORFquant(
    for_ORFquant_file = Sys.getenv("FOR_ORFQUANT"),
    annotation_file   = Sys.getenv("ANNOT"),
    n_cores           = n_cores_requested,
    prefix            = Sys.getenv("SAMPLE"),
    write_temp_files  = TRUE,
    write_GTF_file    = TRUE,
    write_protein_fasta = TRUE,
    interactive       = FALSE
  )
  "success"
}, error = function(e) {
  msg <- conditionMessage(e)
  # Check if it's the known parallel/import collision error
  if (grepl("attempt to apply non-function", msg, fixed = TRUE) && n_cores_requested > 1) {
    message("\n[WARN] ORFquant failed with 'attempt to apply non-function' error in parallel mode.")
    message("[WARN] This is a known R namespace collision issue when using doMC/mclapply.")
    message("[WARN] Retrying with n_cores = 1 (single-threaded mode)...\n")
    tryCatch({
      run_ORFquant(
        for_ORFquant_file = Sys.getenv("FOR_ORFQUANT"),
        annotation_file   = Sys.getenv("ANNOT"),
        n_cores           = 1L,
        prefix            = Sys.getenv("SAMPLE"),
        write_temp_files  = TRUE,
        write_GTF_file    = TRUE,
        write_protein_fasta = TRUE,
        interactive       = FALSE
      )
      "success_single_core"
    }, error = function(e2) {
      # Restore error handler before stopping
      options(error = saved_error_handler)
      stop("ORFquant failed even in single-core mode: ", conditionMessage(e2))
    })
  } else {
    # Restore error handler before stopping
    options(error = saved_error_handler)
    stop(e)
  }
})

# Restore the global error handler
options(error = saved_error_handler)

if (run_orfquant_result == "success_single_core") {
  message("[INFO] ORFquant completed successfully in single-core fallback mode.")
}

writeLines(
  c(
    'ORFQUANT_RUN:',
    paste0('    orfquant: "', packageVersion("ORFquant"), '"'),
    paste0('    r-base: "', R.Version()[["major"]], ".", R.Version()[["minor"]], '"')
  ),
  'versions.yml'
)
RSCRIPTEOF

echo [INFO] Running ORFquant R script...
echo [INFO] Container markers: APPTAINER_NAME=${APPTAINER_NAME:-} SINGULARITY_NAME=${SINGULARITY_NAME:-}
echo [INFO] which Rscript: $(command -v Rscript || true)
set -o pipefail
Rscript --vanilla run_orfquant.R 2>&1 | tee run_orfquant.log

if ls -1 ${SAMPLE}_final_ORFquant_results* >/dev/null 2>&1; then
  echo [INFO] ORFquant results:
  ls -la ${SAMPLE}_final_ORFquant_results* || true
else
  echo [ERROR] ORFquant finished but expected outputs were not found: ${SAMPLE}_final_ORFquant_results* >&2
  echo [ERROR] Contents of current directory: >&2
  pwd >&2 || true
  ls -la >&2 || true
  echo [ERROR] Tail of run_orfquant.log: >&2
  tail -n 200 run_orfquant.log >&2 || true
  exit 2
fi
EOSH

chmod +x "$INNER_SH"

"$RUNTIME" exec \
  --bind "$BIND_SPEC" \
  --pwd "$WORKDIR" \
  "$IMG" \
  bash "$INNER_SH"

echo "[OK] ORFquant outputs in: $OUTDIR (look for ${SAMPLE}_final_ORFquant_results)"
