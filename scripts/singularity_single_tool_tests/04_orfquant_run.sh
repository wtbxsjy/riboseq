#!/usr/bin/env bash
set -euo pipefail

# Avoid noisy locale warnings inside the container (host locales may not exist there)
export SINGULARITYENV_LANG=${SINGULARITYENV_LANG:-C}
export SINGULARITYENV_LC_ALL=${SINGULARITYENV_LC_ALL:-C}
export APPTAINERENV_LANG=${APPTAINERENV_LANG:-C}
export APPTAINERENV_LC_ALL=${APPTAINERENV_LC_ALL:-C}

# Mirrors: modules/local/orfquant/main.nf
# Default container in pipeline points to riboseqc:1.1 (ORFquant may be installed at runtime)
IMG_URL="https://depot.galaxyproject.org/singularity/riboseqc:1.1--r36_1"

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
CPUS=4
OUTDIR="./out_orfquant"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --sample) SAMPLE="$2"; shift 2;;
    --for-orfquant) FOR_ORFQUANT="$2"; shift 2;;
    --annotation) ANNOT="$2"; shift 2;;
    --fasta) FASTA="$2"; shift 2;;
    --orfquant-pkg) ORFQUANT_PKG="$2"; shift 2;;
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

pull_img() {
  local url="$1"
  local base
  base="$(basename "$url")"
  base="${base//:/_}"
  local sif="$(pwd)/containers/${base}.sif"
  if [[ ! -f "$sif" ]]; then
    singularity pull --disable-cache --force "$sif" "$url"
  fi
  echo "$sif"
}

IMG="$(pull_img "$IMG_URL")"

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

export R_LIBS_USER="$(pwd)/Rlibs"
mkdir -p "$R_LIBS_USER"

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

# Ensure RiboseQC is loaded here so we can patch its imports before ORFquant runs.
suppressPackageStartupMessages(library(RiboseQC))

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
repair_imported_symbol("RiboseQC", "Position", "base")

# combine is a BiocGenerics generic; some setups import gridExtra::combine instead.
repair_imported_symbol("ORFquant", "combine", "BiocGenerics")
repair_imported_symbol("RiboseQC", "combine", "BiocGenerics")

# Fail fast with a clear message if the repair did not take effect.
check_import_is_function <- function(pkg, sym) {
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
check_import_is_function("RiboseQC", "Position")

run_ORFquant(
  for_ORFquant_file = Sys.getenv("FOR_ORFQUANT"),
  annotation_file   = Sys.getenv("ANNOT"),
  n_cores           = as.integer(Sys.getenv("CPUS")),
  prefix            = Sys.getenv("SAMPLE"),
  write_temp_files  = TRUE,
  write_GTF_file    = TRUE,
  write_protein_fasta = TRUE,
  interactive       = FALSE
)

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

singularity exec \
  --bind "$WORKDIR:$WORKDIR${BIND_EXTRA:+,$BIND_EXTRA}" \
  --pwd "$WORKDIR" \
  "$IMG" \
  bash "$INNER_SH"

echo "[OK] ORFquant outputs in: $OUTDIR (look for ${SAMPLE}_final_ORFquant_results)"
