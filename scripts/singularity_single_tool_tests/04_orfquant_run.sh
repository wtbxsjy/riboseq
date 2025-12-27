#!/usr/bin/env bash
set -euo pipefail

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

singularity exec \
  --bind "$WORKDIR:$WORKDIR${BIND_EXTRA:+,$BIND_EXTRA}" \
  --pwd "$WORKDIR" \
  "$IMG" \
  bash -lc "
set -euo pipefail
cd '$OUTDIR'

# Handle gz fasta like module does
FASTA_IN='$FASTA'
if [[ \"\$FASTA_IN\" == *.gz ]]; then
  gunzip -c \"\$FASTA_IN\" > \"\$(basename \"\$FASTA_IN\" .gz)\"
  FASTA_IN=\"\$(basename \"\$FASTA_IN\" .gz)\"
fi

export R_LIBS_USER=\"$OUTDIR/Rlibs\"
mkdir -p \"\$R_LIBS_USER\"

cat > run_orfquant.R <<'RSCRIPTEOF'
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

export SAMPLE='$SAMPLE'
export FOR_ORFQUANT='$FOR_ORFQUANT'
export ANNOT='$ANNOT'
export CPUS='$CPUS'
export ORFQUANT_PKG='${ORFQUANT_PKG}'

Rscript run_orfquant.R
"

echo "[OK] ORFquant outputs in: $OUTDIR (look for ${SAMPLE}_final_ORFquant_results)"
