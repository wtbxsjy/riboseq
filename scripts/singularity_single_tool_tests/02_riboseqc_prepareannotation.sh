#!/usr/bin/env bash
set -euo pipefail

# Avoid noisy locale warnings inside the container (host locales may not exist there)
export SINGULARITYENV_LANG=${SINGULARITYENV_LANG:-C}
export SINGULARITYENV_LC_ALL=${SINGULARITYENV_LC_ALL:-C}
export APPTAINERENV_LANG=${APPTAINERENV_LANG:-C}
export APPTAINERENV_LC_ALL=${APPTAINERENV_LC_ALL:-C}

# Mirrors: modules/local/riboseqc/prepareannotation/main.nf
IMG_URL="https://depot.galaxyproject.org/singularity/riboseqc:1.1--r36_1"
SAMTOOLS_IMG_URL="https://depot.galaxyproject.org/singularity/samtools:1.21--h50ea8bc_0"

usage() {
  cat <<'EOF'
Usage:
  02_riboseqc_prepareannotation.sh --gtf annot.gtf --fasta genome.fa [--outdir DIR]

Required:
  --gtf    GTF (can be .gz)
  --fasta  genome FASTA (can be .gz)

Options:
  --outdir output directory (default: ./out_riboseqc_annot)

Env:
  BIND_EXTRA extra singularity binds, comma-separated (e.g. /mnt:/mnt)
EOF
}

GTF=""
FASTA=""
OUTDIR="./out_riboseqc_annot"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --gtf) GTF="$2"; shift 2;;
    --fasta) FASTA="$2"; shift 2;;
    --outdir) OUTDIR="$2"; shift 2;;
    -h|--help) usage; exit 0;;
    *) echo "Unknown arg: $1"; usage; exit 2;;
  esac
done

if [[ -z "$GTF" || -z "$FASTA" ]]; then
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

ensure_fai() {
  local fasta="$1"
  local img="$2"
  if [[ "$fasta" == *.gz ]]; then
    # RiboseQC prepare_annotation_files accepts gz, but faidx typically expects plain.
    # We only ensure .fai when fasta is plain.
    return 0
  fi
  if [[ -f "${fasta}.fai" ]]; then
    return 0
  fi
  echo "[INFO] Missing FASTA .fai; creating with samtools faidx"
  singularity exec \
    --bind "$WORKDIR:$WORKDIR${BIND_EXTRA:+,$BIND_EXTRA}" \
    --pwd "$WORKDIR" \
    "$img" \
    samtools faidx "$fasta"
}

IMG="$(pull_img "$IMG_URL")"
SAMTOOLS_IMG="$(pull_img "$SAMTOOLS_IMG_URL")"

ensure_fai "$FASTA" "$SAMTOOLS_IMG"

singularity exec \
  --bind "$WORKDIR:$WORKDIR${BIND_EXTRA:+,$BIND_EXTRA}" \
  --pwd "$WORKDIR" \
  "$IMG" \
  bash -lc "
set -euo pipefail
cd '$OUTDIR'

cat <<'RSCRIPT' > script.R
library(RiboseQC)

prepare_annotation_files(
  annotation_directory = \".\",
  genome_seq = \"$FASTA\",
  gtf_file = \"$GTF\",
  scientific_name = \"Genome.annotation\",
  annotation_name = \"custom\",
  export_bed_tables_TxDb = FALSE,
  forge_BSgenome = FALSE,
  create_TxDb = TRUE
)

writeLines(
  c(
    'RIBOSEQC_PREPAREANNOTATION:',
    paste0('    riboseqc: \"', packageVersion('RiboseQC'), '\"')
  ),
  'versions.yml'
)
RSCRIPT

Rscript script.R
"

echo "[OK] Annotation output in: $OUTDIR (look for *_Rannot)"
