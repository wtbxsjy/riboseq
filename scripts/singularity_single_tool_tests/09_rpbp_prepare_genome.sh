#!/usr/bin/env bash
set -euo pipefail

# Mirrors: modules/local/rpbp/prepare_genome/main.nf
IMG_URL="https://depot.galaxyproject.org/singularity/rpbp:4.0.1--py312hf731ba3_0"

usage() {
  cat <<'EOF'
Usage:
  09_rpbp_prepare_genome.sh --fasta genome.fa --gtf annot.gtf --rrna ribosomal.fa [--outdir DIR] [--cpus N]

Required:
  --fasta genome FASTA
  --gtf   GTF
  --rrna  ribosomal FASTA

Options:
  --outdir output directory (default: ./out_rpbp_genome)
  --cpus   threads (default: 8)
  --args   extra args passed to prepare-rpbp-genome (quoted string)

Env:
  BIND_EXTRA extra singularity binds, comma-separated (e.g. /mnt:/mnt)
EOF
}

FASTA=""
GTF=""
RRNA=""
OUTDIR="./out_rpbp_genome"
CPUS=8
EXTRA_ARGS=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --fasta) FASTA="$2"; shift 2;;
    --gtf) GTF="$2"; shift 2;;
    --rrna) RRNA="$2"; shift 2;;
    --outdir) OUTDIR="$2"; shift 2;;
    --cpus) CPUS="$2"; shift 2;;
    --args) EXTRA_ARGS="$2"; shift 2;;
    -h|--help) usage; exit 0;;
    *) echo "Unknown arg: $1"; usage; exit 2;;
  esac
done

if [[ -z "$FASTA" || -z "$GTF" || -z "$RRNA" ]]; then
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

mkdir -p ribosomal_index star_index

cat <<EOF > config.yaml
genome_base_path: .
genome_name: genome
gtf: $GTF
fasta: $FASTA
ribosomal_fasta: $RRNA
ribosomal_index: ./ribosomal_index/rRNA
star_index: ./star_index
EOF

prepare-rpbp-genome \
  config.yaml \
  --num-cpus $CPUS \
  --overwrite \
  $EXTRA_ARGS

python3 -c 'import rpbp; print(rpbp.__version__)' > versions.rpbp.txt
"

echo "[OK] Outputs: $OUTDIR/transcript-index/genome.orfs-genomic.bed.gz and *.orfs-exons.bed.gz"
