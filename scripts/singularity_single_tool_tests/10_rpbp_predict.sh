#!/usr/bin/env bash
set -euo pipefail

# Mirrors: modules/local/rpbp/predict/main.nf
IMG_URL="https://depot.galaxyproject.org/singularity/rpbp:4.0.1--py312hf731ba3_0"
SAMTOOLS_IMG_URL="https://depot.galaxyproject.org/singularity/samtools:1.21--h50ea8bc_0"

usage() {
  cat <<'EOF'
Usage:
  10_rpbp_predict.sh --sample ID --bam in.bam --orfs-genomic genome.orfs-genomic.bed.gz --orfs-exons genome.orfs-exons.bed.gz [--outdir DIR] [--cpus N]

Required:
  --sample       sample ID (prefix)
  --bam          input BAM (filtered BAM recommended)
  --orfs-genomic rp-bp orfs genomic bed.gz
  --orfs-exons   rp-bp orfs exons bed.gz

Options:
  --outdir output directory (default: ./out_rpbp_predict)
  --cpus   threads (default: 8)
  --args   extra args passed to rp-bp steps that support it (quoted string)

Env:
  BIND_EXTRA extra singularity binds, comma-separated (e.g. /mnt:/mnt)
EOF
}

SAMPLE=""
BAM=""
ORFS_GEN=""
ORFS_EX=""
OUTDIR="./out_rpbp_predict"
CPUS=8
EXTRA_ARGS=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --sample) SAMPLE="$2"; shift 2;;
    --bam) BAM="$2"; shift 2;;
    --orfs-genomic) ORFS_GEN="$2"; shift 2;;
    --orfs-exons) ORFS_EX="$2"; shift 2;;
    --outdir) OUTDIR="$2"; shift 2;;
    --cpus) CPUS="$2"; shift 2;;
    --args) EXTRA_ARGS="$2"; shift 2;;
    -h|--help) usage; exit 0;;
    *) echo "Unknown arg: $1"; usage; exit 2;;
  esac
done

if [[ -z "$SAMPLE" || -z "$BAM" || -z "$ORFS_GEN" || -z "$ORFS_EX" ]]; then
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

ensure_bai() {
  local bam="$1"
  local img="$2"
  if [[ -f "${bam}.bai" || -f "${bam%.bam}.bai" || -f "${bam}.csi" ]]; then
    return 0
  fi
  echo "[INFO] Missing BAM index; creating with samtools index"
  singularity exec \
    --bind "$WORKDIR:$WORKDIR${BIND_EXTRA:+,$BIND_EXTRA}" \
    --pwd "$WORKDIR" \
    "$img" \
    samtools index -@ 2 "$bam"
}

IMG="$(pull_img "$IMG_URL")"
SAMTOOLS_IMG="$(pull_img "$SAMTOOLS_IMG_URL")"

ensure_bai "$BAM" "$SAMTOOLS_IMG"

singularity exec \
  --bind "$WORKDIR:$WORKDIR${BIND_EXTRA:+,$BIND_EXTRA}" \
  --pwd "$WORKDIR" \
  "$IMG" \
  bash -lc "
set -euo pipefail
cd '$OUTDIR'

# Find models (same as module)
python3 -c '
import os
import rpbp
import sys

pkg_path = os.path.dirname(rpbp.__file__)
models_path = os.path.join(pkg_path, "models")

def get_models(subdir):
    p = os.path.join(models_path, subdir)
    if not os.path.exists(p):
        p = os.path.join(sys.prefix, "share", "rpbp", "models", subdir)

    if not os.path.exists(p):
        sys.stderr.write(f"WARNING: Could not find models directory for {subdir} at {p}\\n")
        return ""

    files = []
    for f in os.listdir(p):
        fp = os.path.join(p, f)
        if os.path.isfile(fp) and not f.endswith(".stan") and not f.endswith(".py") and not f.endswith(".pyc") and not f.startswith("__"):
            files.append(fp)

    if not files:
        sys.stderr.write(f"WARNING: No models found in {p}\\n")

    return " ".join(files)

open("periodic_models.txt","w").write(get_models("periodic"))
open("nonperiodic_models.txt","w").write(get_models("nonperiodic"))
open("translated_models.txt","w").write(get_models("translated"))
open("untranslated_models.txt","w").write(get_models("untranslated"))
'

PERIODIC_MODELS=\$(cat periodic_models.txt)
NONPERIODIC_MODELS=\$(cat nonperiodic_models.txt)
TRANSLATED_MODELS=\$(cat translated_models.txt)
UNTRANSLATED_MODELS=\$(cat untranslated_models.txt)

extract-metagene-profiles \
  '$BAM' \
  '$ORFS_GEN' \
  ${SAMPLE}.metagene-profiles.csv.gz \
  --num-cpus $CPUS

estimate-metagene-profile-bayes-factors \
  ${SAMPLE}.metagene-profiles.csv.gz \
  ${SAMPLE}.metagene-profile-bayes-factors.csv.gz \
  --periodic-models \$PERIODIC_MODELS \
  --nonperiodic-models \$NONPERIODIC_MODELS \
  --num-cpus $CPUS \
  $EXTRA_ARGS

select-periodic-offsets \
  ${SAMPLE}.metagene-profile-bayes-factors.csv.gz \
  ${SAMPLE}.periodic-offsets.csv.gz

ARGS=\$(python3 -c "import pandas as pd; df=pd.read_csv('${SAMPLE}.periodic-offsets.csv.gz'); print('--lengths ' + ' '.join(map(str, df['length'].astype(int))) + ' --offsets ' + ' '.join(map(str, df['highest_peak_offset'].astype(int))))")

extract-orf-profiles \
  '$BAM' \
  '$ORFS_GEN' \
  '$ORFS_EX' \
  ${SAMPLE}.profiles.mtx.gz \
  \$ARGS \
  --num-cpus $CPUS

estimate-orf-bayes-factors \
  ${SAMPLE}.profiles.mtx.gz \
  '$ORFS_GEN' \
  ${SAMPLE}.bayes-factors.bed.gz \
  --translated-models \$TRANSLATED_MODELS \
  --untranslated-models \$UNTRANSLATED_MODELS \
  --num-cpus $CPUS \
  $EXTRA_ARGS

select-final-prediction-set \
  ${SAMPLE}.bayes-factors.bed.gz \
  '$ORFS_GEN' \
  ${SAMPLE}.predicted-orfs.bed.gz \
  ${SAMPLE}.predicted-orfs.dna.fa \
  ${SAMPLE}.predicted-orfs.protein.fa

python3 -c 'import rpbp; print(rpbp.__version__)' > versions.rpbp.txt
"

echo "[OK] Output: $OUTDIR/${SAMPLE}.predicted-orfs.bed.gz"
