#!/usr/bin/env bash
set -euo pipefail

# Mirrors: modules/local/ribocode/detect/main.nf
IMG_URL="https://depot.galaxyproject.org/singularity/ribocode:1.2.11--pyh145b6a8_1"
SAMTOOLS_IMG_URL="https://depot.galaxyproject.org/singularity/samtools:1.21--h50ea8bc_0"

usage() {
  cat <<'EOF'
Usage:
  11_ribocode_detect.sh --sample ID --bam in.bam --gtf annot.gtf --fasta genome.fa [--outdir DIR] [--stranded forward|reverse|unstranded] [--cpus N]

Notes:
  - 更推荐 transcriptome BAM（pipeline 在 FASTQ 模式会额外生成 transcriptome BAM）。
  - 低深度/测试数据可能因 periodicity 不足失败，这是常见现象。

Required:
  --sample sample ID (prefix)
  --bam    input BAM
  --gtf    GTF (can be .gz)
  --fasta  genome FASTA

Options:
  --outdir   output directory (default: ./out_ribocode)
  --stranded forward|reverse|unstranded (default: forward)
  --cpus     threads for samtools index (default: 4)
  --args     extra args passed to RiboCode_onestep (quoted string)

Env:
  BIND_EXTRA extra singularity binds, comma-separated (e.g. /mnt:/mnt)
EOF
}

SAMPLE=""
BAM=""
GTF=""
FASTA=""
OUTDIR="./out_ribocode"
STRANDED="forward"
EXTRA_ARGS=""
CPUS=4

while [[ $# -gt 0 ]]; do
  case "$1" in
    --sample) SAMPLE="$2"; shift 2;;
    --bam) BAM="$2"; shift 2;;
    --gtf) GTF="$2"; shift 2;;
    --fasta) FASTA="$2"; shift 2;;
    --outdir) OUTDIR="$2"; shift 2;;
    --stranded) STRANDED="$2"; shift 2;;
    --cpus) CPUS="$2"; shift 2;;
    --args) EXTRA_ARGS="$2"; shift 2;;
    -h|--help) usage; exit 0;;
    *) echo "Unknown arg: $1"; usage; exit 2;;
  esac
done

if [[ -z "$SAMPLE" || -z "$BAM" || -z "$GTF" || -z "$FASTA" ]]; then
  usage
  exit 2
fi

mkdir -p "$OUTDIR" "./containers"
OUTDIR="$(cd "$OUTDIR" && pwd)"
WORKDIR="$(pwd)"

# Auto-detect and bind mount input file directories
auto_bind_paths() {
  local bam_path="$1"
  local gtf_path="$2"
  local fasta_path="$3"
  local bind_paths=""

  # Convert to absolute paths
  bam_abs="$(cd "$(dirname "$bam_path")" && pwd)/$(basename "$bam_path")"
  gtf_abs="$(cd "$(dirname "$gtf_path")" && pwd)/$(basename "$gtf_path")"
  fasta_abs="$(cd "$(dirname "$fasta_path")" && pwd)/$(basename "$fasta_path")"

  # Extract parent directories
  bam_dir="$(dirname "$bam_abs")"
  gtf_dir="$(dirname "$gtf_abs")"
  fasta_dir="$(dirname "$fasta_abs")"

  # Add unique directories to bind list
  for dir in "$bam_dir" "$gtf_dir" "$fasta_dir" "$OUTDIR"; do
    if [[ ":$bind_paths:" != *":$dir:"* ]]; then
      bind_paths="${bind_paths:+$bind_paths,}$dir:$dir"
    fi
  done

  echo "$bind_paths"
}

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
  local binds="$3"
  if [[ -f "${bam}.bai" || -f "${bam%.bam}.bai" || -f "${bam}.csi" ]]; then
    return 0
  fi
  echo "[INFO] Missing BAM index; creating with samtools index"
  singularity exec \
    --bind "$binds" \
    --pwd "$WORKDIR" \
    "$img" \
    samtools index -@ "$CPUS" "$bam"
}

# Auto-detect bind mounts
AUTO_BINDS="$(auto_bind_paths "$BAM" "$GTF" "$FASTA")"
echo "[INFO] Auto-detected bind mounts: $AUTO_BINDS"

# Convert inputs to absolute paths for container
BAM="$(cd "$(dirname "$BAM")" && pwd)/$(basename "$BAM")"
GTF="$(cd "$(dirname "$GTF")" && pwd)/$(basename "$GTF")"
FASTA="$(cd "$(dirname "$FASTA")" && pwd)/$(basename "$FASTA")"

IMG="$(pull_img "$IMG_URL")"
SAMTOOLS_IMG="$(pull_img "$SAMTOOLS_IMG_URL")"

# Combine auto-detected binds with BIND_EXTRA
ALL_BINDS="$WORKDIR:$WORKDIR,$AUTO_BINDS${BIND_EXTRA:+,$BIND_EXTRA}"

echo "[INFO] Final bind mounts: $ALL_BINDS"

ensure_bai "$BAM" "$SAMTOOLS_IMG" "$ALL_BINDS"

STR_OPT="yes"
case "$STRANDED" in
  forward) STR_OPT="yes";;
  reverse) STR_OPT="reverse";;
  unstranded) STR_OPT="no";;
  *) echo "Invalid --stranded: $STRANDED"; exit 2;;
esac

singularity exec \
  --bind "$ALL_BINDS" \
  --pwd "$WORKDIR" \
  "$IMG" \
  bash -lc "
set -euo pipefail
cd '$OUTDIR'

export TMPDIR=\"$OUTDIR/tmp\"
mkdir -p \"\$TMPDIR\"
export MPLCONFIGDIR=\"$OUTDIR/mplconfig\"
mkdir -p \"\$MPLCONFIGDIR\"

GTF_IN='$GTF'
if [[ \"\$GTF_IN\" == *.gz ]]; then
  gunzip -c \"\$GTF_IN\" > input.gtf
  GTF_IN=input.gtf
fi

python - \"\$GTF_IN\" ribocode.gtf <<'PY'
import re
import sys

gtf_in = sys.argv[1]
gtf_out = sys.argv[2]

level_re = re.compile(r'(?:^|;\\s*)level\\s+\"[^\"]*\"\\s*;')

with open(gtf_in, 'rt', encoding='utf-8', errors='replace') as fin, open(gtf_out, 'wt', encoding='utf-8') as fout:
    for line in fin:
        if line.startswith('#') or not line.strip():
            fout.write(line)
            continue
        parts = line.rstrip('\\n').split('\\t')
        if len(parts) < 9:
            fout.write(line)
            continue
        attrs = parts[8].strip()
        if not level_re.search(attrs):
            if attrs and not attrs.endswith(';'):
                attrs += ';'
            attrs += ' level \"NA\";'
            parts[8] = attrs
        fout.write('\\t'.join(parts) + '\\n')
PY

RiboCode_onestep \
  -g ribocode.gtf \
  -f '$FASTA' \
  -r '$BAM' \
  --stranded '$STR_OPT' \
  -o '$SAMPLE' \
  $EXTRA_ARGS \
  || {
    echo \"RiboCode failed for $SAMPLE (likely insufficient periodicity)\" >&2
    echo \"FAILED: No periodicity detected\" > ${SAMPLE}.ribocode_failed.txt
  }

RiboCode_onestep --version 2>&1 | tail -n1 > versions.ribocode.txt
"

echo "[OK] Output directory: $OUTDIR (prefix: $SAMPLE)"
