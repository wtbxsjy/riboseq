#!/usr/bin/env bash
set -euo pipefail

# Mirrors: modules/local/sorf_bam_filter/main.nf
IMG_URL="https://depot.galaxyproject.org/singularity/samtools:1.21--h50ea8bc_0"

usage() {
  cat <<'EOF'
Usage:
  01_sorf_bam_filter.sh --sample ID --bam in.bam --fai genome.fa.fai [options]

Required:
  --sample   Sample ID (prefix)
  --bam      Input BAM (sorted recommended)
  --fai      FASTA index (.fai) for contig list

Options:
  --unique-mode   auto|nh|mapq   (default: auto)
  --mapq          MAPQ threshold (default: 60)
  --len-min       read length min (SEQ length) (default: 28)
  --len-max       read length max (SEQ length) (default: 30)
  --exclude-regex contig regex to EXCLUDE (default from pipeline in nextflow.config is not auto-applied here)
  --cpus          threads for samtools (index/view) (default: 4)
  --outdir        output directory (default: ./out_sorf_filter)

Env:
  BIND_EXTRA   extra singularity binds, comma-separated (e.g. /mnt:/mnt)
EOF
}

SAMPLE=""
BAM=""
FAI=""
UNIQUE_MODE="auto"
MAPQ=60
LEN_MIN=28
LEN_MAX=30
EXCLUDE_REGEX=""
OUTDIR="./out_sorf_filter"
CPUS=4

while [[ $# -gt 0 ]]; do
  case "$1" in
    --sample) SAMPLE="$2"; shift 2;;
    --bam) BAM="$2"; shift 2;;
    --fai) FAI="$2"; shift 2;;
    --unique-mode) UNIQUE_MODE="$2"; shift 2;;
    --mapq) MAPQ="$2"; shift 2;;
    --len-min) LEN_MIN="$2"; shift 2;;
    --len-max) LEN_MAX="$2"; shift 2;;
    --exclude-regex) EXCLUDE_REGEX="$2"; shift 2;;
    --cpus) CPUS="$2"; shift 2;;
    --outdir) OUTDIR="$2"; shift 2;;
    -h|--help) usage; exit 0;;
    *) echo "Unknown arg: $1"; usage; exit 2;;
  esac
done

if [[ -z "$SAMPLE" || -z "$BAM" || -z "$FAI" ]]; then
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
    samtools index -@ "$CPUS" "$bam"
}

IMG="$(pull_img "$IMG_URL")"

# Make sure BAM is indexable
ensure_bai "$BAM" "$IMG"

# Run filter inside container, matching the pipeline logic
singularity exec \
  --bind "$WORKDIR:$WORKDIR${BIND_EXTRA:+,$BIND_EXTRA}" \
  --pwd "$WORKDIR" \
  "$IMG" \
  bash -lc "
set -euo pipefail
cd '$OUTDIR'

prefix='$SAMPLE'
mode='$UNIQUE_MODE'
mapq=$MAPQ
rlmin=$LEN_MIN
rlmax=$LEN_MAX
re='$EXCLUDE_REGEX'

if [[ -n \"\$re\" ]]; then
  awk -v re=\"\$re\" 'BEGIN{FS=\"\\t\"} \$1 ~ re {print \$1}' '$FAI' | sort -u > \"\${prefix}.sorf.excluded_contigs.txt\" || true
fi

# Exclude: unmapped(0x4), secondary(0x100), duplicate(0x400), supplementary(0x800)
# 0xD04 = 0x4 + 0x100 + 0x400 + 0x800

total_primary_mapped=\`samtools view -@ '$CPUS' -c -F 0xD04 '$BAM'\`

samtools view -h -F 0xD04 '$BAM' \
  | awk -v mode=\"\$mode\" -v mapq=\"\$mapq\" -v rlmin=\"\$rlmin\" -v rlmax=\"\$rlmax\" -v re=\"\$re\" '
      BEGIN{OFS=\"\\t\"}
      /^@/ {print; next}
      {
        rname=\$3
        if (re != \"\" && rname ~ re) next

        seqlen=length(\$10)
        if (rlmin > 0 && seqlen < rlmin) next
        if (rlmax > 0 && seqlen > rlmax) next

        if (mode == \"mapq\") {
          if (\$5 < mapq) next
        } else if (mode == \"nh\" || mode == \"auto\") {
          nh=\"\"
          for (i=12; i<=NF; i++) {
            if (\$i ~ /^NH:i:/) { split(\$i,a,\":\"); nh=a[3]; break }
          }
          if (mode == \"nh\") {
            if (nh == \"\" || nh != 1) next
          } else {
            if (nh != \"\") {
              if (nh != 1) next
            } else {
              if (\$5 < mapq) next
            }
          }
        }

        print
      }
    ' \
      | samtools view -@ '$CPUS' -b -o \"\${prefix}.sorf.filtered.bam\" -

    kept_primary_mapped=\`samtools view -@ '$CPUS' -c -F 0xD04 \"\${prefix}.sorf.filtered.bam\"\`

printf \"sample\\ttotal_primary_mapped\\tkept_primary_mapped\\tpct_kept\\n\" > \"\${prefix}.sorf.filter_stats.tsv\"
awk -v s=\"\${prefix}\" -v t=\"\${total_primary_mapped}\" -v k=\"\${kept_primary_mapped}\" 'BEGIN{pct=(t>0)?(100.0*k/t):0; printf \"%s\\t%d\\t%d\\t%.2f\\n\", s, t, k, pct}' >> \"\${prefix}.sorf.filter_stats.tsv\"

samtools --version | head -n1 > versions.samtools.txt
"

echo "[OK] Output: $OUTDIR/${SAMPLE}.sorf.filtered.bam"
