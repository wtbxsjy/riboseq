#!/usr/bin/env bash
set -euo pipefail

# Mirrors: modules/local/extract_rl_cutoff/main.nf
# Extracts optimal P-site offsets from RiboseQC P_sites_calcs output
# Selects rows with max_coverage=TRUE, extracts read_length, cutoff, comp

# Uses the same container as RiboseQC for R environment
IMG_URL="https://depot.galaxyproject.org/singularity/riboseqc:1.1--r36_1"

usage() {
  cat <<'EOF'
Usage:
  19_extract_rl_cutoff.sh --sample ID --psites-calcs FILE [OPTIONS]

Required:
  --sample         Sample ID (prefix)
  --psites-calcs   RiboseQC P_sites_calcs file (from 03_riboseqc_analysis.sh)

Options:
  --outdir    Output directory (default: ./out_extract_rl_cutoff)

Env:
  BIND_EXTRA  Extra singularity binds (e.g. /mnt:/mnt)

Output:
  {sample}_rl_cutoff.tsv  - Tab-separated file with columns:
                            read_length, cutoff, comp
EOF
}

SAMPLE=""
PSITES_CALCS=""
OUTDIR="./out_extract_rl_cutoff"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --sample) SAMPLE="$2"; shift 2;;
    --psites-calcs) PSITES_CALCS="$2"; shift 2;;
    --outdir) OUTDIR="$2"; shift 2;;
    -h|--help) usage; exit 0;;
    *) echo "Unknown arg: $1"; usage; exit 2;;
  esac
done

if [[ -z "$SAMPLE" || -z "$PSITES_CALCS" ]]; then
  usage
  exit 2
fi

# Create output and container directories
mkdir -p "$OUTDIR" "./containers"
OUTDIR="$(cd "$OUTDIR" && pwd)"
WORKDIR="$(pwd)"

# Get absolute path
abspath() {
  python3 -c "import os,sys; print(os.path.abspath(sys.argv[1]))" "$1"
}

PSITES_CALCS="$(abspath "$PSITES_CALCS")"

# Build bind mounts
BIND_SPEC="$WORKDIR:$WORKDIR,$OUTDIR:$OUTDIR,$(dirname "$PSITES_CALCS"):$(dirname "$PSITES_CALCS")"

# Pull container
pull_img() {
  local url="$1"
  local base="$(basename "$url" | sed 's/:/_/g')"
  local sif="$(pwd)/containers/${base}.sif"
  if [[ ! -f "$sif" ]]; then
    echo "[INFO] Pulling $url..."
    singularity pull --disable-cache --force "$sif" "$url"
  fi
  echo "$sif"
}

IMG="$(pull_img "$IMG_URL")"

echo "[INFO] Extracting P-site offset from: $PSITES_CALCS"
echo "[INFO] Output: $OUTDIR/${SAMPLE}_rl_cutoff.tsv"

singularity exec \
  --bind "$BIND_SPEC${BIND_EXTRA:+,$BIND_EXTRA}" \
  --pwd "$WORKDIR" \
  "$IMG" \
  bash -lc "
set -euo pipefail
cd '$OUTDIR'

cat <<'RSCRIPT' > script.R
# Extract P-site offsets from RiboseQC P_sites_calcs
# Selects rows with max_coverage=TRUE (or 1 for numeric)

psites_calcs <- readRDS('$PSITES_CALCS')

# Debug output
cat('Data structure:\\n')
cat('Rows:', nrow(psites_calcs), '\\n')
cat('Columns:', paste(colnames(psites_calcs), collapse=', '), '\\n')

if ('max_coverage' %in% colnames(psites_calcs)) {
    cat('max_coverage class:', class(psites_calcs\$max_coverage), '\\n')
    cat('max_coverage unique values:', paste(unique(psites_calcs\$max_coverage), collapse=', '), '\\n')
    
    # Handle both logical and numeric max_coverage
    if (is.logical(psites_calcs\$max_coverage)) {
        optimal <- psites_calcs[psites_calcs\$max_coverage == TRUE, ]
    } else if (is.numeric(psites_calcs\$max_coverage)) {
        optimal <- psites_calcs[psites_calcs\$max_coverage == 1, ]
    } else {
        optimal <- psites_calcs[as.logical(psites_calcs\$max_coverage), ]
    }
    
    if (nrow(optimal) == 0) {
        cat('WARNING: No rows with max_coverage indicator. Using all unique read lengths.\\n')
        optimal <- psites_calcs[!duplicated(psites_calcs\$read_length), ]
    }
    
    # Select relevant columns
    result <- optimal[, c('read_length', 'cutoff', 'comp')]
    
    # Write output
    write.table(result, '${SAMPLE}_rl_cutoff.tsv', 
                sep='\\t', row.names=FALSE, quote=FALSE)
    
    cat('Extracted', nrow(result), 'P-site offset entries\\n')
    cat('Read length range:', min(result\$read_length), '-', max(result\$read_length), '\\n')
} else {
    cat('ERROR: max_coverage column not found in P_sites_calcs\\n')
    quit(status=1)
}
RSCRIPT

Rscript script.R
"

echo "[OK] P-site offsets extracted to: $OUTDIR/${SAMPLE}_rl_cutoff.tsv"
cat "$OUTDIR/${SAMPLE}_rl_cutoff.tsv"
