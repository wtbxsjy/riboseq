#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  retrieve_ensembl_data.sh --species SPECIES --release NUM --assembly STR [OPTIONS]

Required:
  --species STR        Ensembl species name (e.g., homo_sapiens, mus_musculus)
  --release NUM        Ensembl release number (e.g., 110)
  --assembly STR       Genome assembly (e.g., GRCh38, GRCm39)

Options:
  --outdir DIR         Output directory (default: ./Ens<release>)
  --version-suffix STR File version suffix used in Ensembl GTF filenames
                       (default: <release>; e.g. use 110.38 for GRCm38)
  --base-url URL       Base URL for Ensembl FTP
                       (default: https://ftp.ensembl.org/pub/release-<release>)
  --keep-downloads     Keep downloaded .gz files (default: delete after extraction)
  -h, --help           Show this help message

Description:
  Downloads and prepares Ensembl annotation files for gencode-riboseqORFs.
  Creates standardized filenames in the output directory:
    PROTEOME_FASTA
    TRANSCRIPTOME_FASTA
    SORTED_TRANSCRIPTOME_GTF
    TRANSCRIPT_SUPPORT
    PSITES_BED
EOF
}

SPECIES=""
RELEASE=""
ASSEMBLY=""
OUTDIR=""
VERSION_SUFFIX=""
BASE_URL=""
KEEP_DOWNLOADS=false

while [[ $# -gt 0 ]]; do
  case "$1" in
    --species) SPECIES="$2"; shift 2;;
    --release) RELEASE="$2"; shift 2;;
    --assembly) ASSEMBLY="$2"; shift 2;;
    --outdir) OUTDIR="$2"; shift 2;;
    --version-suffix) VERSION_SUFFIX="$2"; shift 2;;
    --base-url) BASE_URL="$2"; shift 2;;
    --keep-downloads) KEEP_DOWNLOADS=true; shift;;
    -h|--help) usage; exit 0;;
    *) echo "[ERROR] Unknown argument: $1"; usage; exit 2;;
  esac
done

if [[ -z "$SPECIES" || -z "$RELEASE" || -z "$ASSEMBLY" ]]; then
  echo "[ERROR] --species, --release, and --assembly are required."
  usage
  exit 2
fi

if [[ -z "$OUTDIR" ]]; then
  OUTDIR="./Ens${RELEASE}"
fi

if [[ -z "$VERSION_SUFFIX" ]]; then
  VERSION_SUFFIX="${RELEASE}"
fi

species_lower="$(echo "$SPECIES" | tr '[:upper:]' '[:lower:]' | tr ' ' '_' | tr -s '_')"
species_caps="$(echo "$species_lower" | awk -F_ 'BEGIN{OFS="_"} {$1=toupper(substr($1,1,1)) substr($1,2); $1=$1; print}')"

if [[ -z "$BASE_URL" ]]; then
  BASE_URL="https://ftp.ensembl.org/pub/release-${RELEASE}"
fi
BASE_URL="${BASE_URL%/}"

mkdir -p "$OUTDIR"
OUTDIR="$(cd "$OUTDIR" && pwd)"
cd "$OUTDIR"

download_file() {
  local url="$1"
  local out="$2"
  local required="${3:-true}"
  local urls=("$url")
  if [[ "$url" == https://* ]]; then
    urls+=("http://${url#https://}")
    urls+=("ftp://${url#https://}")
  elif [[ "$url" == http://* ]]; then
    urls+=("https://${url#http://}")
    urls+=("ftp://${url#http://}")
  elif [[ "$url" == ftp://* ]]; then
    urls+=("https://${url#ftp://}")
    urls+=("http://${url#ftp://}")
  fi

  for try_url in "${urls[@]}"; do
    if wget -q -c --no-check-certificate "$try_url" -O "$out" && [[ -s "$out" ]]; then
      return 0
    fi
  done

  rm -f "$out"
  if [[ "$required" == "true" ]]; then
    echo "[ERROR] Failed to download from Ensembl FTP: ${urls[0]}"
    exit 1
  fi
  echo "[WARN] Optional file not found or download failed: ${urls[0]}"
  return 1
}

strip_fasta_versions() {
  awk 'BEGIN{FS="."} /^>/{print $1; next} {print}' "$1"
}

GTF_FILE="${species_caps}.${ASSEMBLY}.${VERSION_SUFFIX}.gtf.gz"
CDNA_FILE="${species_caps}.${ASSEMBLY}.cdna.all.fa.gz"
NCRNA_FILE="${species_caps}.${ASSEMBLY}.ncrna.fa.gz"
PEP_FILE="${species_caps}.${ASSEMBLY}.pep.all.fa.gz"

download_file "${BASE_URL}/gtf/${species_lower}/${GTF_FILE}" "$GTF_FILE" true
download_file "${BASE_URL}/fasta/${species_lower}/cdna/${CDNA_FILE}" "$CDNA_FILE" true
NCRNA_AVAILABLE=true
if ! download_file "${BASE_URL}/fasta/${species_lower}/ncrna/${NCRNA_FILE}" "$NCRNA_FILE" false; then
  NCRNA_AVAILABLE=false
fi
download_file "${BASE_URL}/fasta/${species_lower}/pep/${PEP_FILE}" "$PEP_FILE" true

gunzip -c "$GTF_FILE" > "${GTF_FILE%.gz}"
gunzip -c "$CDNA_FILE" > "${CDNA_FILE%.gz}"
if [[ "$NCRNA_AVAILABLE" == "true" ]]; then
  gunzip -c "$NCRNA_FILE" > "${NCRNA_FILE%.gz}"
fi
gunzip -c "$PEP_FILE" > "${PEP_FILE%.gz}"

GTF_SORTED="${species_caps}.${ASSEMBLY}.${VERSION_SUFFIX}.sorted.gtf"
(grep "^#" "${GTF_FILE%.gz}" || true; grep -v "^#" "${GTF_FILE%.gz}" | sort -k1,1 -k4,4n -k5,5n) > "$GTF_SORTED"

TRANSCRIPTOME="${species_caps}.${ASSEMBLY}.transcriptome.fa"
if [[ "$NCRNA_AVAILABLE" == "true" ]]; then
  cat <(strip_fasta_versions "${CDNA_FILE%.gz}") <(strip_fasta_versions "${NCRNA_FILE%.gz}") > "$TRANSCRIPTOME"
else
  strip_fasta_versions "${CDNA_FILE%.gz}" > "$TRANSCRIPTOME"
fi

tmp_pep="${PEP_FILE%.gz}.tmp"
strip_fasta_versions "${PEP_FILE%.gz}" > "$tmp_pep"
mv "$tmp_pep" "${PEP_FILE%.gz}"

PSITES_BED="psites.bed"
python3 - <<'PYTHON_SCRIPT' "$GTF_SORTED" "$PSITES_BED"
import sys

gtf_file = sys.argv[1]
output_bed = sys.argv[2]

with open(gtf_file, 'r') as gtf, open(output_bed, 'w') as bed:
    for line in gtf:
        if line.startswith('#'):
            continue
        fields = line.rstrip().split('\t')
        if len(fields) < 9 or fields[2] != 'start_codon':
            continue
        chrom = fields[0]
        start = int(fields[3]) - 1  # BED is 0-based
        end = int(fields[4])
        strand = fields[6]
        attrs = fields[8]
        transcript_id = ""
        for attr in attrs.split(';'):
            attr = attr.strip()
            if attr.startswith('transcript_id'):
                transcript_id = attr.split('"')[1]
                break
        if transcript_id:
            bed.write(f"{chrom}\t{start}\t{end}\t{transcript_id}\t.\t{strand}\n")
PYTHON_SCRIPT

TSL_FILE="transcript_support_level.txt"
BIOMART_URL="http://ensembl.org/biomart/martservice?query=<Query virtualSchemaName=\"default\" formatter=\"TSV\" header=\"1\" uniqueRows=\"0\" count=\"\" datasetConfigVersion=\"0.6\"><Dataset name=\"${species_lower}_gene_ensembl\" interface=\"default\"><Attribute name=\"ensembl_transcript_id\"/><Attribute name=\"transcript_tsl\"/><Attribute name=\"transcript_appris\"/></Dataset></Query>"
TSL_OK=false
if wget -q -O "${TSL_FILE}.tmp" "$BIOMART_URL" && \
   [[ -s "${TSL_FILE}.tmp" ]] && \
   head -1 "${TSL_FILE}.tmp" | grep -qiE "transcript|ensembl_transcript_id" && \
   ! head -1 "${TSL_FILE}.tmp" | grep -qiE "^Query ERROR|^Error|Exception"; then
  mv "${TSL_FILE}.tmp" "$TSL_FILE"
  TSL_OK=true
fi
rm -f "${TSL_FILE}.tmp"
if [[ "$TSL_OK" == "false" ]]; then
  echo "[WARN] BioMart download failed or returned invalid data; generating basic transcript support from GTF."
  python3 - <<'PYTHON_SCRIPT' "$GTF_SORTED" "$TSL_FILE"
import re
import sys

gtf_file = sys.argv[1]
output_file = sys.argv[2]
transcripts = {}

with open(gtf_file, 'r') as gtf:
    for line in gtf:
        if line.startswith('#'):
            continue
        fields = line.rstrip().split('\t')
        if len(fields) < 9 or fields[2] != 'transcript':
            continue
        attrs = fields[8]
        transcript_id = ''
        tsl = 'NA'
        appris = 'NA'
        for attr in attrs.split(';'):
            attr = attr.strip()
            if attr.startswith('transcript_id'):
                transcript_id = attr.split('"')[1]
            elif 'transcript_support_level' in attr or 'tsl' in attr.lower():
                match = re.search(r'(\d+)', attr)
                if match:
                    tsl = match.group(1)
        if transcript_id:
            transcripts[transcript_id] = {'tsl': tsl, 'appris': appris}

with open(output_file, 'w') as out:
    out.write("Transcript stable ID\tTSL\tAPPRIS\n")
    for tid, data in transcripts.items():
        out.write(f"{tid}\t{data['tsl']}\t{data['appris']}\n")
PYTHON_SCRIPT
fi

ln -sf "$GTF_SORTED" "SORTED_TRANSCRIPTOME_GTF"
ln -sf "$TRANSCRIPTOME" "TRANSCRIPTOME_FASTA"
ln -sf "${PEP_FILE%.gz}" "PROTEOME_FASTA"
ln -sf "$TSL_FILE" "TRANSCRIPT_SUPPORT"
ln -sf "$PSITES_BED" "PSITES_BED"

if [[ "$KEEP_DOWNLOADS" == "false" ]]; then
  rm -f *.gz
fi

echo "[OK] Ensembl annotation prepared at: $OUTDIR"
ls -lh SORTED_TRANSCRIPTOME_GTF TRANSCRIPTOME_FASTA PROTEOME_FASTA TRANSCRIPT_SUPPORT PSITES_BED
