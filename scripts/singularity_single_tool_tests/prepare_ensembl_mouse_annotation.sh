#!/usr/bin/env bash
set -euo pipefail

################################################################################
# Script: prepare_ensembl_mouse_annotation.sh
# Purpose: Download and prepare Ensembl annotation for Mus musculus (mouse)
#          for use with gencode-riboseqORFs
################################################################################

usage() {
  cat <<'EOF'
Usage:
  prepare_ensembl_mouse_annotation.sh [OPTIONS]

Options:
  --release NUM      Ensembl release number (default: 110)
  --assembly STR     Genome assembly (default: GRCm39)
  --outdir DIR       Output directory (default: ./Ensembl_mouse)
  --keep-downloads   Keep downloaded gz files (default: delete after extraction)
  -h, --help         Show this help message

Description:
  Downloads and prepares Ensembl annotation files for Mus musculus (mouse)
  in the format required by gencode-riboseqORFs.

  Creates standardized symlinks:
    PROTEOME_FASTA           -> Mus_musculus.*.pep.all.fa
    TRANSCRIPTOME_FASTA      -> Mus_musculus.*.cdna.all.fa + ncrna.fa
    SORTED_TRANSCRIPTOME_GTF -> Mus_musculus.*.gtf (sorted)
    TRANSCRIPT_SUPPORT       -> transcript_support_level.txt (from BioMart)
    PSITES_BED               -> psites.bed (generated)

Examples:
  # Download latest Ensembl 110 for GRCm39
  bash prepare_ensembl_mouse_annotation.sh --release 110 --assembly GRCm39

  # Download to specific directory
  bash prepare_ensembl_mouse_annotation.sh --outdir /data/ensembl_mouse_110

  # Download Ensembl 109 for GRCm38
  bash prepare_ensembl_mouse_annotation.sh --release 109 --assembly GRCm38

EOF
}

RELEASE=110
ASSEMBLY="GRCm39"
OUTDIR="./Ensembl_mouse"
KEEP_DOWNLOADS=false

while [[ $# -gt 0 ]]; do
  case "$1" in
    --release) RELEASE="$2"; shift 2;;
    --assembly) ASSEMBLY="$2"; shift 2;;
    --outdir) OUTDIR="$2"; shift 2;;
    --keep-downloads) KEEP_DOWNLOADS=true; shift;;
    -h|--help) usage; exit 0;;
    *) echo "[ERROR] Unknown argument: $1"; usage; exit 2;;
  esac
done

# Create output directory
mkdir -p "$OUTDIR"
OUTDIR="$(cd "$OUTDIR" && pwd)"

echo "========================================================================"
echo "  Preparing Ensembl Annotation for Mus musculus (Mouse)"
echo "========================================================================"
echo ""
echo "Release: $RELEASE"
echo "Assembly: $ASSEMBLY"
echo "Output: $OUTDIR"
echo ""

cd "$OUTDIR"

SPECIES="mus_musculus"
BASE_URL="https://ftp.ensembl.org/pub/release-${RELEASE}"

# Determine version suffix based on assembly
# GRCm39 is the current assembly, GRCm38 is older
if [[ "$ASSEMBLY" == "GRCm39" ]]; then
  VERSION_SUFFIX="${RELEASE}"
elif [[ "$ASSEMBLY" == "GRCm38" ]]; then
  VERSION_SUFFIX="${RELEASE}.38"
else
  echo "[WARNING] Unknown assembly: $ASSEMBLY, using default versioning"
  VERSION_SUFFIX="${RELEASE}"
fi

echo "------------------------------------------------------------------------"
echo "Step 1: Downloading GTF annotation"
echo "------------------------------------------------------------------------"
GTF_URL="${BASE_URL}/gtf/${SPECIES}/Mus_musculus.${ASSEMBLY}.${VERSION_SUFFIX}.gtf.gz"
GTF_FILE="Mus_musculus.${ASSEMBLY}.${VERSION_SUFFIX}.gtf.gz"

if [[ ! -f "$GTF_FILE" ]]; then
  echo "Downloading: $GTF_URL"
  wget -c "$GTF_URL" -O "$GTF_FILE"
else
  echo "[INFO] GTF already exists: $GTF_FILE"
fi

# Extract and sort GTF
GTF_EXTRACTED="${GTF_FILE%.gz}"
if [[ ! -f "$GTF_EXTRACTED" ]]; then
  echo "Extracting GTF..."
  gunzip -c "$GTF_FILE" > "$GTF_EXTRACTED"
fi

echo "Sorting GTF by chromosome and position..."
GTF_SORTED="Mus_musculus.${ASSEMBLY}.${VERSION_SUFFIX}.sorted.gtf"
if [[ ! -f "$GTF_SORTED" ]]; then
  (grep "^#" "$GTF_EXTRACTED" || true; grep -v "^#" "$GTF_EXTRACTED" | sort -k1,1 -k4,4n) > "$GTF_SORTED"
fi

echo ""
echo "------------------------------------------------------------------------"
echo "Step 2: Downloading cDNA (transcriptome)"
echo "------------------------------------------------------------------------"
CDNA_URL="${BASE_URL}/fasta/${SPECIES}/cdna/Mus_musculus.${ASSEMBLY}.cdna.all.fa.gz"
CDNA_FILE="Mus_musculus.${ASSEMBLY}.cdna.all.fa.gz"

if [[ ! -f "$CDNA_FILE" ]]; then
  echo "Downloading: $CDNA_URL"
  wget -c "$CDNA_URL" -O "$CDNA_FILE"
else
  echo "[INFO] cDNA already exists: $CDNA_FILE"
fi

if [[ ! -f "${CDNA_FILE%.gz}" ]]; then
  echo "Extracting cDNA..."
  gunzip -c "$CDNA_FILE" > "${CDNA_FILE%.gz}"
fi

echo ""
echo "------------------------------------------------------------------------"
echo "Step 3: Downloading ncRNA"
echo "------------------------------------------------------------------------"
NCRNA_URL="${BASE_URL}/fasta/${SPECIES}/ncrna/Mus_musculus.${ASSEMBLY}.ncrna.fa.gz"
NCRNA_FILE="Mus_musculus.${ASSEMBLY}.ncrna.fa.gz"

if [[ ! -f "$NCRNA_FILE" ]]; then
  echo "Downloading: $NCRNA_URL"
  wget -c "$NCRNA_URL" -O "$NCRNA_FILE"
else
  echo "[INFO] ncRNA already exists: $NCRNA_FILE"
fi

if [[ ! -f "${NCRNA_FILE%.gz}" ]]; then
  echo "Extracting ncRNA..."
  gunzip -c "$NCRNA_FILE" > "${NCRNA_FILE%.gz}"
fi

# Merge cDNA and ncRNA into transcriptome
TRANSCRIPTOME="Mus_musculus.${ASSEMBLY}.transcriptome.fa"
echo "Merging cDNA and ncRNA into transcriptome..."
if [[ ! -f "$TRANSCRIPTOME" ]]; then
  cat "${CDNA_FILE%.gz}" "${NCRNA_FILE%.gz}" > "$TRANSCRIPTOME"
fi

echo ""
echo "------------------------------------------------------------------------"
echo "Step 4: Downloading proteome"
echo "------------------------------------------------------------------------"
PEP_URL="${BASE_URL}/fasta/${SPECIES}/pep/Mus_musculus.${ASSEMBLY}.pep.all.fa.gz"
PEP_FILE="Mus_musculus.${ASSEMBLY}.pep.all.fa.gz"

if [[ ! -f "$PEP_FILE" ]]; then
  echo "Downloading: $PEP_URL"
  wget -c "$PEP_URL" -O "$PEP_FILE"
else
  echo "[INFO] Proteome already exists: $PEP_FILE"
fi

if [[ ! -f "${PEP_FILE%.gz}" ]]; then
  echo "Extracting proteome..."
  gunzip -c "$PEP_FILE" > "${PEP_FILE%.gz}"
fi

echo ""
echo "------------------------------------------------------------------------"
echo "Step 5: Generating P-sites BED file"
echo "------------------------------------------------------------------------"
PSITES_BED="psites.bed"

echo "Extracting start codons from GTF to generate P-sites..."
python3 - <<'PYTHON_SCRIPT' "$GTF_SORTED" "$PSITES_BED"
import sys

gtf_file = sys.argv[1]
output_bed = sys.argv[2]

with open(gtf_file, 'r') as gtf, open(output_bed, 'w') as bed:
    for line in gtf:
        if line.startswith('#'):
            continue

        fields = line.strip().split('\t')
        if len(fields) < 9:
            continue

        chrom = fields[0]
        feature = fields[2]
        start = int(fields[3])
        end = int(fields[4])
        strand = fields[6]

        # Only process start_codon features
        if feature != 'start_codon':
            continue

        # Extract transcript_id from attributes
        attrs = fields[8]
        transcript_id = ''
        for attr in attrs.split(';'):
            attr = attr.strip()
            if attr.startswith('transcript_id'):
                transcript_id = attr.split('"')[1]
                break

        if not transcript_id:
            continue

        # Write BED entry (convert to 0-based for BED)
        bed_start = start - 1
        bed_end = end
        bed.write(f"{chrom}\t{bed_start}\t{bed_end}\t{transcript_id}\t.\t{strand}\n")

print(f"Generated {output_bed}")
PYTHON_SCRIPT

echo ""
echo "------------------------------------------------------------------------"
echo "Step 6: Downloading transcript support levels (TSL/APPRIS)"
echo "------------------------------------------------------------------------"
echo "[INFO] Generating basic transcript support file..."
echo "[NOTE] For full TSL/APPRIS data, download from Ensembl BioMart manually"

# Generate a basic transcript support file from GTF
TSL_FILE="transcript_support_level.txt"
python3 - <<'PYTHON_SCRIPT' "$GTF_SORTED" "$TSL_FILE"
import sys
import re

gtf_file = sys.argv[1]
output_file = sys.argv[2]

transcripts = {}

with open(gtf_file, 'r') as gtf:
    for line in gtf:
        if line.startswith('#'):
            continue

        fields = line.strip().split('\t')
        if len(fields) < 9:
            continue

        feature = fields[2]
        if feature != 'transcript':
            continue

        attrs = fields[8]

        # Extract transcript_id
        transcript_id = ''
        tsl = 'NA'
        appris = 'NA'

        for attr in attrs.split(';'):
            attr = attr.strip()
            if attr.startswith('transcript_id'):
                transcript_id = attr.split('"')[1]
            elif 'transcript_support_level' in attr or 'tsl' in attr.lower():
                # Try to extract TSL value
                match = re.search(r'(\d+)', attr)
                if match:
                    tsl = match.group(1)

        if transcript_id:
            transcripts[transcript_id] = {'tsl': tsl, 'appris': appris}

# Write output
with open(output_file, 'w') as out:
    out.write("transcript_id\tTSL\tAPPRIS\n")
    for tid, data in transcripts.items():
        out.write(f"{tid}\t{data['tsl']}\t{data['appris']}\n")

print(f"Generated {output_file} with {len(transcripts)} transcripts")
PYTHON_SCRIPT

echo ""
echo "------------------------------------------------------------------------"
echo "Step 7: Creating standardized symlinks"
echo "------------------------------------------------------------------------"
# Create symlinks with standard names expected by gencode-riboseqORFs

ln -sf "$GTF_SORTED" "SORTED_TRANSCRIPTOME_GTF"
ln -sf "$TRANSCRIPTOME" "TRANSCRIPTOME_FASTA"
ln -sf "${PEP_FILE%.gz}" "PROTEOME_FASTA"
ln -sf "$TSL_FILE" "TRANSCRIPT_SUPPORT"
ln -sf "$PSITES_BED" "PSITES_BED"

echo "[OK] Created standardized symlinks:"
ls -lh SORTED_TRANSCRIPTOME_GTF TRANSCRIPTOME_FASTA PROTEOME_FASTA TRANSCRIPT_SUPPORT PSITES_BED

echo ""
echo "------------------------------------------------------------------------"
echo "Step 8: Cleanup"
echo "------------------------------------------------------------------------"
if [[ "$KEEP_DOWNLOADS" == "false" ]]; then
  echo "Removing compressed files..."
  rm -f *.gz
  echo "[OK] Cleanup complete"
else
  echo "[INFO] Keeping compressed files as requested"
fi

echo ""
echo "========================================================================"
echo "  ✓ Ensembl annotation preparation complete!"
echo "========================================================================"
echo ""
echo "Output directory: $OUTDIR"
echo ""
echo "Contents:"
ls -lh "$OUTDIR"
echo ""
echo "------------------------------------------------------------------------"
echo "Usage with gencode-riboseqORFs:"
echo "------------------------------------------------------------------------"
echo ""
echo "bash 16_gencode_orf_mapper.sh \\"
echo "  --project YOUR_PROJECT \\"
echo "  --fasta merged_orfs.fa \\"
echo "  --bed merged_orfs.bed \\"
echo "  --ensembl-dir $OUTDIR \\"
echo "  --image containers/biopython_1.81.sif"
echo ""
echo "========================================================================"
