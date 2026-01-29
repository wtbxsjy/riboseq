#!/usr/bin/env bash
set -euo pipefail

# Quick diagnostic check for gencode converter output

usage() {
  cat <<'EOF'
Quick check for gencode-riboseqORFs format files

Usage:
  quick_check.sh --bed FILE.bed --fasta FILE.fa

This script performs basic validation:
1. Checks BED format (6 columns)
2. Verifies ORF naming consistency
3. Validates FASTA sequence lengths
4. Reports summary statistics

For detailed validation, use validate_output.py
EOF
}

BED=""
FASTA=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --bed) BED="$2"; shift 2;;
    --fasta) FASTA="$2"; shift 2;;
    -h|--help) usage; exit 0;;
    *) echo "Unknown arg: $1"; usage; exit 2;;
  esac
done

if [[ -z "$BED" || -z "$FASTA" ]]; then
  usage
  exit 2
fi

echo "=================================================="
echo "  Quick Check: GENCODE Converter Output"
echo "=================================================="
echo ""

# Check BED file
echo "📄 BED File: $BED"
echo "---"

if [[ ! -f "$BED" ]]; then
  echo "❌ File not found!"
  exit 1
fi

# Count lines
BED_LINES=$(wc -l < "$BED")
echo "   Total lines: $BED_LINES"

# Check column count
echo -n "   Checking format... "
BAD_COLS=$(awk -F'\t' 'NF != 6 {print NR}' "$BED" | head -5)
if [[ -n "$BAD_COLS" ]]; then
  echo "❌ Found lines with wrong column count:"
  for line in $BAD_COLS; do
    echo "      Line $line: $(sed -n "${line}p" "$BED" | awk -F'\t' '{print NF " columns"}')"
  done
else
  echo "✅ All lines have 6 columns"
fi

# Check ORF name format in column 4
echo -n "   Checking ORF names... "
BAD_NAMES=$(awk -F'\t' '$4 !~ /^[A-Z0-9]+_[0-9]+_[0-9]+aa$/ {print NR": "$4}' "$BED" | head -5)
if [[ -n "$BAD_NAMES" ]]; then
  echo "❌ Found invalid ORF names:"
  echo "$BAD_NAMES" | while read line; do
    echo "      $line"
  done
else
  echo "✅ All ORF names follow GENE_START_LENaa format"
fi

# Check for duplicates
echo -n "   Checking for duplicates... "
DUPLICATES=$(awk -F'\t' '{print $4}' "$BED" | sort | uniq -d | head -5)
if [[ -n "$DUPLICATES" ]]; then
  echo "⚠️  Found duplicate ORF names:"
  echo "$DUPLICATES" | while read orf; do
    count=$(grep -c "\\s${orf}\\s" "$BED")
    echo "      $orf (appears $count times)"
  done
else
  echo "✅ No duplicate ORF names"
fi

# Extract unique ORF names for comparison
awk -F'\t' '{print $4}' "$BED" | sort > /tmp/bed_orfs.txt

echo ""

# Check FASTA file
echo "📄 FASTA File: $FASTA"
echo "---"

if [[ ! -f "$FASTA" ]]; then
  echo "❌ File not found!"
  exit 1
fi

# Count entries
FASTA_ORFS=$(grep -c "^>" "$FASTA")
echo "   Total ORFs: $FASTA_ORFS"

# Check header format
echo -n "   Checking headers... "
BAD_HEADERS=$(grep "^>" "$FASTA" | grep -vP '^>[A-Z0-9]+_[0-9]+_[0-9]+aa--\S+$' | head -5)
if [[ -n "$BAD_HEADERS" ]]; then
  echo "❌ Found invalid headers:"
  echo "$BAD_HEADERS" | while read line; do
    echo "      $line"
  done
else
  echo "✅ All headers follow GENE_START_LENaa--STUDY format"
fi

# Check stop codons
echo -n "   Checking stop codons... "
MISSING_STOP=$(awk '/^>/ {if (seq && seq !~ /\*$/) print header; header=$0; seq=""} !/^>/ {seq=seq$0} END {if (seq && seq !~ /\*$/) print header}' "$FASTA" | head -5)
if [[ -n "$MISSING_STOP" ]]; then
  echo "⚠️  Found sequences without stop codon:"
  echo "$MISSING_STOP" | while read line; do
    echo "      $line"
  done
else
  echo "✅ All sequences end with stop codon (*)"
fi

# Check length consistency
echo -n "   Checking sequence lengths... "
python3 - <<'PYCHECK' "$FASTA"
import sys
import re

issues = []
with open(sys.argv[1], 'r') as f:
    header = None
    seq = []
    for line in f:
        line = line.strip()
        if line.startswith('>'):
            if header and seq:
                seq_str = ''.join(seq)
                match = re.search(r'_(\d+)aa--', header)
                if match:
                    declared = int(match.group(1))
                    actual = len(seq_str.rstrip('*'))
                    if declared != actual:
                        issues.append(f"{header}: declared {declared}aa, actual {actual}aa")
                        if len(issues) >= 5:
                            break
            header = line[1:]
            seq = []
        else:
            seq.append(line)
    
    # Check last entry
    if header and seq:
        seq_str = ''.join(seq)
        match = re.search(r'_(\d+)aa--', header)
        if match:
            declared = int(match.group(1))
            actual = len(seq_str.rstrip('*'))
            if declared != actual and len(issues) < 5:
                issues.append(f"{header}: declared {declared}aa, actual {actual}aa")

if issues:
    print("❌ Found length mismatches:")
    for issue in issues:
        print(f"      {issue}")
    sys.exit(1)
else:
    print("✅ All sequences match their declared lengths")
PYCHECK

if [[ $? -ne 0 ]]; then
  echo "   ⚠️  Run with validate_output.py for details"
fi

# Extract ORF names from FASTA
grep "^>" "$FASTA" | sed 's/^>//;s/--.*$//' | sort > /tmp/fasta_orfs.txt

echo ""

# Cross-check
echo "🔗 Cross-Validation"
echo "---"

ONLY_BED=$(comm -23 /tmp/bed_orfs.txt /tmp/fasta_orfs.txt | wc -l)
ONLY_FASTA=$(comm -13 /tmp/bed_orfs.txt /tmp/fasta_orfs.txt | wc -l)

if [[ $ONLY_BED -gt 0 ]]; then
  echo "   ⚠️  $ONLY_BED ORFs only in BED file"
  comm -23 /tmp/bed_orfs.txt /tmp/fasta_orfs.txt | head -3 | while read orf; do
    echo "      - $orf"
  done
  if [[ $ONLY_BED -gt 3 ]]; then
    echo "      ... and $((ONLY_BED - 3)) more"
  fi
fi

if [[ $ONLY_FASTA -gt 0 ]]; then
  echo "   ⚠️  $ONLY_FASTA ORFs only in FASTA file"
  comm -13 /tmp/bed_orfs.txt /tmp/fasta_orfs.txt | head -3 | while read orf; do
    echo "      - $orf"
  done
  if [[ $ONLY_FASTA -gt 3 ]]; then
    echo "      ... and $((ONLY_FASTA - 3)) more"
  fi
fi

if [[ $ONLY_BED -eq 0 && $ONLY_FASTA -eq 0 ]]; then
  echo "   ✅ BED and FASTA have matching ORF sets"
fi

# Cleanup
rm -f /tmp/bed_orfs.txt /tmp/fasta_orfs.txt

echo ""
echo "=================================================="
echo "Summary:"
echo "  BED entries: $BED_LINES"
echo "  FASTA ORFs: $FASTA_ORFS"
echo ""

if [[ $ONLY_BED -eq 0 && $ONLY_FASTA -eq 0 ]]; then
  echo "✅ Quick check PASSED"
  echo ""
  echo "For detailed validation, run:"
  echo "  python3 validate_output.py --bed $BED --fasta $FASTA"
  exit 0
else
  echo "⚠️  Issues detected"
  echo ""
  echo "For detailed analysis, run:"
  echo "  python3 validate_output.py --bed $BED --fasta $FASTA"
  exit 1
fi
