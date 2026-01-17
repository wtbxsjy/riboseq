#!/bin/bash
# Quick test script to validate all GENCODE format converters
# This script runs all converters with test data and validates outputs

set -e  # Exit on error

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

echo "╔════════════════════════════════════════════════════════════╗"
echo "║  GENCODE Format Converters - Integration Test Suite       ║"
echo "╚════════════════════════════════════════════════════════════╝"
echo ""

# Colors for output
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

success_count=0
fail_count=0

# Function to print test status
print_status() {
    if [ $1 -eq 0 ]; then
        echo -e "${GREEN}✓ PASS${NC}: $2"
        ((success_count++))
    else
        echo -e "${RED}✗ FAIL${NC}: $2"
        ((fail_count++))
    fi
}

# Function to validate FASTA format
validate_fasta() {
    local file=$1
    local min_orfs=$2

    if [ ! -f "$file" ]; then
        echo "File not found: $file"
        return 1
    fi

    # Check for FASTA headers
    local header_count=$(grep -c "^>" "$file" || true)
    if [ "$header_count" -lt "$min_orfs" ]; then
        echo "Expected at least $min_orfs ORFs, found $header_count"
        return 1
    fi

    # Check header format: >{NAME}--{STUDY_ID}
    if ! grep -q "^>[^-]*--[^-]*$" "$file"; then
        echo "Invalid FASTA header format (expected: >NAME--STUDY_ID)"
        return 1
    fi

    # Check for stop codons
    if ! grep -v "^>" "$file" | grep -q "\*$"; then
        echo "Sequences must end with stop codon (*)"
        return 1
    fi

    return 0
}

# Function to validate BED format
validate_bed() {
    local file=$1
    local min_orfs=$2

    if [ ! -f "$file" ]; then
        echo "File not found: $file"
        return 1
    fi

    # Check line count
    local line_count=$(wc -l < "$file")
    if [ "$line_count" -lt "$min_orfs" ]; then
        echo "Expected at least $min_orfs ORFs, found $line_count"
        return 1
    fi

    # Check 6 columns (chr, start, end, name, study_id, strand)
    local col_count=$(awk '{print NF}' "$file" | sort -u)
    if [ "$col_count" != "6" ]; then
        echo "Expected 6 columns, found $col_count"
        return 1
    fi

    # Check strand column (+ or -)
    if ! awk '{print $6}' "$file" | grep -q "^[+-]$"; then
        echo "Invalid strand values (expected: + or -)"
        return 1
    fi

    return 0
}

echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "Test 1: Ribo-TISH to GENCODE Converter"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

cd "$PROJECT_ROOT/test_data/ribotish_to_gencode"

echo "Running ribotish_to_gencode.py..."
python3 ../../bin/ribotish_to_gencode.py \
    --predict ribotish_predict.txt \
    --fasta test_genome.fa \
    --study_id TEST_RIBOTISH \
    --output_prefix test_ribotish \
    --min_length 16 > /dev/null 2>&1

validate_fasta "test_ribotish.gencode.fa" 8
print_status $? "Ribo-TISH FASTA output validation"

validate_bed "test_ribotish.gencode.bed" 8
print_status $? "Ribo-TISH BED output validation"

# Check FASTA-BED consistency
fasta_count=$(grep -c "^>" test_ribotish.gencode.fa)
bed_count=$(wc -l < test_ribotish.gencode.bed)
if [ "$fasta_count" -eq "$bed_count" ]; then
    print_status 0 "Ribo-TISH FASTA-BED consistency ($fasta_count ORFs)"
else
    print_status 1 "Ribo-TISH FASTA-BED consistency (FASTA: $fasta_count, BED: $bed_count)"
fi

echo ""

echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "Test 2: Ribotricer to GENCODE Converter"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

cd "$PROJECT_ROOT/test_data/ribotricer_to_gencode"

echo "Running ribotricer_to_gencode.py..."
python3 ../../bin/ribotricer_to_gencode.py \
    --tsv ribotricer_translating_ORFs.tsv \
    --fasta test_genome.fa \
    --study_id TEST_RIBOTRICER \
    --output_prefix test_ribotricer \
    --min_length 16 \
    --min_phase_score 0.5 > /dev/null 2>&1

validate_fasta "test_ribotricer.gencode.fa" 8
print_status $? "Ribotricer FASTA output validation"

validate_bed "test_ribotricer.gencode.bed" 8
print_status $? "Ribotricer BED output validation"

# Check FASTA-BED consistency
fasta_count=$(grep -c "^>" test_ribotricer.gencode.fa)
bed_count=$(wc -l < test_ribotricer.gencode.bed)
if [ "$fasta_count" -eq "$bed_count" ]; then
    print_status 0 "Ribotricer FASTA-BED consistency ($fasta_count ORFs)"
else
    print_status 1 "Ribotricer FASTA-BED consistency (FASTA: $fasta_count, BED: $bed_count)"
fi

echo ""

echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "Test 3: Multi-tool Integration Test"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

TEMP_DIR="$PROJECT_ROOT/test_data/integration_test"
mkdir -p "$TEMP_DIR"

echo "Merging outputs from both converters..."

# Merge FASTA files
cat "$PROJECT_ROOT/test_data/ribotish_to_gencode/test_ribotish.gencode.fa" \
    "$PROJECT_ROOT/test_data/ribotricer_to_gencode/test_ribotricer.gencode.fa" \
    > "$TEMP_DIR/merged.gencode.fa"

# Merge BED files
cat "$PROJECT_ROOT/test_data/ribotish_to_gencode/test_ribotish.gencode.bed" \
    "$PROJECT_ROOT/test_data/ribotricer_to_gencode/test_ribotricer.gencode.bed" \
    > "$TEMP_DIR/merged.gencode.bed"

# Validate merged files
validate_fasta "$TEMP_DIR/merged.gencode.fa" 16
print_status $? "Merged FASTA validation (multi-tool)"

validate_bed "$TEMP_DIR/merged.gencode.bed" 16
print_status $? "Merged BED validation (multi-tool)"

# Check for study ID separation
if grep -q "TEST_RIBOTISH" "$TEMP_DIR/merged.gencode.fa" && \
   grep -q "TEST_RIBOTRICER" "$TEMP_DIR/merged.gencode.fa"; then
    print_status 0 "Study ID preservation in merged file"
else
    print_status 1 "Study ID preservation in merged file"
fi

# Check coordinate ranges
echo "Checking coordinate ranges..."
awk '{print $2, $3}' "$TEMP_DIR/merged.gencode.bed" | \
    awk '{if ($1 >= $2) print "Invalid coordinates: " $1 " >= " $2}' > "$TEMP_DIR/coord_check.txt"

if [ ! -s "$TEMP_DIR/coord_check.txt" ]; then
    print_status 0 "BED coordinate validation (start < end)"
else
    cat "$TEMP_DIR/coord_check.txt"
    print_status 1 "BED coordinate validation (start < end)"
fi

echo ""

echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "Test 4: Format Compliance Tests"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

# Test FASTA header format compliance
echo "Validating FASTA header format..."
if grep "^>" "$TEMP_DIR/merged.gencode.fa" | grep -qE "^>[A-Za-z0-9_]+-*[A-Za-z0-9_]*--[A-Za-z0-9_]+$"; then
    print_status 0 "FASTA header format compliance"
else
    print_status 1 "FASTA header format compliance"
fi

# Test BED strand values
echo "Validating BED strand values..."
if awk '{print $6}' "$TEMP_DIR/merged.gencode.bed" | grep -vqE "^[+-]$"; then
    print_status 1 "BED strand value compliance"
else
    print_status 0 "BED strand value compliance"
fi

# Test chromosome naming
echo "Validating chromosome naming..."
if awk '{print $1}' "$TEMP_DIR/merged.gencode.bed" | grep -qE "^chr[0-9XYM]+$"; then
    print_status 0 "Chromosome naming convention"
else
    print_status 1 "Chromosome naming convention"
fi

echo ""

echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "Summary Statistics"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

echo ""
echo "Ribo-TISH Output:"
echo "  - ORFs detected: $(grep -c "^>" "$PROJECT_ROOT/test_data/ribotish_to_gencode/test_ribotish.gencode.fa")"
echo "  - Chromosomes: $(awk '{print $1}' "$PROJECT_ROOT/test_data/ribotish_to_gencode/test_ribotish.gencode.bed" | sort -u | wc -l)"
echo "  - Strand distribution: +$(awk '$6=="+"' "$PROJECT_ROOT/test_data/ribotish_to_gencode/test_ribotish.gencode.bed" | wc -l) / -$(awk '$6=="-"' "$PROJECT_ROOT/test_data/ribotish_to_gencode/test_ribotish.gencode.bed" | wc -l)"

echo ""
echo "Ribotricer Output:"
echo "  - ORFs detected: $(grep -c "^>" "$PROJECT_ROOT/test_data/ribotricer_to_gencode/test_ribotricer.gencode.fa")"
echo "  - Chromosomes: $(awk '{print $1}' "$PROJECT_ROOT/test_data/ribotricer_to_gencode/test_ribotricer.gencode.bed" | sort -u | wc -l)"
echo "  - Strand distribution: +$(awk '$6=="+"' "$PROJECT_ROOT/test_data/ribotricer_to_gencode/test_ribotricer.gencode.bed" | wc -l) / -$(awk '$6=="-"' "$PROJECT_ROOT/test_data/ribotricer_to_gencode/test_ribotricer.gencode.bed" | wc -l)"

echo ""
echo "Merged Output:"
echo "  - Total ORFs: $(grep -c "^>" "$TEMP_DIR/merged.gencode.fa")"
echo "  - Total chromosomes: $(awk '{print $1}' "$TEMP_DIR/merged.gencode.bed" | sort -u | wc -l)"
echo "  - Studies represented: $(grep "^>" "$TEMP_DIR/merged.gencode.fa" | sed 's/.*--//' | sort -u | wc -l)"

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "Test Results"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

total_tests=$((success_count + fail_count))
success_rate=$(awk "BEGIN {printf \"%.1f\", ($success_count / $total_tests) * 100}")

echo -e "Tests Passed: ${GREEN}$success_count${NC} / $total_tests"
echo -e "Tests Failed: ${RED}$fail_count${NC} / $total_tests"
echo -e "Success Rate: ${YELLOW}$success_rate%${NC}"
echo ""

if [ $fail_count -eq 0 ]; then
    echo -e "${GREEN}╔════════════════════════════════════════════╗${NC}"
    echo -e "${GREEN}║  ✓ ALL TESTS PASSED                       ║${NC}"
    echo -e "${GREEN}║  Converters are ready for integration!    ║${NC}"
    echo -e "${GREEN}╚════════════════════════════════════════════╝${NC}"
    exit 0
else
    echo -e "${RED}╔════════════════════════════════════════════╗${NC}"
    echo -e "${RED}║  ✗ SOME TESTS FAILED                       ║${NC}"
    echo -e "${RED}║  Please review the output above            ║${NC}"
    echo -e "${RED}╚════════════════════════════════════════════╝${NC}"
    exit 1
fi
