#!/usr/bin/env bash
set -euo pipefail

################################################################################
# Script: patch_functions_for_single_sample.sh
# Purpose: Patch functions.py to work in single-sample mode
################################################################################

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FUNCTIONS_FILE="$SCRIPT_DIR/../gencode-riboseqORFs/functions.py"

if [[ ! -f "$FUNCTIONS_FILE" ]]; then
  echo "[ERROR] functions.py not found at: $FUNCTIONS_FILE"
  exit 1
fi

echo "========================================================================"
echo "  Patching functions.py for single-sample mode"
echo "========================================================================"
echo ""
echo "File: $FUNCTIONS_FILE"
echo ""

# Create backup
BACKUP_FILE="${FUNCTIONS_FILE}.backup_$(date +%Y%m%d_%H%M%S)"
cp "$FUNCTIONS_FILE" "$BACKUP_FILE"
echo "[INFO] Created backup: $BACKUP_FILE"
echo ""

# Apply patch using sed
echo "[INFO] Applying patch..."

# Find line 508 and replace the problematic code block
sed -i.tmp '508s/.*/\t# Modified to handle missing file (single-sample mode)\n\triboseq_list_file = "list_riboseqs\/list_riboseq_orfs.txt"\n\tif os.path.exists(riboseq_list_file):\n\t\tfor line in open(riboseq_list_file):/' "$FUNCTIONS_FILE"

# Add proper indentation to the subsequent lines (509-540)
# All lines inside the for loop need extra indentation
sed -i.tmp '509,540s/^\t/\t\t/' "$FUNCTIONS_FILE"

# Add the else clause after the for loop block
sed -i.tmp '541s/.*/\telse:\n\t\tprint("  [INFO] Riboseq ORF list file not found, running in single-sample mode")\n\n&/' "$FUNCTIONS_FILE"

# Clean up temp files
rm -f "${FUNCTIONS_FILE}.tmp"

echo "[OK] Patch applied successfully!"
echo ""
echo "Changes made:"
echo "  - Added file existence check for list_riboseqs/list_riboseq_orfs.txt"
echo "  - Script will now work in single-sample mode without this file"
echo ""
echo "If you want to revert, restore from: $BACKUP_FILE"
echo ""
