#!/bin/bash
# Quick Start Example for prepare_workflow.py
# This example demonstrates basic usage with test data

set -euo pipefail

echo "========================================"
echo "Prepare Workflow - Quick Start Example"
echo "========================================"
echo ""

# Configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORKDIR="$HOME/riboseq_test_project"
DATA_DIR="/path/to/your/fastq_files"  # CHANGE THIS
REFERENCE_DIR="/path/to/your/reference"  # CHANGE THIS (optional)
CONTAINER_DIR="/path/to/your/containers"  # CHANGE THIS (optional)

# Detect OS
if [[ "$OSTYPE" == "linux-gnu"* ]]; then
    PLATFORM="Linux"
elif [[ "$OSTYPE" == "darwin"* ]]; then
    PLATFORM="macOS"
else
    PLATFORM="Unknown"
fi

echo "Platform: $PLATFORM"
echo "Script directory: $SCRIPT_DIR"
echo ""

# Check if prepare_workflow.py exists
PREPARE_SCRIPT="$SCRIPT_DIR/prepare_workflow.py"
if [[ ! -f "$PREPARE_SCRIPT" ]]; then
    echo "ERROR: prepare_workflow.py not found at $PREPARE_SCRIPT"
    exit 1
fi

# Check if data directory exists
if [[ ! -d "$DATA_DIR" ]]; then
    echo "WARNING: Data directory not found: $DATA_DIR"
    echo "Please edit this script and set DATA_DIR to your FASTQ directory"
    echo ""
    echo "Example usage:"
    echo "  1. Edit this script: nano $0"
    echo "  2. Change DATA_DIR='/path/to/your/fastq_files'"
    echo "  3. Run again: bash $0"
    exit 1
fi

# Count FASTQ files
FASTQ_COUNT=$(find "$DATA_DIR" -maxdepth 1 -name "*.fastq.gz" -o -name "*.fq.gz" 2>/dev/null | wc -l)
echo "Found $FASTQ_COUNT FASTQ files in $DATA_DIR"
echo ""

if [[ $FASTQ_COUNT -eq 0 ]]; then
    echo "WARNING: No FASTQ files found in data directory"
    echo "Please check your DATA_DIR setting"
    exit 1
fi

# Ask for confirmation
echo "This will create a new workflow in: $WORKDIR"
echo ""
read -p "Continue? (y/N): " -n 1 -r
echo ""

if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    echo "Cancelled."
    exit 0
fi

# Build command
CMD="python3 $PREPARE_SCRIPT -w $WORKDIR -d $DATA_DIR --genome GRCh38 --species human"

# Add optional parameters if directories exist
if [[ -d "$REFERENCE_DIR" ]]; then
    CMD="$CMD -r $REFERENCE_DIR"
    echo "Adding reference directory: $REFERENCE_DIR"
fi

if [[ -d "$CONTAINER_DIR" ]]; then
    CMD="$CMD -c $CONTAINER_DIR"
    echo "Adding container directory: $CONTAINER_DIR"
fi

echo ""
echo "========================================"
echo "Running prepare_workflow.py"
echo "========================================"
echo ""
echo "Command: $CMD"
echo ""

# Execute
$CMD

# Check result
if [[ $? -eq 0 ]]; then
    echo ""
    echo "========================================"
    echo "SUCCESS!"
    echo "========================================"
    echo ""
    echo "Workflow prepared in: $WORKDIR"
    echo ""
    echo "Next steps:"
    echo "  1. Review sample sheet: $WORKDIR/scripts/samplesheet.csv"
    echo "  2. Review configuration: $WORKDIR/scripts/workflow_config.json"
    echo "  3. Run pipeline:"
    echo "     cd $WORKDIR"
    echo "     bash run_pipeline.sh"
    echo ""
    echo "Directory structure:"
    ls -la "$WORKDIR"
else
    echo ""
    echo "ERROR: prepare_workflow.py failed"
    echo "Please check the error messages above"
    exit 1
fi
