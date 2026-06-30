#!/bin/bash

# ==============================================================================
# Script: sra_to_fastq_gz_v2.sh
# Description: Converts SRA files to compressed FASTQ (.fastq.gz) files
#              with flexible output directory control.
# Author: Gemini
# Date: 2025-11-09
#
# Usage:
#   ./sra_to_fastq_gz_v2.sh [-t <dump_threads>] [-p <pigz_threads>] [-o <out_dir>] <SRA_FILE_1> ...
#
# Dependencies: sra-tools (for fasterq-dump) and pigz
# ==============================================================================

# --- Default Settings ---
THREADS_FASTERQ=8
THREADS_PIGZ=8
# Global output directory (empty by default, meaning use SRA file's dir)
GLOBAL_OUTPUT_DIR=""

# --- Functions ---

# Function to display usage information
show_usage() {
    echo "Usage: $0 [-t <dump_threads>] [-p <pigz_threads>] [-o <out_dir>] <SRA_FILE_1> [SRA_FILE_2] ..."
    echo
    echo "Options:"
    echo "  -t <threads>   Number of threads for fasterq-dump (default: $THREADS_FASTERQ)"
    echo "  -p <threads>   Number of threads for pigz (default: $THREADS_PIGZ)"
    echo "  -o <out_dir>   Specify a single output directory for all files."
    echo "                 (Default: Output each fastq.gz to its corresponding SRA file's directory)"
    echo "  -h             Show this help message"
    echo
    echo "Example:"
    echo "  # Process file, output to a specific dir '/data/fastq_out'"
    echo "  $0 -t 16 -o /data/fastq_out SRR1234567.sra"
    echo
    echo "  # Process multiple files, outputting each to its own source directory"
    echo "  $0 -p 4 /path/to/SRR123.sra /another/path/SRR456.sra"
}

# --- Parse Command-Line Options ---

# Note the added 'o:' to parse the -o option
while getopts ":t:p:o:h" opt; do
  case $opt in
    t)
      THREADS_FASTERQ="$OPTARG"
      ;;
    p)
      THREADS_PIGZ="$OPTARG"
      ;;
    o)
      GLOBAL_OUTPUT_DIR="$OPTARG"
      ;;
    h)
      show_usage
      exit 0
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      show_usage
      exit 1
      ;;
    :)
      echo "Option -$OPTARG requires an argument." >&2
      show_usage
      exit 1
      ;;
  esac
done

# Remove parsed options from the argument list
shift $((OPTIND-1))

# --- Check for Input Files ---
if [ $# -eq 0 ]; then
    echo "Error: No SRA files provided."
    show_usage
    exit 1
fi

# --- Check for Dependencies ---
# (Dependency checks remain the same...)
#if ! command -v fasterq-dump &> /dev/null; then
#    echo "Error: 'fasterq-dump' not found. Please install NCBI SRA Toolkit (sra-tools)."
#    exit 1
#fi

#if ! command -v pigz &> /dev/null; then
#    echo "Error: 'pigz' not found. Please install pigz."
#    exit 1
#fi

# --- Main Processing Loop ---

echo "--- Starting SRA to FASTQ.gz Conversion ---"
echo "fasterq-dump threads: $THREADS_FASTERQ"
echo "pigz threads: $THREADS_PIGZ"
if [ -n "$GLOBAL_OUTPUT_DIR" ]; then
    echo "Global output directory: $GLOBAL_OUTPUT_DIR"
    # Create the global output directory if it doesn't exist
    mkdir -p "$GLOBAL_OUTPUT_DIR"
    if [ $? -ne 0 ]; then
        echo "Error: Could not create output directory '$GLOBAL_OUTPUT_DIR'. Aborting."
        exit 1
    fi
else
    echo "Output directory: Defaulting to each SRA file's source directory."
fi
echo "-----------------------------------------------"

# Loop through all SRA files provided as arguments
for sra_file in "$@"; do
    
    if [ ! -f "$sra_file" ]; then
        echo "Warning: File '$sra_file' not found. Skipping."
        continue
    fi

    # --- Determine Final Output Directory ---
    final_output_dir=""
    if [ -n "$GLOBAL_OUTPUT_DIR" ]; then
        # Use the user-specified global directory
        final_output_dir="$GLOBAL_OUTPUT_DIR"
    else
        # Default: Use the directory containing the SRA file
        final_output_dir=$(dirname "$sra_file")
        
        # Handle case where dirname might be "." (current dir)
        if [ "$final_output_dir" = "." ]; then
            final_output_dir=$(pwd)
        fi
    fi

    # Get the base accession name (e.g., "SRR1234567" from "SRR1234567.sra")
    accession=$(basename "$sra_file" .sra)

    echo "[1/2] Processing '$sra_file'..."
    echo "      Outputting to: '$final_output_dir'"
    
    # Run fasterq-dump
    # -e : threads
    # -S : Split files for paired-end
    # -p : Show progress
    # -O : Specify output directory <-- NEW
    fasterq-dump -e "$THREADS_FASTERQ" -S -p -O "$final_output_dir" "$sra_file"

    # Check if fasterq-dump was successful
    if [ $? -ne 0 ]; then
        echo "Error: fasterq-dump failed for '$sra_file'."
        continue # Skip to the next file
    fi

    echo "[2/2] Compressing output for '$accession' in '$final_output_dir'..."
    
    # Find all .fastq files generated FOR THIS ACCESSION in the output dir
    # This is more robust, as it looks in the correct output directory
    fastq_files=$(find "$final_output_dir" -maxdepth 1 -type f -name "${accession}*.fastq")

    if [ -z "$fastq_files" ]; then
        echo "Warning: No .fastq files found for $accession in $final_output_dir. Dump may have failed."
        continue
    fi

    # Compress each generated file
    for fq in $fastq_files; do
        echo "  -> Compressing '$fq'..."
        pigz -p "$THREADS_PIGZ" "$fq"
        
        if [ $? -eq 0 ]; then
            echo "  -> Successfully created '$fq.gz'"
        else
            echo "  -> Error: pigz failed for '$fq'"
        fi
    done

    echo "--- Finished processing '$sra_file' ---"
    echo
done

echo "-----------------------------------------------"
echo "All processing complete."

