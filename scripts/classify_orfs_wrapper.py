#!/usr/bin/env python3
import argparse
import sys
import os
import subprocess
from pathlib import Path

def run_command(cmd, desc):
    print(f"Running {desc}...", file=sys.stderr)
    print(" ".join(cmd), file=sys.stderr)
    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running {desc}: {e}", file=sys.stderr)
        sys.exit(1)

def get_script_path(script_name):
    # Assume scripts are in the same directory as this wrapper or subdirectories
    wrapper_dir = Path(__file__).parent.absolute()
    
    # Check common locations
    candidates = [
        wrapper_dir / script_name,
        wrapper_dir / "class_orf" / script_name,
        wrapper_dir / "gencode-riboseqORFs" / script_name,
        wrapper_dir.parent / "scripts" / "gencode-riboseqORFs" / script_name # In case we are in scripts/
    ]
    
    for c in candidates:
        if c.exists():
            return str(c)
            
    print(f"Error: Could not find script {script_name}", file=sys.stderr)
    sys.exit(1)

def main():
    parser = argparse.ArgumentParser(description="Unified ORF Classifier Wrapper")
    parser.add_argument("--mode", required=True, choices=['gencode', 'orfquant', 'orf_type'], help="Classification mode")
    parser.add_argument("--input", required=True, help="Input file prefix (from unify_orf_predictions.py output) or full path")
    parser.add_argument("--output_dir", required=True, help="Output directory")
    parser.add_argument("--gtf", help="Reference GTF file (Required for orfquant and orf_type)")
    parser.add_argument("--fasta", help="Reference FASTA file")
    parser.add_argument("--ensembl_dir", help="Ensembl directory (Required for gencode mode)")
    parser.add_argument("--cpus", type=str, default="1", help="Number of CPUs")
    
    args = parser.parse_args()
    
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Determine input files based on prefix or full path
    input_base = args.input
    if input_base.endswith('.gtf') or input_base.endswith('.bed') or input_base.endswith('.metadata.tsv'):
        # User provided full path, strip extension for base if needed, but keep path
        pass
    
    # Mode handling
    if args.mode == 'gencode':
        if not args.ensembl_dir:
            print("Error: --ensembl_dir is required for gencode mode", file=sys.stderr)
            sys.exit(1)
            
        script = get_script_path("ORF_mapper_to_GENCODE_v1.1.py")
        
        # Gencode mapper expects FASTA and BED.
        # Assume input is the prefix from unify_orf_predictions
        # If input is a file, try to deduce others or complain.
        
        # Let's assume input is the prefix "X" and we look for X.bed and X.metadata.tsv (to extract fasta?)
        # Actually gencode mapper needs FASTA with sequences.
        # unify_orf_predictions produces .bed and .metadata.tsv (which has sequences).
        # We might need to generate a fasta from metadata if not present.
        # Or if unify_orf_predictions produces sequences in metadata, we can generate it.
        # But wait, unify_orf_predictions doesn't produce .fa output by default in my new script?
        # Ah, I missed adding .fa output to unify_orf_predictions.
        # However, gencode mapper needs -f (fasta) and -b (bed).
        # The unified BED12 is suitable for -b.
        # The FASTA needs to be generated.
        
        # Helper: Generate FASTA from metadata if input is prefix
        # If input is a prefix like "out/unified", we look for "out/unified.metadata.tsv"
        metadata_file = f"{input_base}.metadata.tsv"
        fasta_file = f"{input_base}.orfs.fa"
        bed_file = f"{input_base}.bed"
        
        if not os.path.exists(fasta_file) and os.path.exists(metadata_file):
            print(f"Generating FASTA from {metadata_file}...", file=sys.stderr)
            with open(metadata_file, 'r') as m, open(fasta_file, 'w') as f:
                header = m.readline().strip().split('\t')
                try:
                    seq_idx = header.index('sequence')
                    id_idx = header.index('orf_id')
                    for line in m:
                        parts = line.strip().split('\t')
                        f.write(f">{parts[id_idx]}\n{parts[seq_idx]}\n")
                except ValueError:
                    print("Error: Metadata missing sequence or orf_id column", file=sys.stderr)
                    sys.exit(1)
        
        cmd = [
            "python3", script,
            "-d", args.ensembl_dir,
            "-f", fasta_file,
            "-b", bed_file,
            "-o", os.path.join(args.output_dir, "gencode_results")
        ]
        run_command(cmd, "GENCODE Classifier")

    elif args.mode == 'orfquant':
        if not args.gtf:
            print("Error: --gtf is required for orfquant mode", file=sys.stderr)
            sys.exit(1)
            
        script = get_script_path("run_orfquant_classify.R")
        
        # Input should be the unified GTF
        input_gtf = f"{input_base}.gtf" if not input_base.endswith('.gtf') else input_base
        
        output_file = os.path.join(args.output_dir, "orfquant_classification.tsv")
        
        cmd = [
            "Rscript", script,
            "--input", input_gtf,
            "--annotation", args.gtf,
            "--output", output_file
        ]
        run_command(cmd, "ORFquant Classifier")

    elif args.mode == 'orf_type':
        if not args.gtf:
            print("Error: --gtf is required for orf_type mode", file=sys.stderr)
            sys.exit(1)
            
        script = get_script_path("class_ORFtype.py")
        
        # Input should be metadata TSV
        input_tsv = f"{input_base}.metadata.tsv" if not input_base.endswith('.tsv') else input_base
        output_file = os.path.join(args.output_dir, "orftype_classification.tsv")
        
        cmd = [
            "python3", script,
            "--input", input_tsv,
            "--gtf", args.gtf,
            "--output", output_file,
            "--cpus", args.cpus
        ]
        run_command(cmd, "ORFtype Classifier")

if __name__ == "__main__":
    main()
