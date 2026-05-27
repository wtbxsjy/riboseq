#!/usr/bin/env python3
import argparse
import sys
import os
import gzip
import subprocess
from pathlib import Path

try:
    from orfont.core.utils import chrom_aliases, CODON_TABLE, translate_nt as _translate_nt
    def translate_nt(nt_seq): return _translate_nt(nt_seq)
except ImportError:
    # Fallback definitions when orfont package is not installed
    def chrom_aliases(chrom):
        if not chrom: return []
        aliases = [chrom]
        if chrom.startswith("chr"):
            base = chrom[3:]; aliases.append(base)
            if base == "M": aliases.append("MT")
            elif base == "MT": aliases.append("M")
        else:
            aliases.append(f"chr{chrom}")
            if chrom == "MT": aliases.extend(["chrM", "M"])
            elif chrom == "M": aliases.extend(["chrM", "MT"])
        out = []; seen = set()
        for a in aliases:
            if a and a not in seen: seen.add(a); out.append(a)
        return out

    CODON_TABLE = {
        'TTT':'F','TTC':'F','TTA':'L','TTG':'L','CTT':'L','CTC':'L','CTA':'L','CTG':'L',
        'ATT':'I','ATC':'I','ATA':'I','ATG':'M','GTT':'V','GTC':'V','GTA':'V','GTG':'V',
        'TCT':'S','TCC':'S','TCA':'S','TCG':'S','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
        'ACT':'T','ACC':'T','ACA':'T','ACG':'T','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
        'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
        'AAT':'N','AAC':'N','AAA':'K','AAG':'K','GAT':'D','GAC':'D','GAA':'E','GAG':'E',
        'TGT':'C','TGC':'C','TGA':'*','TGG':'W','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
        'AGT':'S','AGC':'S','AGA':'R','AGG':'R','GGT':'G','GGC':'G','GGA':'G','GGG':'G',
    }
    def translate_nt(nt_seq):
        nt_seq = nt_seq.upper().replace('U', 'T')
        aa = []
        for i in range(0, len(nt_seq) - 2, 3):
            codon = nt_seq[i:i+3]; aa.append(CODON_TABLE.get(codon, 'X'))
        return ''.join(aa)


def _open(path, mode='r'):
    if path.endswith('.gz'):
        return gzip.open(path, mode + 't')
    return open(path, mode)


def build_reference_chrom_map(gtf_path):
    ref_chroms = set()
    with _open(gtf_path) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t", 1)
            if parts and parts[0]:
                ref_chroms.add(parts[0])

    chrom_map = {chrom: chrom for chrom in ref_chroms}
    for chrom in ref_chroms:
        for alias in chrom_aliases(chrom):
            chrom_map.setdefault(alias, chrom)
    return chrom_map

def _resolve(prefix, suffix):
    """Return the file path, trying .gz extension first."""
    gz = f"{prefix}{suffix}.gz"
    if os.path.exists(gz):
        return gz
    return f"{prefix}{suffix}"

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

def main(argv=None):
    parser = argparse.ArgumentParser(description="Unified ORF Classifier Wrapper")
    parser.add_argument("--mode", required=True, choices=['gencode', 'orfquant', 'orf_type'], help="Classification mode")
    parser.add_argument("--input", required=True, help="Input file prefix (from unify_orf_predictions.py output) or full path")
    parser.add_argument("--output_dir", required=True, help="Output directory")
    parser.add_argument("--gtf", help="Reference GTF file (Required for orfquant and orf_type)")
    parser.add_argument("--fasta", help="Reference FASTA file")
    parser.add_argument("--ensembl_dir", help="Ensembl directory (Required for gencode mode)")
    parser.add_argument("--gencode_impl", choices=["original", "fast", "indexed_fast"], default="original", help="Implementation to use for gencode mode")
    parser.add_argument("--cpus", type=str, default="1", help="Number of CPUs")

    args = parser.parse_args(argv)

    os.makedirs(args.output_dir, exist_ok=True)
    
    # Determine input files based on prefix or full path
    input_base = args.input
    if input_base.endswith('.gtf') or input_base.endswith('.gtf.gz') or input_base.endswith('.bed') or input_base.endswith('.bed.gz') or input_base.endswith('.metadata.tsv'):
        # User provided full path, strip extension for base if needed, but keep path
        pass
    
    # Mode handling
    if args.mode == 'gencode':
        if not args.ensembl_dir:
            print("Error: --ensembl_dir is required for gencode mode", file=sys.stderr)
            sys.exit(1)

        if args.gencode_impl == "fast":
            script = get_script_path("run_gencode_classify_fast.py")
            cmd = [
                "python3", script,
                "--input_prefix", input_base,
                "--ensembl_dir", args.ensembl_dir,
                "--output_dir", args.output_dir,
                "--cpus", args.cpus
            ]
            run_command(cmd, "GENCODE Classifier (fast)")
            return

        if args.gencode_impl == "indexed_fast":
            script = get_script_path("run_gencode_classify_indexed.py")
            cmd = [
                "python3", script,
                "--input_prefix", input_base,
                "--ensembl_dir", args.ensembl_dir,
                "--output_dir", args.output_dir,
                "--cpus", args.cpus
            ]
            run_command(cmd, "GENCODE Classifier (indexed_fast)")
            return

        script = get_script_path("ORF_mapper_to_GENCODE_v1.1.py")

        metadata_file = f"{input_base}.metadata.tsv"
        fasta_file    = f"{input_base}.orfs.fa"   # protein (AA) FASTA for the mapper
        bed12_file    = _resolve(input_base, '.bed')        # our BED12 output (may be .bed.gz)
        bed6_file     = f"{input_base}.orfs.bed6"  # BED6 required by the mapper
        chrom_map = build_reference_chrom_map(os.path.join(args.ensembl_dir, "SORTED_TRANSCRIPTOME_GTF"))

        # ── Step 1: read metadata to build orf_id → study_id mapping ──────────
        # The GENCODE mapper uses BED col[4] as the "study" identifier and
        # constructs FASTA keys as "{orf_id}--{study_id}".  We derive study_id
        # from the first sample in the metadata `samples` column.
        orf_to_study = {}
        if os.path.exists(metadata_file):
            with _open(metadata_file) as mf:
                hdr = mf.readline().strip().split('\t')
                id_idx  = hdr.index('orf_id')
                try:
                    smp_idx = hdr.index('samples')
                except ValueError:
                    smp_idx = -1
                for line in mf:
                    parts = line.strip().split('\t')
                    if len(parts) <= id_idx:
                        continue
                    orf_id   = parts[id_idx]
                    study_id = (parts[smp_idx].split(',')[0]
                                if smp_idx >= 0 and smp_idx < len(parts) and parts[smp_idx]
                                else 'unified')
                    orf_to_study[orf_id] = study_id

        # ── Step 2: generate BED6 from BED12 ──────────────────────────────────
        # The mapper expects BED6: chr start end name study_id strand
        # (col[4] = study ID, not numeric score; col[5] = strand).
        # Our BED12 has col[5]=strand but col[4]="0" (score); we replace it.
        print(f"Generating BED6 for GENCODE mapper from {bed12_file}...", file=sys.stderr)
        remapped_chroms = 0
        with _open(bed12_file) as bf, open(bed6_file, 'w') as b6:
            for line in bf:
                parts = line.rstrip('\n').split('\t')
                if len(parts) < 6:
                    continue
                chrom, start, end, name, _, strand = parts[:6]
                mapped_chrom = chrom_map.get(chrom, chrom)
                if mapped_chrom != chrom:
                    remapped_chroms += 1
                study_id = orf_to_study.get(name, 'unified')
                b6.write(f"{mapped_chrom}\t{start}\t{end}\t{name}\t{study_id}\t{strand}\n")
        if remapped_chroms:
            print(f"Remapped {remapped_chroms} BED records to reference chromosome names", file=sys.stderr)

        # ── Step 3: generate protein FASTA ────────────────────────────────────
        # The mapper looks up sequences as orfs_fa["{orf_id}--{study_id}"] and
        # compares them against translated transcript sequences, so the FASTA
        # must contain amino-acid (protein) sequences.
        # Uses canonical translate_nt from orfont.core.utils (or built-in fallback).

        if os.path.exists(fasta_file):
            print(f"Using existing protein FASTA (skipping generation): {fasta_file}", file=sys.stderr)
        elif os.path.exists(metadata_file):
            print(f"Generating protein FASTA for GENCODE mapper from {metadata_file}...", file=sys.stderr)
            with _open(metadata_file) as mf, open(fasta_file, 'w') as ff:
                hdr = mf.readline().strip().split('\t')
                id_idx = hdr.index('orf_id')
                try:
                    seq_idx = hdr.index('sequence')
                except ValueError:
                    print("Error: Metadata missing sequence column", file=sys.stderr)
                    sys.exit(1)
                for line in mf:
                    parts = line.strip().split('\t')
                    if len(parts) <= max(id_idx, seq_idx):
                        continue
                    orf_id   = parts[id_idx]
                    nt_seq   = parts[seq_idx]
                    study_id = orf_to_study.get(orf_id, 'unified')
                    aa_seq   = translate_nt(nt_seq)
                    # Key format expected by the mapper: "{orf_id}--{study_id}"
                    ff.write(f">{orf_id}--{study_id}\n{aa_seq}\n")

        cmd = [
            "python3", script,
            "-d", args.ensembl_dir,
            "-f", fasta_file,
            "-b", bed6_file,
            "-o", os.path.join(args.output_dir, "gencode_results")
        ]
        run_command(cmd, "GENCODE Classifier")

    elif args.mode == 'orfquant':
        if not args.gtf:
            print("Error: --gtf is required for orfquant mode", file=sys.stderr)
            sys.exit(1)
            
        script = get_script_path("run_orfquant_classify.R")
        
        # Input should be the unified GTF (may be .gtf.gz)
        input_gtf = _resolve(input_base, '.gtf') if not (input_base.endswith('.gtf') or input_base.endswith('.gtf.gz')) else input_base
        
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
            "--output_prefix", os.path.join(args.output_dir, "orftype_results"),
            "--cpus", args.cpus
        ]
        run_command(cmd, "ORFtype Classifier")

if __name__ == "__main__":
    main(sys.argv[1:])
