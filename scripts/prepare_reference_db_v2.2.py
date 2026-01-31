#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Prepare Reference Database for Ribosome Profiling (Ribo-seq) Analysis
(Updated Version: v2.2)

Updates in v2.2:
- Added 'reference' directory output for alignment usage.
- Added Genome FASTA downloads to SPECIES_CONFIG.
- Automatically saves consistent Genome FASTA, Transcripts, and GTF to the 'reference' folder.

Key Features:
1.  **Species-Specific Logic**: Handles Animals (Human/Mouse) and Plants.
2.  **Layered Output**: Generates separate FASTA files for Nuclear rRNA, Mito rRNA, Chloroplast, tRNA, etc.
3.  **Reference Management**: Downloads and organizes Genome FASTA, Transcripts, and GTF for aligners.
4.  **Microbial Contamination**: Automatically extracts common bacterial rRNA.

Dependencies:
    - seqkit, wget, gzip, grep, awk, requests
"""

import os
import sys
import argparse
import subprocess
import logging
import shutil
import gzip
import requests
from pathlib import Path

# --- Configuration ---

logging.basicConfig(
    level=logging.INFO,
    format='[%(asctime)s] %(levelname)s: %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)

# SILVA Database URLs (Release 138.2)
SILVA_URLS = {
    "ssu": "https://www.arb-silva.de/fileadmin/silva_databases/release_138_2/Exports/SILVA_138.2_SSURef_NR99_tax_silva.fasta.gz",
    "lsu": "https://www.arb-silva.de/fileadmin/silva_databases/release_138_2/Exports/SILVA_138.2_LSURef_NR99_tax_silva.fasta.gz"
}

# Common Microbial Contaminants
MICROBIAL_CONTAMINANTS = [
    "Escherichia coli",
    "Staphylococcus aureus",
    "Bacillus subtilis",
    "Mycoplasma",
    "Pseudomonas aeruginosa"
]

# Species Configuration (Added 'genome_url' for Reference)
SPECIES_CONFIG = {
    "human": {
        "tax_name": "Homo sapiens",
        # GENCODE Release 49
        "gtf_url": "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.annotation.gtf.gz",
        "transcripts_url": "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.transcripts.fa.gz",
        "lnc_url": "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.lncRNA_transcripts.fa.gz",
        "genome_url": "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/GRCh38.primary_assembly.genome.fa.gz",
        "gtrnadb_url": "https://gtrnadb.org/genomes/eukaryota/Hsapi38/hg38-mature-tRNAs.fa"
    },
    "mouse": {
        "tax_name": "Mus musculus",
        # GENCODE Release M38
        "gtf_url": "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M38/gencode.vM38.annotation.gtf.gz",
        "transcripts_url": "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M38/gencode.vM38.transcripts.fa.gz",
        "lnc_url": "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M38/gencode.vM38.lncRNA_transcripts.fa.gz",
        "genome_url": "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M38/GRCm39.primary_assembly.genome.fa.gz",
        "gtrnadb_url": "https://gtrnadb.org/genomes/eukaryota/Mmusc39/mm39-mature-tRNAs.fa"
    },
    # --- Plants (Ensembl Plants Release 62) ---
    "rice": {
        "tax_name": "Oryza sativa",
        "gtf_url": "https://ftp.ensemblgenomes.org/pub/plants/release-62/gtf/oryza_sativa/Oryza_sativa.IRGSP-1.0.62.gtf.gz",
        "transcripts_url": "https://ftp.ensemblgenomes.org/pub/plants/release-62/fasta/oryza_sativa/cdna/Oryza_sativa.IRGSP-1.0.cdna.all.fa.gz",
        "lnc_url": "https://ftp.ensemblgenomes.org/pub/plants/release-62/fasta/oryza_sativa/ncrna/Oryza_sativa.IRGSP-1.0.ncrna.fa.gz",
        "genome_url": "https://ftp.ensemblgenomes.org/pub/plants/release-62/fasta/oryza_sativa/dna/Oryza_sativa.IRGSP-1.0.dna.toplevel.fa.gz",
        "gtrnadb_url": "https://gtrnadb.org/genomes/eukaryota/Osati7/orySat7-mature-tRNAs.fa",
        "chloroplast_acc": "NC_001320"
    },
    "maize": {
        "tax_name": "Zea mays",
        "gtf_url": "https://ftp.ensemblgenomes.org/pub/plants/release-62/gtf/zea_mays/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.62.gtf.gz",
        "transcripts_url": "https://ftp.ensemblgenomes.org/pub/plants/release-62/fasta/zea_mays/cdna/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.cdna.all.fa.gz",
        "lnc_url": "https://ftp.ensemblgenomes.org/pub/plants/release-62/fasta/zea_mays/ncrna/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.ncrna.fa.gz",
        "genome_url": "https://ftp.ensemblgenomes.org/pub/plants/release-62/fasta/zea_mays/dna/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.dna.toplevel.fa.gz",
        "gtrnadb_url": "https://gtrnadb.org/genomes/eukaryota/Zmays8/zeaMay8-mature-tRNAs.fa",
        "chloroplast_acc": "NC_001666"
    },
    "wheat": {
        "tax_name": "Triticum aestivum",
        "gtf_url": "https://ftp.ensemblgenomes.org/pub/plants/release-62/gtf/triticum_aestivum/Triticum_aestivum.IWGSC.62.gtf.gz",
        "transcripts_url": "https://ftp.ensemblgenomes.org/pub/plants/release-62/fasta/triticum_aestivum/cdna/Triticum_aestivum.IWGSC.cdna.all.fa.gz",
        "lnc_url": "https://ftp.ensemblgenomes.org/pub/plants/release-62/fasta/triticum_aestivum/ncrna/Triticum_aestivum.IWGSC.ncrna.fa.gz",
        "genome_url": "https://ftp.ensemblgenomes.org/pub/plants/release-62/fasta/triticum_aestivum/dna/Triticum_aestivum.IWGSC.dna.toplevel.fa.gz",
        "gtrnadb_url": "https://gtrnadb.org/genomes/eukaryota/Taest2/triAes2-mature-tRNAs.fa",
        "chloroplast_acc": "NC_002762"
    },
    "soybean": {
        "tax_name": "Glycine max",
        "gtf_url": "https://ftp.ensemblgenomes.org/pub/plants/release-62/gtf/glycine_max/Glycine_max.Glycine_max_v2.1.62.gtf.gz",
        "transcripts_url": "https://ftp.ensemblgenomes.org/pub/plants/release-62/fasta/glycine_max/cdna/Glycine_max.Glycine_max_v2.1.cdna.all.fa.gz",
        "lnc_url": "https://ftp.ensemblgenomes.org/pub/plants/release-62/fasta/glycine_max/ncrna/Glycine_max.Glycine_max_v2.1.ncrna.fa.gz",
        "genome_url": "https://ftp.ensemblgenomes.org/pub/plants/release-62/fasta/glycine_max/dna/Glycine_max.Glycine_max_v2.1.dna.toplevel.fa.gz",
        "gtrnadb_url": "https://gtrnadb.org/genomes/eukaryota/Gmax2.1/glyMax2.1-mature-tRNAs.fa",
        "chloroplast_acc": "NC_007942"
    }
}

# --- Helper Functions ---

def run_cmd(cmd, shell=False, check=True):
    try:
        if shell:
            subprocess.run(cmd, shell=True, check=check, executable='/bin/bash')
        else:
            subprocess.run(cmd, check=check)
    except subprocess.CalledProcessError as e:
        logger.error(f"Command failed: {cmd}")
        raise e

def ensure_dir(path):
    Path(path).mkdir(parents=True, exist_ok=True)

def parse_gtf_attributes(attr_string):
    attributes = {}
    if not attr_string: return attributes
    for block in attr_string.strip().split(';'):
        block = block.strip()
        if not block: continue
        if '=' in block and ' ' not in block: key, value = block.split('=', 1)
        else:
            if ' ' not in block: continue
            key, value = block.split(' ', 1)
        attributes[key] = value.strip().strip('"')
    return attributes

def download_file_wget(url, output_file):
    """Downloads a file using wget if it doesn't exist."""
    if not os.path.exists(output_file):
        logger.info(f"Downloading {url}...")
        # Using wget allows easy handling of FTP/HTTP and progress bars
        run_cmd(["wget", "--no-check-certificate", "-q", "-O", output_file, url])
    else:
        logger.info(f"File already exists (skipping download): {output_file}")

# --- Core Logic ---

def download_silva(ref_dir):
    ensure_dir(ref_dir)
    for key, url in SILVA_URLS.items():
        filename = f"silva_{key}.fasta.gz"
        filepath = os.path.join(ref_dir, filename)
        download_file_wget(url, filepath)

def extract_silva_taxa(tax_names, output_file, ref_dir):
    ssu_path = os.path.join(ref_dir, "silva_ssu.fasta.gz")
    lsu_path = os.path.join(ref_dir, "silva_lsu.fasta.gz")
    pattern = "|".join(tax_names)
    logger.info(f"Extracting SILVA sequences matching: {pattern}")
    cmd_ssu = f"seqkit grep -n -r -p \"{pattern}\" \"{ssu_path}\" -o \"{output_file}\""
    cmd_lsu = f"seqkit grep -n -r -p \"{pattern}\" \"{lsu_path}\" >> \"{output_file}\""
    run_cmd(cmd_ssu, shell=True)
    run_cmd(cmd_lsu, shell=True)

def create_microbial_contamination_db(ref_dir, contam_dir):
    output_file = os.path.join(contam_dir, "common_microbial_rRNA.fasta")
    if os.path.exists(output_file) and os.path.getsize(output_file) > 0:
        logger.info("Microbial contamination DB already exists.")
        return output_file
    logger.info("Creating Common Microbial Contamination DB...")
    extract_silva_taxa(MICROBIAL_CONTAMINANTS, output_file, ref_dir)
    return output_file

def download_chloroplast_ncbi(accession, species, contam_dir):
    output_file = os.path.join(contam_dir, f"{species}_chloroplast_genome.fasta")
    if os.path.exists(output_file) and os.path.getsize(output_file) > 0:
        return output_file
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={accession}&rettype=fasta&retmode=text"
    logger.info(f"[{species}] Downloading Chloroplast Genome {accession}...")
    try:
        run_cmd(["wget", "--no-check-certificate", "-q", "-O", output_file, url])
        return output_file
    except Exception as e:
        logger.error(f"[{species}] Failed to download chloroplast: {e}")
        return None

def fetch_external_trna(species, contam_dir):
    """Downloads GtRNAdb using python requests to avoid 403 Forbidden."""
    url = SPECIES_CONFIG[species].get('gtrnadb_url')
    if not url: return None
    
    out_file = os.path.join(contam_dir, f"{species}_gtrnadb_tRNA.fasta")
    
    # Check for valid existing file
    if os.path.exists(out_file):
         with open(out_file, 'r') as f:
            line = f.readline()
            if not line or not line.startswith('>'):
                logger.warning(f"Found invalid GtRNAdb file, removing: {out_file}")
                os.remove(out_file)

    if not os.path.exists(out_file):
        logger.info(f"[{species}] Downloading GtRNAdb...")
        headers = {'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36'}
        try:
            r = requests.get(url, headers=headers, stream=True, timeout=30)
            r.raise_for_status()
            with open(out_file, 'wb') as f:
                for chunk in r.iter_content(chunk_size=8192):
                    f.write(chunk)
            # Verify
            with open(out_file, 'r') as f:
                if not f.readline().startswith('>'): raise ValueError("Content is not FASTA")
        except Exception as e:
            logger.error(f"Failed to download {url}: {e}")
            if os.path.exists(out_file): os.remove(out_file)
            raise e
    return out_file

def prepare_alignment_references(species, ref_dir, final_ref_dir):
    """
    Downloads Genome FASTA, Transcripts, and GTF to staging,
    then copies them to the final 'reference' directory for alignment usage.
    """
    conf = SPECIES_CONFIG[species]
    logger.info(f"[{species}] Preparing alignment references...")

    # Define Local Staging Paths (Source)
    gtf_local = os.path.join(ref_dir, f"{species}.gtf.gz")
    cdna_local = os.path.join(ref_dir, f"{species}.transcripts.fa.gz")
    genome_local = os.path.join(ref_dir, f"{species}.genome.fa.gz")

    # Define Final Reference Paths (Destination)
    final_gtf = os.path.join(final_ref_dir, f"{species}.gtf.gz")
    final_cdna = os.path.join(final_ref_dir, f"{species}.transcripts.fa.gz")
    final_genome = os.path.join(final_ref_dir, f"{species}.genome.fa.gz")

    # 1. GTF
    download_file_wget(conf['gtf_url'], gtf_local)
    if not os.path.exists(final_gtf):
        shutil.copy(gtf_local, final_gtf)

    # 2. Transcripts (cDNA)
    download_file_wget(conf['transcripts_url'], cdna_local)
    if not os.path.exists(final_cdna):
        shutil.copy(cdna_local, final_cdna)

    # 3. Genome FASTA
    if conf.get('genome_url'):
        download_file_wget(conf['genome_url'], genome_local)
        if not os.path.exists(final_genome):
            shutil.copy(genome_local, final_genome)
    else:
        logger.warning(f"[{species}] No genome_url configured! Skipping genome download.")

def process_genome_annotation_split(species, ref_dir, contam_dir):
    """
    Parses GTF and splits extraction into separate files.
    Note: Files are assumed to be downloaded by prepare_alignment_references already.
    """
    conf = SPECIES_CONFIG[species]
    
    # Paths must match what prepare_alignment_references uses
    gtf_local = os.path.join(ref_dir, f"{species}.gtf.gz")
    cdna_local = os.path.join(ref_dir, f"{species}.transcripts.fa.gz")
    
    # lncRNA is optional and might not be in the alignment ref set, so we handle it here
    ncrna_local = os.path.join(ref_dir, f"{species}.ncrna.fa.gz") if conf.get('lnc_url') else None
    if ncrna_local: download_file_wget(conf['lnc_url'], ncrna_local)

    # 2. Categorize IDs
    categories = {
        'nuclear_rRNA': {'rrna', 'ribozyme', 'srna'},
        'mito_rRNA': {'mt_rrna'},
        'chloroplast_rRNA': {'chloroplast_rrna'}, 
        'tRNA': {'trna', 'mt_trna', 'chloroplast_trna'},
        'other_ncRNA': {'snrna', 'snorna', 'scrna', 'scarna', 'misc_rna'}
    }
    category_ids = {k: set() for k in categories}
    
    logger.info(f"[{species}] Parsing GTF for contamination sequences...")
    try:
        with gzip.open(gtf_local, 'rt') as f:
            for line in f:
                if line.startswith("#"): continue
                parts = line.strip().split('\t')
                if len(parts) < 9: continue
                if parts[2] not in ['transcript', 'exon', 'gene']: continue
                
                attr = parse_gtf_attributes(parts[8])
                biotype = None
                for key in ['transcript_type', 'gene_type', 'transcript_biotype', 'gene_biotype', 'biotype']:
                    if key in attr: biotype = attr[key].lower(); break
                
                if not biotype: continue
                t_id = attr.get('transcript_id')
                if not t_id: continue

                for cat, keywords in categories.items():
                    if biotype in keywords or any(k in biotype for k in keywords):
                        category_ids[cat].add(t_id)
                        break
    except Exception as e:
        logger.error(f"Error parsing GTF: {e}")
        return {}

    # 3. Clean FASTA headers (Removed | and everything after)
    fasta_sources = [cdna_local]
    if ncrna_local and os.path.exists(ncrna_local): fasta_sources.append(ncrna_local)
    
    cleaned_fasta = os.path.join(ref_dir, f"{species}.all_transcripts.clean.fa")
    # Using double backslash \\| to ensure it works in shell via python
    cmd_cat = f"cat {' '.join(fasta_sources)} | seqkit replace -p \"\\|.*\" -r \"\" > \"{cleaned_fasta}\""
    run_cmd(cmd_cat, shell=True)

    # 4. Extract per category
    result_files = {}
    for cat, ids in category_ids.items():
        if not ids: continue
        id_file = os.path.join(ref_dir, f"{species}_{cat}_ids.txt")
        with open(id_file, 'w') as f:
            for i in ids: f.write(f"{i}\n")
            
        out_file = os.path.join(contam_dir, f"{species}_genome_{cat}.fasta")
        if not os.path.exists(out_file):
            cmd_extract = f"seqkit grep -f \"{id_file}\" \"{cleaned_fasta}\" -o \"{out_file}\""
            run_cmd(cmd_extract, shell=True)
        
        if os.path.exists(out_file) and os.path.getsize(out_file) > 0:
            result_files[cat] = out_file
            
    if os.path.exists(cleaned_fasta): os.remove(cleaned_fasta)
    return result_files

def tag_and_merge(file_map, output_file):
    temp_files = []
    for source, fpath in file_map:
        if not fpath or not os.path.exists(fpath): continue
        tagged_file = fpath + ".tagged"
        cmd = f"seqkit replace -p \"^>\" -r \">{source}|\" \"{fpath}\" -o \"{tagged_file}\""
        run_cmd(cmd, shell=True)
        temp_files.append(tagged_file)
    
    if not temp_files:
        logger.warning("No files to merge!")
        return

    logger.info(f"Merging {len(temp_files)} files into {output_file}...")
    run_cmd(f"cat {' '.join(temp_files)} > \"{output_file}\"", shell=True)
    
    dedup_file = output_file.replace(".fasta", ".dedup.fasta")
    run_cmd(f"seqkit rmdup -s -i \"{output_file}\" -o \"{dedup_file}\"", shell=True)
    os.rename(dedup_file, output_file)

    for f in temp_files:
        if os.path.exists(f): os.remove(f)

def main():
    parser = argparse.ArgumentParser(description="Prepare Ribosome Profiling Reference & Contamination DB")
    parser.add_argument("-o", "--output-dir", default="./reference_data_project", help="Base directory")
    parser.add_argument("-s", "--species", nargs="+", default=list(SPECIES_CONFIG.keys()))
    args = parser.parse_args()
    
    check_dependencies = ['seqkit', 'wget', 'gzip', 'grep']
    for tool in check_dependencies:
        if not shutil.which(tool):
            logger.error(f"Missing tool: {tool}"); sys.exit(1)

    base_dir = Path(args.output_dir)
    ref_dir = base_dir / "reference_source"       # Raw downloads / staging
    contam_dir = base_dir / "contamination_indices" # Final Contamination DB
    final_ref_dir = base_dir / "reference"        # Final Alignment References (Genome/GTF)
    
    ensure_dir(ref_dir)
    ensure_dir(contam_dir)
    ensure_dir(final_ref_dir)

    # 1. Prepare Shared Resources
    download_silva(str(ref_dir))
    microbial_db = create_microbial_contamination_db(str(ref_dir), str(contam_dir))

    # 2. Process Species
    for sp in args.species:
        if sp not in SPECIES_CONFIG: continue
        logger.info(f"=== Processing {sp} ===")
        
        # --- NEW: Prepare Alignment References (Genome, GTF, Transcripts) ---
        # This downloads files to ref_dir AND copies them to the 'reference' folder
        prepare_alignment_references(sp, str(ref_dir), str(final_ref_dir))

        files_to_merge = []

        # A. SILVA specific rRNA
        silva_out = os.path.join(contam_dir, f"{sp}_silva_rRNA.fasta")
        extract_silva_taxa([SPECIES_CONFIG[sp]['tax_name']], silva_out, str(ref_dir))
        files_to_merge.append(("SILVA_rRNA", silva_out))

        # B. Genome Annotations (Splitting logic uses files downloaded in step above)
        genome_files = process_genome_annotation_split(sp, str(ref_dir), str(contam_dir))
        for cat, fpath in genome_files.items():
            files_to_merge.append((f"GENOME_{cat}", fpath))

        # C. Independent Chloroplast
        if 'chloroplast_acc' in SPECIES_CONFIG[sp]:
            chloro_file = download_chloroplast_ncbi(SPECIES_CONFIG[sp]['chloroplast_acc'], sp, str(contam_dir))
            if chloro_file:
                files_to_merge.append(("NCBI_Chloroplast_Genome", chloro_file))

        # D. External tRNA
        trna_file = fetch_external_trna(sp, str(contam_dir))
        if trna_file:
            files_to_merge.append(("GtRNAdb", trna_file))

        # E. Microbial Contaminants
        if microbial_db:
            files_to_merge.append(("Ext_Microbial", microbial_db))

        # F. Final Merge
        final_out = os.path.join(contam_dir, f"{sp}_final_contamination.fasta")
        tag_and_merge(files_to_merge, final_out)

        logger.info(f"=== Finished {sp}. ===")
        logger.info(f"Alignment Refs: {final_ref_dir}")
        logger.info(f"Contamination DB: {final_out}")

if __name__ == "__main__":
    main()