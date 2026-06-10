#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Prepare Workflow for Ribo-seq Analysis
========================================

This script sets up a complete working environment for nf-core/riboseq pipeline:
1. Creates standardized directory structure (data, reference, containers, process, result)
2. Creates symbolic links for data files and container images
3. Generates sample sheet CSV using get_sample_sheet.py
4. Prepares Nextflow execution shell script with recommended parameters

Author: Automated Workflow Setup
Date: 2026-01-31
Version: 1.0

Dependencies:
    - Python 3.6+
    - scripts/get_sample_sheet.py (for sample sheet generation)
"""

import os
import sys
import argparse
import subprocess
import json
import logging
from pathlib import Path
from datetime import datetime
import shutil

# --- Configuration ---

logging.basicConfig(
    level=logging.INFO,
    format='[%(asctime)s] %(levelname)s: %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)

# Directory structure
STANDARD_DIRS = {
    'data': 'Raw sequencing data (FASTQ files)',
    'reference': 'Reference genome, GTF, and FASTA files',
    'containers': 'Singularity/Docker container images',
    'process': 'Nextflow work directory and logs',
    'result': 'Pipeline output results',
    'scripts': 'Analysis scripts and configurations'
}

# Container mapping
DEFAULT_CONTAINERS = {
    'orfquant': 'orfquant_patched.sif',
    'rpbp': 'rpbp.sif',
    'python': 'python_3.9.sif'
}

# Reference files mapping
# Patterns based on prepare_reference_db_v2.2.py output:
#   - reference/: {species}.genome.fa.gz, {species}.gtf.gz, {species}.transcripts.fa.gz
#   - contamination_indices/: {species}_final_contamination.fasta
REFERENCE_FILES = {
    'genome_fasta': [
        '*.genome.fa', '*.genome.fasta', '*.genome.fa.gz', '*.genome.fasta.gz',  # species.genome.fa
        '*.dna.toplevel.fa*', '*.dna.primary_assembly.fa*',  # Ensembl naming
        '*primary_assembly.genome.fa*',  # GENCODE naming (GRCh38.primary_assembly.genome.fa.gz)
    ],
    'gtf': ['*.gtf', '*.gtf.gz'],
    'transcripts': [
        '*.transcripts.fa', '*.transcripts.fa.gz',  # species.transcripts.fa
        '*.cdna.all.fa*', '*.cdna.fa*',  # Ensembl naming
    ],
    'contaminant': [
        '*_final_contamination.fasta', '*_final_contamination.fa',  # Primary: from prepare_reference_db
        '*contamination*.fa*', '*contaminant*.fa*',
    ]
}

ENSEMBL_SPECIES_MAP = {
    'human': 'homo_sapiens',
    'mouse': 'mus_musculus',
    'rice': 'oryza_sativa',
    'maize': 'zea_mays',
    'wheat': 'triticum_aestivum',
    'soybean': 'glycine_max',
    'arabidopsis': 'arabidopsis_thaliana',
}

ENSEMBL_DIVISION_MAP = {
    'human': 'ensembl',
    'mouse': 'ensembl',
    'rice': 'plants',
    'maize': 'plants',
    'wheat': 'plants',
    'soybean': 'plants',
    'arabidopsis': 'plants',
}

DEFAULT_ENSEMBL_RELEASE = {
    'ensembl': 110,
    'plants': 58,
}

DEFAULT_ENSEMBL_ASSEMBLY = {
    'human':   'GRCh38',
    'mouse':   'GRCm39',
    'rice':    'IRGSP-1.0',
    'maize':   'Zm-B73-REFERENCE-NAM-5.0',
    'wheat':   'IWGSC',
    'soybean': 'Soybean_JD17',
    'arabidopsis': 'TAIR10',
}


def parse_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description="Prepare complete working environment for Ribo-seq analysis",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic setup with data directory
  %(prog)s -w /path/to/workdir -d /path/to/fastq_dir

  # Full setup with references and containers
  %(prog)s -w /path/to/workdir \\
           -d /path/to/fastq_dir \\
           -r /path/to/reference_dir \\
           --contaminant-dir /path/to/contamination_indices \\
           --species human

  # Using prepare_reference_db_v2.2.py output structure
  %(prog)s -w /path/to/workdir \\
           -d /path/to/fastq_dir \\
           -r /path/to/reference_data_project/reference \\
           --contaminant-dir /path/to/reference_data_project/contamination_indices \\
           --species rice

  # Setup with custom containers
  %(prog)s -w /path/to/workdir \\
           -d /path/to/fastq_dir \\
           --orfquant-container /path/to/orfquant_patched.sif \\
           --rpbp-container /path/to/rpbp.sif
        """
    )
    
    # Required arguments
    parser.add_argument('-w', '--workdir', required=True,
                        help='Working directory (will be created if not exists)')
    parser.add_argument('-d', '--data-dir', required=True,
                        help='Directory containing FASTQ files')
    
    # Optional data sources
    parser.add_argument('-r', '--reference-dir', default=None,
                        help='Directory containing reference files (genome FASTA, GTF, transcripts). '
                             'Expected files: {species}.genome.fa, {species}.gtf, {species}.transcripts.fa. '
                             'If not provided, prepare_reference_db_v2.2.py will be used to generate references.')
    parser.add_argument('--contaminant-dir', default=None,
                        help='Directory containing contamination database files. '
                             'Expected file: {species}_final_contamination.fasta. '
                             'If not provided, will search in reference-dir or auto-generate.')
    parser.add_argument('-c', '--container-dir', default=None,
                        help='(Deprecated) Directory containing Singularity container images. '
                             'Prefer specifying containers directly with --orfquant-container and --rpbp-container.')
    
    # Reference genome settings
    parser.add_argument('--genome', default='GRCh38',
                        help='Genome name (default: GRCh38). Example: GRCh38/GRCm39/IRGSP-1.0')
    parser.add_argument('--species', default='human',
                        choices=['human', 'mouse', 'rice', 'maize', 'wheat', 'soybean', 'arabidopsis'],
                        help='Species name (default: human). Choices: human, mouse, rice, maize, wheat, soybean')

    # ORF classification Ensembl annotation
    parser.add_argument('--orf-classify-ensembl-dir', default=None,
                        help='Existing Ensembl annotation directory for ORF classification')
    parser.add_argument('--ensembl-release', type=int, default=None,
                        help='Ensembl release for ORF classification (default depends on species)')
    parser.add_argument('--ensembl-assembly', default=None,
                        help='Genome assembly for Ensembl annotation (default: --genome)')
    parser.add_argument('--ensembl-version-suffix', default=None,
                        help='Ensembl GTF version suffix override (e.g., 110.38)')
    parser.add_argument('--ensembl-base-url', default=None,
                        help='Override Ensembl FTP base URL')
    parser.add_argument('--skip-orf-classify-ensembl', action='store_false',
                        dest='prepare_orf_classify_ensembl',
                        help='Skip preparing Ensembl files for ORF classification (default: prepare)')
    parser.set_defaults(prepare_orf_classify_ensembl=True)
    
    # Sample sheet options
    parser.add_argument('--strandedness', default='auto',
                        help='Strandedness for sample sheet (default: auto)')
    parser.add_argument('--sample-type', default='riboseq',
                        help='Sample type (riboseq/tiseq) (default: riboseq)')
    parser.add_argument('--group-map', default=None,
                        help=(
                            'Optional path to a JSON or two-column CSV file mapping sample IDs to '
                            'group names, used to populate the "group" column in the samplesheet. '
                            'JSON format: {"sample_id": "group_name", ...}. '
                            'CSV format (no header): sample_id,group_name. '
                            'Samples not listed in the map will have an empty group.'
                        ))
    parser.add_argument('--merge-replicates', action='store_true', default=False,
                        help='Add --merge_replicates to the generated Nextflow run script (requires --group-map to be useful)')
    
    # Ribotricer options
    parser.add_argument('--ribotricer-phase-score-cutoff', type=float, default=None,
                        help='Ribotricer phase score cutoff (default: None uses ribotricer built-in default ~0.429). '
                             'Lower this (e.g. 0.1) for plant/low-depth samples with weak metagene periodicity.')
    
    # sORF filter options
    parser.add_argument('--sorf-unique-mode', default=None,
                        choices=['auto', 'nh', 'mapq'],
                        help='Mode for unique-mapping filter: auto (detect NH tag), nh, or mapq (default: auto)')
    
    # QC options
    parser.add_argument('--skip-collect-qc-stats', action='store_true', default=False,
                        help='Skip the COLLECT_QC_STATS step that aggregates per-sample QC metrics')
    
    # Container options
    parser.add_argument('--orfquant-container', default=None,
                        help='Path to ORFquant patched container')
    parser.add_argument('--rpbp-container', default=None,
                        help='Path to RPBP container')
    parser.add_argument('--unify-orf-container', default=None,
                        help='Path to container for unify_orf_predictions and classify_orfs (Python/biopython)')
    parser.add_argument('--gencode-orf-mapper-container', default=None,
                        help='Path to gencode-orf-mapper container for GENCODE ORF classification')
    parser.add_argument('--ribowaltz-container', default=None,
                        help='Path to riboWaltz container for P-site analysis')
    
    # Pipeline options
    parser.add_argument('--run-prefilter-qc', action='store_true',
                        help='Enable prefilter QC for comparison')
    parser.add_argument('--unify-orf-min-len', type=int, default=6,
                        help='Minimum ORF length in amino acids for unification (default: 6)')
    parser.add_argument('--profile', default='singularity',
                        help='Nextflow profile (default: singularity)')
    
    # P-site correction options
    parser.add_argument('--orfquant-psite-correction', action='store_true', default=True,
                        help='Enable ORFquant P-site offset correction (default: enabled)')
    parser.add_argument('--no-orfquant-psite-correction', action='store_false', dest='orfquant_psite_correction',
                        help='Disable ORFquant P-site offset correction')
    
    # Advanced ORF unification options
    parser.add_argument('--unify-orf-frame-merge-min-overlap', type=float, default=0.9,
                        help='Minimum overlap fraction of shorter ORF for single-exon frame-aware merging (default: 0.9)')
    parser.add_argument('--unify-orf-no-frame-merge', action='store_true',
                        help='Disable frame-aware merging (only use exact matches)')
    parser.add_argument('--unify-orf-seq-cluster', action='store_true',
                        help='Enable Stage 3 sequence-similarity clustering (disabled by default)')
    
    # SRA conversion options
    parser.add_argument('--sra-threads', type=int, default=8,
                        help='Number of threads for fasterq-dump when converting SRA files (default: 8)')
    parser.add_argument('--pigz-threads', type=int, default=8,
                        help='Number of threads for pigz compression when converting SRA files (default: 8)')
    
    # Resource allocation options
    parser.add_argument('--aligner', default='star',
                        choices=['star', 'hisat2'],
                        help='Genome alignment tool to use (default: star)')
    parser.add_argument('--max-memory', default='150.GB',
                        help='Maximum memory for Nextflow (default: 150.GB)')
    parser.add_argument('--max-cpus', type=int, default=42,
                        help='Maximum CPUs for Nextflow (default: 42)')
    parser.add_argument('--max-time', default='48.h',
                        help='Maximum time for Nextflow (default: 48.h)')
    parser.add_argument('--save-reference', action='store_true', default=True,
                        help='Save reference files for reuse (default: True)')
    parser.add_argument('--no-save-reference', action='store_false', dest='save_reference',
                        help='Do not save reference files')
    parser.add_argument('--bg', action='store_true', default=False,
                        help='Run Nextflow in background mode with -bg (default: False, run interactively)')
    
    # Script generation
    parser.add_argument('--script-name', default='run_pipeline.sh',
                        help='Name of generated execution script (default: run_pipeline.sh)')
    parser.add_argument('--dry-run', action='store_true',
                        help='Dry run mode (show actions without executing)')
    
    return parser.parse_args()


def create_directory_structure(workdir, dry_run=False):
    """Create standardized directory structure"""
    logger.info("=" * 60)
    logger.info("Creating directory structure")
    logger.info("=" * 60)
    
    workdir = Path(workdir).resolve()
    created_dirs = []
    
    for dir_name, description in STANDARD_DIRS.items():
        dir_path = workdir / dir_name
        
        if dry_run:
            logger.info(f"[DRY RUN] Would create: {dir_path}")
            logger.info(f"          Purpose: {description}")
        else:
            dir_path.mkdir(parents=True, exist_ok=True)
            logger.info(f"✓ Created: {dir_path}")
            logger.info(f"  Purpose: {description}")
            created_dirs.append(dir_path)
    
    return created_dirs


def create_symlinks(source_dir, target_dir, patterns=None, dry_run=False):
    """Create symbolic links from source to target directory"""
    source_dir = Path(source_dir).resolve()
    target_dir = Path(target_dir).resolve()
    
    if not source_dir.exists():
        logger.warning(f"Source directory does not exist: {source_dir}")
        return []
    
    linked_files = []
    
    if patterns:
        # Search by patterns (use ** for recursive search)
        for pattern in patterns:
            # Try both direct and recursive search
            files = list(source_dir.glob(pattern))
            # Also try recursive search if no direct matches
            if not files:
                files = list(source_dir.glob(f"**/{pattern}"))
            
            for file_path in files:
                if file_path.is_file():
                    link_path = target_dir / file_path.name
                    
                    if dry_run:
                        logger.info(f"[DRY RUN] Would link: {file_path} -> {link_path}")
                        linked_files.append(link_path)
                    else:
                        if link_path.exists() or link_path.is_symlink():
                            logger.info(f"  Skipping (exists): {link_path.name}")
                            linked_files.append(link_path)  # Still track existing files
                        else:
                            link_path.symlink_to(file_path)
                            logger.info(f"✓ Linked: {file_path.name}")
                            linked_files.append(link_path)
    else:
        # Link all files
        for file_path in source_dir.iterdir():
            if file_path.is_file():
                link_path = target_dir / file_path.name
                
                if dry_run:
                    logger.info(f"[DRY RUN] Would link: {file_path} -> {link_path}")
                    linked_files.append(link_path)
                else:
                    if link_path.exists() or link_path.is_symlink():
                        logger.info(f"  Skipping (exists): {link_path.name}")
                        linked_files.append(link_path)  # Still track existing files
                    else:
                        link_path.symlink_to(file_path)
                        logger.info(f"✓ Linked: {file_path.name}")
                        linked_files.append(link_path)
    
    return linked_files


def convert_sra_to_fastq(source_dir, target_dir, dry_run=False, threads_dump=8, threads_pigz=8):
    """Convert SRA files to FASTQ.gz using scripts/sra2fq.sh"""
    source_dir = Path(source_dir).resolve()
    target_dir = Path(target_dir).resolve()

    sra_files = list(source_dir.glob('*.sra'))
    if not sra_files:
        return []

    # Find sra2fq.sh script
    script_dir = Path(__file__).parent
    sra2fq_script = script_dir / 'sra2fq.sh'

    if not sra2fq_script.exists():
        logger.error(f"sra2fq.sh not found at: {sra2fq_script}")
        return []

    # Build command: sra2fq.sh -t <threads> -p <pigz_threads> -o <output_dir> <sra_files...>
    cmd = [
        'bash',
        str(sra2fq_script),
        '-t', str(threads_dump),
        '-p', str(threads_pigz),
        '-o', str(target_dir)
    ] + [str(f) for f in sra_files]

    if dry_run:
        logger.info(f"[DRY RUN] Would convert SRA: {' '.join(cmd)}")
        return [target_dir / f"{f.stem}.fastq.gz" for f in sra_files]

    try:
        logger.info(f"Converting {len(sra_files)} SRA file(s) to FASTQ.gz...")
        subprocess.run(cmd, check=True)
        generated = list(target_dir.glob('*.fastq.gz')) + list(target_dir.glob('*.fq.gz'))
        logger.info(f"Generated {len(generated)} FASTQ.gz file(s)")
        return generated
    except subprocess.CalledProcessError as e:
        logger.error(f"Failed to convert SRA files: {e}")
        return []


def prepare_reference_db_if_missing(workdir, species, dry_run=False):
    """Prepare reference database using prepare_reference_db_v2.2.py when not provided"""
    script_dir = Path(__file__).parent
    prep_script = script_dir / 'prepare_reference_db_v2.2.py'
    if not prep_script.exists():
        logger.error(f"prepare_reference_db_v2.2.py not found at: {prep_script}")
        return None

    base_dir = Path(workdir).resolve() / 'reference_data_project'
    reference_dir = base_dir / 'reference'
    contaminant_dir = base_dir / 'contamination_indices'

    cmd = [
        'python3',
        str(prep_script),
        '-o', str(base_dir),
        '-s', str(species)
    ]

    if dry_run:
        logger.info(f"[DRY RUN] Would execute: {' '.join(cmd)}")
        return {
            'base_dir': base_dir,
            'reference_dir': reference_dir,
            'contaminant_dir': contaminant_dir
        }

    try:
        logger.info(f"Executing: {' '.join(cmd)}")
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        logger.error(f"Failed to prepare reference DB: {e}")
        return None

    return {
        'base_dir': base_dir,
        'reference_dir': reference_dir,
        'contaminant_dir': contaminant_dir
    }


def setup_data_directory(workdir, data_dir, dry_run=False, sra_threads=8, pigz_threads=8):
    """Setup data directory with symbolic links to FASTQ files"""
    logger.info("\n" + "=" * 60)
    logger.info("Setting up data directory")
    logger.info("=" * 60)
    
    target_dir = Path(workdir) / 'data'

    # Convert SRA files to FASTQ.gz using sra2fq.sh
    sra_converted = convert_sra_to_fastq(data_dir, target_dir, dry_run, 
                                          threads_dump=sra_threads, threads_pigz=pigz_threads)
    if sra_converted:
        logger.info(f"\nTotal FASTQ.gz files converted from SRA: {len(sra_converted)}")
    
    # Link FASTQ files
    fastq_patterns = ['*.fastq.gz', '*.fq.gz', '*.fastq', '*.fq']
    linked = create_symlinks(data_dir, target_dir, fastq_patterns, dry_run)

    # Keep only FASTQ files in workdir/data
    if not dry_run:
        for file_path in target_dir.iterdir():
            if file_path.is_file() and not file_path.name.endswith(('.fastq', '.fq', '.fastq.gz', '.fq.gz')):
                file_path.unlink()
                logger.info(f"  Removed non-FASTQ file: {file_path.name}")
    
    logger.info(f"\nTotal FASTQ files linked: {len(linked)}")
    return linked


def setup_reference_directory(workdir, reference_dirs, species, dry_run=False):
    """Setup reference directory with symbolic links for specified species"""
    logger.info("\n" + "=" * 60)
    logger.info(f"Setting up reference directory for species: {species}")
    logger.info("=" * 60)

    if not reference_dirs:
        logger.warning("No reference directory specified, skipping")
        return {}

    if isinstance(reference_dirs, (str, Path)):
        reference_dirs = [reference_dirs]
    else:
        reference_dirs = list(reference_dirs)
    
    target_dir = Path(workdir) / 'reference'
    linked_refs = {}
    
    for ref_type, patterns in REFERENCE_FILES.items():
        logger.info(f"\nSearching for {ref_type} files...")
        files = []
        for ref_dir in reference_dirs:
            all_files = create_symlinks(ref_dir, target_dir, patterns, dry_run)
            # Filter files to only include the specified species
            species_files = [f for f in all_files if species in f.name]
            files.extend(species_files)
        if files:
            linked_refs[ref_type] = files
            logger.info(f"  Found {len(files)} {ref_type} file(s) for {species}")

    # Decompress only the linked files (species-specific)
    if linked_refs:
        logger.info("\nDecompressing linked reference .gz files...")
        files_to_decompress = []
        for ref_type, file_list in linked_refs.items():
            files_to_decompress.extend(file_list)
        
        for gz_path in files_to_decompress:
            if str(gz_path).endswith('.gz'):
                out_path = Path(str(gz_path).rsplit('.gz', 1)[0])
                if out_path.exists() and out_path.stat().st_size > 0:
                    logger.info(f"  Skipping (already decompressed): {out_path.name}")
                    continue
                
                if dry_run:
                    logger.info(f"[DRY RUN] Would decompress: {gz_path.name} -> {out_path.name}")
                    continue
                
                try:
                    import gzip
                    with gzip.open(gz_path, 'rb') as f_in, open(out_path, 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)
                    logger.info(f"✓ Decompressed: {out_path.name}")
                except Exception as e:
                    logger.error(f"Failed to decompress {gz_path}: {e}")
    
    return linked_refs


def prepare_orf_classify_ensembl_dir(workdir, args, dry_run=False):
    """Prepare Ensembl annotation files for ORF classification"""
    if args.orf_classify_ensembl_dir:
        ensembl_dir = Path(args.orf_classify_ensembl_dir).resolve()
        if not ensembl_dir.exists():
            logger.error(f"Ensembl directory not found: {ensembl_dir}")
            return None
        return ensembl_dir

    if not getattr(args, 'prepare_orf_classify_ensembl', True):
        return None

    if args.species not in ENSEMBL_SPECIES_MAP:
        logger.warning(f"No Ensembl species mapping for: {args.species}")
        return None

    ensembl_species = ENSEMBL_SPECIES_MAP[args.species]
    division = ENSEMBL_DIVISION_MAP.get(args.species, 'ensembl')
    release = args.ensembl_release or DEFAULT_ENSEMBL_RELEASE.get(division)
    assembly = args.ensembl_assembly or DEFAULT_ENSEMBL_ASSEMBLY.get(args.species) or args.genome

    if not release:
        logger.error("Ensembl release not specified and no default available.")
        return None
    if not assembly:
        logger.error("Ensembl assembly not specified.")
        return None

    if args.ensembl_base_url:
        base_url = args.ensembl_base_url
    elif division == 'plants':
        base_url = f"https://ftp.ensemblgenomes.org/pub/plants/release-{release}"
    else:
        base_url = f"https://ftp.ensembl.org/pub/release-{release}"

    ensembl_dir = Path(workdir) / 'reference' / 'ensembl' / f"Ens{release}_{ensembl_species}"
    prep_script = Path(__file__).parent / 'retrieve_ensembl_data.sh'
    if not prep_script.exists():
        logger.error(f"retrieve_ensembl_data.sh not found at: {prep_script}")
        return None

    cmd = [
        'bash', str(prep_script),
        '--species', ensembl_species,
        '--release', str(release),
        '--assembly', str(assembly),
        '--outdir', str(ensembl_dir),
        '--base-url', base_url,
    ]
    if args.ensembl_version_suffix:
        cmd.extend(['--version-suffix', args.ensembl_version_suffix])

    logger.info("\n" + "=" * 60)
    logger.info("Preparing Ensembl annotation for ORF classification")
    logger.info("=" * 60)
    logger.info(f"Ensembl species: {ensembl_species}")
    logger.info(f"Release: {release}")
    logger.info(f"Assembly: {assembly}")
    logger.info(f"Base URL: {base_url}")
    logger.info(f"Output: {ensembl_dir}")

    if dry_run:
        logger.info(f"[DRY RUN] Would execute: {' '.join(cmd)}")
        return ensembl_dir

    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        logger.error(f"Failed to prepare Ensembl annotation: {e}")
        return None

    return ensembl_dir


def setup_container_directory(workdir, container_dir, orfquant_container=None,
                              rpbp_container=None, unify_orf_container=None,
                              gencode_orf_mapper_container=None, ribowaltz_container=None,
                              dry_run=False):
    """Setup container directory by copying specified images into workdir"""
    logger.info("\n" + "=" * 60)
    logger.info("Setting up container directory")
    logger.info("=" * 60)
    
    target_dir = Path(workdir) / 'containers'
    linked_containers = {}
    
    # Copy specified containers (overwrite symlinks to avoid NFS squashfuse timeout)
    if orfquant_container:
        orfquant_path = Path(orfquant_container).resolve()
        if orfquant_path.exists():
            link_path = target_dir / 'orfquant_patched.sif'
            if dry_run:
                logger.info(f"[DRY RUN] Would copy ORFquant: {orfquant_path} -> {link_path}")
            else:
                if link_path.is_symlink():
                    link_path.unlink()
                if not link_path.exists():
                    shutil.copy2(orfquant_path, link_path)
                    logger.info(f"✓ Copied ORFquant container: {orfquant_path.name}")
                else:
                    logger.info(f"  Skipping (exists): {link_path.name}")
            linked_containers['orfquant'] = link_path

    if rpbp_container:
        rpbp_path = Path(rpbp_container).resolve()
        if rpbp_path.exists():
            link_path = target_dir / 'rpbp.sif'
            if dry_run:
                logger.info(f"[DRY RUN] Would copy RPBP: {rpbp_path} -> {link_path}")
            else:
                if link_path.is_symlink(): link_path.unlink()
                if not link_path.exists():
                    shutil.copy2(rpbp_path, link_path)
                    logger.info(f"✓ Copied RPBP container: {rpbp_path.name}")
                else:
                    logger.info(f"  Skipping (exists): {link_path.name}")
            linked_containers['rpbp'] = link_path

    if unify_orf_container:
        unify_path = Path(unify_orf_container).resolve()
        if unify_path.exists():
            link_path = target_dir / 'unify_orf.sif'
            if dry_run:
                logger.info(f"[DRY RUN] Would copy Unify ORF: {unify_path} -> {link_path}")
            else:
                if link_path.is_symlink(): link_path.unlink()
                if not link_path.exists():
                    shutil.copy2(unify_path, link_path)
                    logger.info(f"✓ Copied Unify ORF container: {unify_path.name}")
                else:
                    logger.info(f"  Skipping (exists): {link_path.name}")
            linked_containers['unify_orf'] = link_path

    if gencode_orf_mapper_container:
        gencode_path = Path(gencode_orf_mapper_container).resolve()
        if gencode_path.exists():
            link_path = target_dir / 'gencode_orf_mapper.sif'
            if dry_run:
                logger.info(f"[DRY RUN] Would copy gencode-orf-mapper: {gencode_path} -> {link_path}")
            else:
                if link_path.is_symlink(): link_path.unlink()
                if not link_path.exists():
                    shutil.copy2(gencode_path, link_path)
                    logger.info(f"✓ Copied gencode-orf-mapper container: {gencode_path.name}")
                else:
                    logger.info(f"  Skipping (exists): {link_path.name}")
            linked_containers['gencode_orf_mapper'] = link_path

    if ribowaltz_container:
        rw_path = Path(ribowaltz_container).resolve()
        if rw_path.exists():
            link_path = target_dir / 'ribowaltz.sif'
            if dry_run:
                logger.info(f"[DRY RUN] Would copy riboWaltz: {rw_path} -> {link_path}")
            else:
                if not link_path.exists():
                    shutil.copy2(rw_path, link_path)
                    logger.info(f"✓ Copied riboWaltz container: {rw_path.name}")
            linked_containers['ribowaltz'] = link_path

    # Backward-compatible: copy all containers from container_dir
    if container_dir:
        logger.warning("container-dir is deprecated. Prefer specifying containers directly.")
        container_patterns = ['*.sif', '*.img']
        files = []
        for pattern in container_patterns:
            for file_path in Path(container_dir).resolve().glob(pattern):
                if not file_path.is_file():
                    continue
                dest_path = target_dir / file_path.name
                if dry_run:
                    logger.info(f"[DRY RUN] Would copy container: {file_path} -> {dest_path}")
                    files.append(dest_path)
                    continue
                if not dest_path.exists():
                    shutil.copy2(file_path, dest_path)
                files.append(dest_path)
        logger.info(f"\nTotal containers copied from directory: {len(files)}")
        for f in files:
            linked_containers[f.stem] = f
    
    return linked_containers


def generate_sample_sheet(workdir, data_dir, strandedness='auto', 
                          sample_type='riboseq', group_map=None, dry_run=False):
    """Generate sample sheet CSV using get_sample_sheet.py"""
    logger.info("\n" + "=" * 60)
    logger.info("Generating sample sheet")
    logger.info("=" * 60)
    
    # Find get_sample_sheet.py script
    script_dir = Path(__file__).parent
    sample_sheet_script = script_dir / 'get_sample_sheet.py'
    
    if not sample_sheet_script.exists():
        logger.error(f"get_sample_sheet.py not found at: {sample_sheet_script}")
        return None
    
    # Use data directory in workdir
    data_path = Path(workdir) / 'data'
    output_csv = Path(workdir) / 'scripts' / 'samplesheet.csv'
    
    cmd = [
        'python3',
        str(sample_sheet_script),
        '-i', str(data_path),
        '-o', str(output_csv),
        '--strandedness', strandedness,
        '--type', sample_type
    ]

    if group_map:
        group_map_path = Path(group_map).resolve()
        if group_map_path.exists():
            cmd.extend(['--group-map', str(group_map_path)])
            logger.info(f"Using group map: {group_map_path}")
        else:
            logger.warning(f"Group map file not found, ignoring: {group_map_path}")
    
    if dry_run:
        logger.info(f"[DRY RUN] Would execute: {' '.join(cmd)}")
        return output_csv
    
    try:
        logger.info(f"Executing: {' '.join(cmd)}")
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        logger.info(result.stdout)
        logger.info(f"✓ Sample sheet generated: {output_csv}")
        return output_csv
    except subprocess.CalledProcessError as e:
        logger.error(f"Failed to generate sample sheet: {e}")
        logger.error(e.stderr)
        return None


def generate_nextflow_script(workdir, args, sample_sheet, containers,
                             references, orf_classify_ensembl_dir=None, dry_run=False):
    """Generate Nextflow execution shell script"""
    logger.info("\n" + "=" * 60)
    logger.info("Generating Nextflow execution script")
    logger.info("=" * 60)
    
    workdir = Path(workdir).resolve()
    # Place script in process directory for better organization
    script_path = workdir / 'process' / args.script_name
    
    # Determine the pipeline directory (parent of scripts directory)
    pipeline_dir = Path(__file__).resolve().parent.parent
    main_nf = pipeline_dir / 'main.nf'
    
    # Build Nextflow command (based on validated reference parameters)
    nf_cmd_parts = [
        f"nextflow run {main_nf}",
        f"-profile {args.profile}",
        f"-w {workdir / 'process' / 'work'}",
        f"--input {sample_sheet}",
        f"--outdir {workdir / 'result'}",
        f"--aligner {args.aligner}",
    ]
    
    # Add resource limits
    nf_cmd_parts.extend([
        f"--max_memory '{args.max_memory}'",
        f"--max_cpus {args.max_cpus}",
        f"--max_time '{args.max_time}'",
    ])
    
    # Background mode and save reference
    if getattr(args, 'bg', False):
        nf_cmd_parts.append("-bg")
    if args.save_reference:
        nf_cmd_parts.append("--save_reference true")
    
    # Add container paths
    if 'orfquant' in containers:
        nf_cmd_parts.append(f"--orfquant_container {containers['orfquant']}")
    
    if 'rpbp' in containers:
        nf_cmd_parts.append(f"--rpbp_container {containers['rpbp']}")
    
    if 'unify_orf' in containers:
        nf_cmd_parts.append(f"--unify_orf_container {containers['unify_orf']}")
    
    if 'gencode_orf_mapper' in containers:
        nf_cmd_parts.append(f"--gencode_orf_mapper_container {containers['gencode_orf_mapper']}")

    if 'ribowaltz' in containers:
        nf_cmd_parts.append(f"--ribowaltz_container {containers['ribowaltz']}")

    if 'deseq2' in containers:
        nf_cmd_parts.append(f"--deseq2_container {containers['deseq2']}")

    if 'price' in containers:
        nf_cmd_parts.append(f"--price_container {containers['price']}")

    if 'orf_qc' in containers:
        nf_cmd_parts.append(f"--orf_qc_container {containers['orf_qc']}")

    if 'expression_quant' in containers:
        nf_cmd_parts.append(f"--expression_quant_container {containers['expression_quant']}")

    def prefer_uncompressed(path_obj):
        path_obj = Path(path_obj)
        if path_obj.suffix == '.gz':
            uncompressed = path_obj.with_suffix('')
            if uncompressed.exists():
                return uncompressed
        return path_obj

    # Add reference files if provided
    # Use species-specific file naming: {species}.genome.fa, {species}.gtf, {species}.transcripts.fa
    contaminant_provided = False
    if references:
        if 'genome_fasta' in references and references['genome_fasta']:
            fasta_path = prefer_uncompressed(references['genome_fasta'][0])
            nf_cmd_parts.append(f"--fasta {fasta_path}")
        if 'gtf' in references and references['gtf']:
            gtf_path = prefer_uncompressed(references['gtf'][0])
            nf_cmd_parts.append(f"--gtf {gtf_path}")
        if 'transcripts' in references and references['transcripts']:
            transcript_path = prefer_uncompressed(references['transcripts'][0])
            nf_cmd_parts.append(f"--transcript_fasta {transcript_path}")
        if 'contaminant' in references and references['contaminant']:
            # Prefer *_final_contamination.fasta from prepare_reference_db_v2.2.py
            contam_files = references['contaminant']
            # Sort to prioritize *_final_contamination.fasta
            final_contam = [f for f in contam_files if 'final_contamination' in str(f)]
            if final_contam:
                contam_path = prefer_uncompressed(final_contam[0])
            else:
                contam_path = prefer_uncompressed(contam_files[0])
            nf_cmd_parts.append(f"--contaminant_fasta {contam_path}")
            contaminant_provided = True
    
    # Skip contaminant filter if no contaminant FASTA provided
    if not contaminant_provided:
        nf_cmd_parts.append("--skip_contaminant_filter")

    if orf_classify_ensembl_dir:
        nf_cmd_parts.append(f"--orf_classify_ensembl_dir {orf_classify_ensembl_dir}")
    
    # Add pipeline options
    if args.run_prefilter_qc:
        nf_cmd_parts.append("--run_prefilter_qc true")
    
    # Skip RPBP by default (not fully validated yet)
    nf_cmd_parts.append("--skip_rpbp true")
    
    # Replicate BAM merging
    if getattr(args, 'merge_replicates', False):
        nf_cmd_parts.append("--merge_replicates true")

    # Ribotricer phase score cutoff (only set when explicitly provided)
    if getattr(args, 'ribotricer_phase_score_cutoff', None) is not None:
        nf_cmd_parts.append(f"--ribotricer_phase_score_cutoff {args.ribotricer_phase_score_cutoff}")

    # sORF unique mapping mode (only set when explicitly provided)
    if getattr(args, 'sorf_unique_mode', None) is not None:
        nf_cmd_parts.append(f"--sorf_unique_mode {args.sorf_unique_mode}")

    # Skip collect QC stats
    if getattr(args, 'skip_collect_qc_stats', False):
        nf_cmd_parts.append("--skip_collect_qc_stats true")

    # P-site offset correction for ORFquant (enabled by default)
    if hasattr(args, 'orfquant_psite_correction'):
        if args.orfquant_psite_correction:
            nf_cmd_parts.append("--orfquant_psite_correction true")
        else:
            nf_cmd_parts.append("--orfquant_psite_correction false")
    
    # ORF unification and classification options (default: run, so no need to set false)
    # Note: orf_classify_mode is deprecated; all three classifiers now run simultaneously.
    nf_cmd_parts.extend([
        f"--unify_orf_min_len {args.unify_orf_min_len}",
    ])
    
    # Advanced ORF merging options
    if hasattr(args, 'unify_orf_frame_merge_min_overlap'):
        nf_cmd_parts.append(f"--unify_orf_frame_merge_min_overlap {args.unify_orf_frame_merge_min_overlap}")
    if hasattr(args, 'unify_orf_no_frame_merge') and args.unify_orf_no_frame_merge:
        nf_cmd_parts.append("--unify_orf_no_frame_merge true")
    if hasattr(args, 'unify_orf_seq_cluster') and args.unify_orf_seq_cluster:
        nf_cmd_parts.append("--unify_orf_seq_cluster true")
    
    # Generate script content
    script_content = f"""#!/bin/bash
# =============================================================================
# Ribo-seq Pipeline Execution Script
# Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}
# Working Directory: {workdir}
# =============================================================================

set -euo pipefail

# --- Environment Setup ---
WORKDIR="{workdir}"
PROCESS_DIR="${{WORKDIR}}/process"
RESULT_DIR="${{WORKDIR}}/result"

# Create process directory if not exists
mkdir -p "${{PROCESS_DIR}}"
mkdir -p "${{RESULT_DIR}}"

# Change to process directory (where this script is located)
cd "${{PROCESS_DIR}}"

# --- Pipeline Execution ---
echo "=========================================="
echo "Starting Ribo-seq Pipeline"
echo "=========================================="
echo "Working directory: ${{WORKDIR}}"
echo "Process directory: ${{PROCESS_DIR}}"
echo "Sample sheet: {sample_sheet}"
echo "Output directory: ${{RESULT_DIR}}"
echo "Work directory: ${{PROCESS_DIR}}/work"
echo "=========================================="
echo ""

# Nextflow command
# WARNING: Do NOT comment out lines with # inside this multi-line command.
# Due to bash line-continuation (backslash), a # in the middle of the command
# will comment out ALL subsequent arguments. To disable an argument, DELETE the line.
{' \\\n    '.join(nf_cmd_parts)} \\
    -resume \\
    -with-report "${{RESULT_DIR}}/pipeline_report.html" \\
    -with-timeline "${{RESULT_DIR}}/timeline.html" \\
    -with-dag "${{RESULT_DIR}}/flowchart.html"

# --- Completion ---
echo ""
echo "=========================================="
echo "Pipeline execution completed!"
echo "Results: ${{RESULT_DIR}}"
echo "Reports: ${{RESULT_DIR}}/pipeline_report.html"
echo "=========================================="
"""
    
    if dry_run:
        logger.info(f"[DRY RUN] Would create script: {script_path}")
        logger.info("\n" + "=" * 60)
        logger.info("Script content preview:")
        logger.info("=" * 60)
        print(script_content)
    else:
        # Ensure process directory exists
        script_path.parent.mkdir(parents=True, exist_ok=True)
        
        with open(script_path, 'w') as f:
            f.write(script_content)
        
        # Make executable
        script_path.chmod(0o755)
        logger.info(f"✓ Execution script created: {script_path}")
        logger.info(f"  Run with: bash {script_path}")
        logger.info(f"  Or from workdir: bash process/{args.script_name}")
    
    return script_path


def create_config_summary(workdir, args, sample_sheet, containers,
                          references, script_path, auto_reference_info=None,
                          orf_classify_ensembl_dir=None, dry_run=False):
    """Create a summary configuration file"""
    workdir = Path(workdir).resolve()
    summary_path = workdir / 'scripts' / 'workflow_config.json'
    
    config = {
        'workflow_setup': {
            'created_at': datetime.now().isoformat(),
            'working_directory': str(workdir),
            'dry_run': dry_run
        },
        'data': {
            'source_directory': args.data_dir,
            'sample_sheet': str(sample_sheet) if sample_sheet else None,
            'strandedness': args.strandedness,
            'sample_type': args.sample_type
        },
        'reference': {
            'genome': args.genome,
            'species': args.species,
            'source_directory': args.reference_dir or (str(auto_reference_info.get('reference_dir')) if auto_reference_info else None),
            'auto_generated': bool(auto_reference_info),
            'auto_reference_base_dir': str(auto_reference_info.get('base_dir')) if auto_reference_info else None,
            'auto_reference_dirs': [str(auto_reference_info.get('reference_dir')), str(auto_reference_info.get('contaminant_dir'))] if auto_reference_info else None,
            'files': {k: [str(f) for f in v] for k, v in references.items()} if references else {}
        },
        'orf_classify_ensembl': {
            'ensembl_dir': str(orf_classify_ensembl_dir) if orf_classify_ensembl_dir else None,
            'ensembl_release': args.ensembl_release,
            'ensembl_assembly': args.ensembl_assembly or DEFAULT_ENSEMBL_ASSEMBLY.get(args.species) or args.genome,
            'ensembl_base_url': args.ensembl_base_url,
            'auto_prepared': bool(orf_classify_ensembl_dir),
        },
        'containers': {
            'source_directory': args.container_dir,
            'orfquant': str(containers.get('orfquant', '')),
            'rpbp': str(containers.get('rpbp', '')),
            'unify_orf': str(containers.get('unify_orf', '')),
            'gencode_orf_mapper': str(containers.get('gencode_orf_mapper', ''))
        },
        'pipeline_options': {
            'run_prefilter_qc': args.run_prefilter_qc,
            'unify_orf_min_len': args.unify_orf_min_len,
            'unify_orf_frame_merge_min_overlap': getattr(args, 'unify_orf_frame_merge_min_overlap', 0.9),
            'unify_orf_no_frame_merge': getattr(args, 'unify_orf_no_frame_merge', False),
            'unify_orf_seq_cluster': getattr(args, 'unify_orf_seq_cluster', False),
            'profile': args.profile
        },
        'scripts': {
            'execution_script': str(script_path) if script_path else None
        }
    }
    
    if dry_run:
        logger.info(f"\n[DRY RUN] Would create config: {summary_path}")
    else:
        with open(summary_path, 'w') as f:
            json.dump(config, f, indent=2)
        logger.info(f"\n✓ Configuration summary saved: {summary_path}")
    
    return summary_path


def print_final_summary(workdir, sample_sheet, script_path):
    """Print final summary"""
    logger.info("\n" + "=" * 60)
    logger.info("WORKFLOW PREPARATION COMPLETE")
    logger.info("=" * 60)
    logger.info(f"\nWorking directory: {workdir}")
    logger.info(f"\nDirectory structure:")
    for dir_name in STANDARD_DIRS.keys():
        logger.info(f"  - {dir_name}/")
    
    if sample_sheet:
        logger.info(f"\nSample sheet: {sample_sheet}")
    
    if script_path:
        logger.info(f"\nExecution script: {script_path}")
        logger.info(f"\nTo start the pipeline, run:")
        logger.info(f"  bash {script_path}")
    
    logger.info("\n" + "=" * 60)


def main():
    """Main execution function"""
    args = parse_args()
    
    # Validate inputs
    if not os.path.exists(args.data_dir):
        logger.error(f"Data directory does not exist: {args.data_dir}")
        sys.exit(1)
    
    workdir = Path(args.workdir).resolve()
    
    logger.info("=" * 60)
    logger.info("Ribo-seq Workflow Preparation")
    logger.info("=" * 60)
    logger.info(f"Working directory: {workdir}")
    logger.info(f"Data directory: {args.data_dir}")
    if args.reference_dir:
        logger.info(f"Reference directory: {args.reference_dir}")
    if args.contaminant_dir:
        logger.info(f"Contaminant directory: {args.contaminant_dir}")
    if args.container_dir:
        logger.info(f"Container directory: {args.container_dir}")
    logger.info(f"Genome: {args.genome}")
    logger.info(f"Species: {args.species}")
    if args.dry_run:
        logger.info("\n*** DRY RUN MODE - No changes will be made ***")
    logger.info("")
    
    # Step 1: Create directory structure
    create_directory_structure(workdir, args.dry_run)
    
    # Step 2: Setup data directory (with SRA conversion support)
    setup_data_directory(workdir, args.data_dir, args.dry_run,
                         sra_threads=args.sra_threads, pigz_threads=args.pigz_threads)
    
    # Step 3: Setup reference directory
    reference_dirs = []
    auto_reference_info = None
    
    # Collect user-specified directories
    if args.reference_dir:
        reference_dirs.append(args.reference_dir)
        logger.info(f"Using reference directory: {args.reference_dir}")
    if args.contaminant_dir:
        reference_dirs.append(args.contaminant_dir)
        logger.info(f"Using contaminant directory: {args.contaminant_dir}")
    
    # If no directories specified, use auto-generation
    if not reference_dirs:
        logger.info("No reference directory provided. Preparing reference database...")
        auto_reference_info = prepare_reference_db_if_missing(workdir, args.species, args.dry_run)
        if auto_reference_info:
            logger.info(f"Auto reference base dir: {auto_reference_info['base_dir']}")
            reference_dirs = [auto_reference_info['reference_dir']]
            if auto_reference_info.get('contaminant_dir'):
                reference_dirs.append(auto_reference_info['contaminant_dir'])
        else:
            logger.error("Failed to prepare reference database.")
            sys.exit(1)

    references = setup_reference_directory(workdir, reference_dirs, args.species, args.dry_run)
    
    # Validate required references exist, regenerate if missing
    required_refs = ['genome_fasta', 'gtf', 'contaminant']
    missing_refs = [r for r in required_refs if r not in references or not references[r]]
    
    if missing_refs and not auto_reference_info:
        logger.warning(f"Missing required reference files: {missing_refs}")
        logger.info("Attempting to prepare missing references using prepare_reference_db_v2.2.py...")
        auto_reference_info = prepare_reference_db_if_missing(workdir, args.species, args.dry_run)
        if auto_reference_info:
            # Add auto-generated directories to reference search
            if auto_reference_info.get('reference_dir'):
                reference_dirs.append(auto_reference_info['reference_dir'])
            if auto_reference_info.get('contaminant_dir'):
                reference_dirs.append(auto_reference_info['contaminant_dir'])
            # Re-run reference setup
            references = setup_reference_directory(workdir, reference_dirs, args.species, args.dry_run)
    
    # Final validation
    final_missing = [r for r in required_refs if r not in references or not references[r]]
    if final_missing:
        logger.warning(f"Could not find all required references. Missing: {final_missing}")
        logger.warning("Pipeline may fail. Consider providing explicit reference paths.")

    # Step 4: Prepare Ensembl annotation for ORF classification
    orf_classify_ensembl_dir = prepare_orf_classify_ensembl_dir(
        workdir,
        args,
        args.dry_run
    )

    # Step 5: Setup container directory
    containers = setup_container_directory(
        workdir,
        args.container_dir,
        args.orfquant_container,
        args.rpbp_container,
        args.unify_orf_container,
        args.gencode_orf_mapper_container,
        args.ribowaltz_container,
        args.dry_run
    )
    
    # Step 6: Generate sample sheet
    sample_sheet = generate_sample_sheet(
        workdir,
        args.data_dir,
        args.strandedness,
        args.sample_type,
        getattr(args, 'group_map', None),
        args.dry_run
    )
    
    # Step 7: Generate execution script
    script_path = generate_nextflow_script(
        workdir,
        args,
        sample_sheet,
        containers,
        references,
        orf_classify_ensembl_dir,
        args.dry_run
    )
    
    # Step 8: Create configuration summary
    create_config_summary(
        workdir,
        args,
        sample_sheet,
        containers,
        references,
        script_path,
        auto_reference_info,
        orf_classify_ensembl_dir,
        args.dry_run
    )
    
    # Final summary
    if not args.dry_run:
        print_final_summary(workdir, sample_sheet, script_path)


if __name__ == "__main__":
    main()
