#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Test Script for prepare_workflow.py
===================================

This test validates that the generated run_pipeline.sh script contains
correct Nextflow parameters based on the expected format:

Expected parameters:
- --fasta: Points to {species}.genome.fa (NOT transcripts)
- --gtf: Points to {species}.gtf
- --transcript_fasta: Points to {species}.transcripts.fa
- --contaminant_fasta: Points to {species}_final_contamination.fasta
- --aligner: star (default)
- --max_memory, --max_cpus, --max_time: Resource limits
- -bg: Background mode
- --save_reference: Save references
- --skip_rpbp: Skip RPBP by default
- NO --genome parameter (should not be present)

Usage:
    python3 test_prepare_workflow.py --test-dir /path/to/test/workdir
"""

import os
import sys
import argparse
import tempfile
import shutil
import subprocess
import re
from pathlib import Path


def create_mock_environment(test_dir):
    """Create mock data and reference files for testing"""
    test_dir = Path(test_dir)
    
    # Create directories
    data_dir = test_dir / 'mock_data'
    ref_dir = test_dir / 'mock_reference'
    contam_dir = test_dir / 'mock_contamination'
    
    for d in [data_dir, ref_dir, contam_dir]:
        d.mkdir(parents=True, exist_ok=True)
    
    # Create mock FASTQ files
    (data_dir / 'sample1_R1.fastq.gz').touch()
    (data_dir / 'sample2_R1.fastq.gz').touch()
    
    # Create mock reference files (simulating prepare_reference_db_v2.2.py output)
    species = 'rice'
    (ref_dir / f'{species}.genome.fa').write_text('>chr1\nACGT\n')
    (ref_dir / f'{species}.gtf').write_text('# GTF mock\n')
    (ref_dir / f'{species}.transcripts.fa').write_text('>transcript1\nACGT\n')
    
    # Create mock contamination file
    (contam_dir / f'{species}_final_contamination.fasta').write_text('>contam1\nACGT\n')
    
    return {
        'data_dir': data_dir,
        'ref_dir': ref_dir,
        'contam_dir': contam_dir,
        'species': species
    }


def run_prepare_workflow(workdir, mock_env, dry_run=True):
    """Run prepare_workflow.py with mock environment"""
    script_dir = Path(__file__).parent
    prepare_script = script_dir / 'prepare_workflow.py'
    
    if not prepare_script.exists():
        raise FileNotFoundError(f"prepare_workflow.py not found at {prepare_script}")
    
    cmd = [
        'python3', str(prepare_script),
        '-w', str(workdir),
        '-d', str(mock_env['data_dir']),
        '-r', str(mock_env['ref_dir']),
        '--species', mock_env['species'],
        '--profile', 'test_local_singularity',
        '--skip-orf-classify-ensembl',
    ]
    
    # Add contamination directory as additional reference
    # The script should find files in both ref_dir and contam_dir
    # For this test, we simulate by creating symlinks or passing both
    
    if dry_run:
        cmd.append('--dry-run')
    
    print(f"Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    if result.returncode != 0:
        print(f"STDOUT: {result.stdout}")
        print(f"STDERR: {result.stderr}")
        raise RuntimeError(f"prepare_workflow.py failed with exit code {result.returncode}")
    
    return result


def validate_generated_script(script_path, mock_env):
    """Validate the generated run_pipeline.sh script"""
    if not script_path.exists():
        return False, "Script file does not exist"
    
    content = script_path.read_text()
    errors = []
    warnings = []
    passed = []
    
    species = mock_env['species']
    
    # === CRITICAL CHECKS ===
    
    # 1. --fasta should point to genome.fa, NOT transcripts.fa
    fasta_match = re.search(r'--fasta\s+(\S+)', content)
    if fasta_match:
        fasta_path = fasta_match.group(1)
        if 'genome' in fasta_path:
            passed.append(f"✓ --fasta correctly points to genome file: {fasta_path}")
        elif 'transcript' in fasta_path.lower():
            errors.append(f"✗ --fasta incorrectly points to transcripts: {fasta_path}")
        else:
            warnings.append(f"? --fasta path unclear: {fasta_path}")
    else:
        errors.append("✗ --fasta parameter not found")
    
    # 2. --genome should NOT be present (we use --fasta instead)
    if re.search(r'--genome\s+\w+', content):
        errors.append("✗ --genome parameter found (should not be present)")
    else:
        passed.append("✓ --genome parameter correctly absent")
    
    # 3. --gtf should be present
    if re.search(r'--gtf\s+\S+', content):
        passed.append("✓ --gtf parameter found")
    else:
        errors.append("✗ --gtf parameter not found")
    
    # 4. --transcript_fasta should be present
    transcript_match = re.search(r'--transcript_fasta\s+(\S+)', content)
    if transcript_match:
        transcript_path = transcript_match.group(1)
        if 'transcript' in transcript_path.lower() or 'cdna' in transcript_path.lower():
            passed.append(f"✓ --transcript_fasta correctly set: {transcript_path}")
        else:
            warnings.append(f"? --transcript_fasta path may be wrong: {transcript_path}")
    else:
        warnings.append("? --transcript_fasta parameter not found (optional)")
    
    # 5. --contaminant_fasta should point to *_final_contamination.fasta
    contam_match = re.search(r'--contaminant_fasta\s+(\S+)', content)
    if contam_match:
        contam_path = contam_match.group(1)
        if 'final_contamination' in contam_path:
            passed.append(f"✓ --contaminant_fasta correctly points to final_contamination: {contam_path}")
        else:
            warnings.append(f"? --contaminant_fasta may not be optimal: {contam_path}")
    else:
        # Check if skip_contaminant_filter is set instead
        if '--skip_contaminant_filter' in content:
            passed.append("✓ --skip_contaminant_filter set (no contaminant provided)")
        else:
            errors.append("✗ Neither --contaminant_fasta nor --skip_contaminant_filter found")
    
    # 6. --aligner should be present (default: star)
    if re.search(r'--aligner\s+(star|hisat2)', content):
        passed.append("✓ --aligner parameter found")
    else:
        errors.append("✗ --aligner parameter not found")
    
    # 7. Resource limits should be present
    if re.search(r"--max_memory\s+'[\d.]+\.GB'", content):
        passed.append("✓ --max_memory parameter found")
    else:
        warnings.append("? --max_memory parameter format may be wrong")
    
    if re.search(r'--max_cpus\s+\d+', content):
        passed.append("✓ --max_cpus parameter found")
    else:
        warnings.append("? --max_cpus parameter not found")
    
    if re.search(r"--max_time\s+'[\d.]+\.h'", content):
        passed.append("✓ --max_time parameter found")
    else:
        warnings.append("? --max_time parameter format may be wrong")
    
    # 8. -bg (background) should be present
    if re.search(r'\s-bg\s', content) or content.strip().endswith('-bg'):
        passed.append("✓ -bg (background) flag found")
    else:
        warnings.append("? -bg flag not found")
    
    # 9. --save_reference should be present
    if '--save_reference' in content:
        passed.append("✓ --save_reference flag found")
    else:
        warnings.append("? --save_reference flag not found")
    
    # 10. --skip_rpbp should be present (default behavior)
    if '--skip_rpbp' in content:
        passed.append("✓ --skip_rpbp flag found")
    else:
        warnings.append("? --skip_rpbp flag not found")
    
    # 11. -resume should be present
    if '-resume' in content:
        passed.append("✓ -resume flag found")
    else:
        warnings.append("? -resume flag not found")
    
    return errors, warnings, passed


def test_dry_run(test_dir, mock_env):
    """Test dry-run mode and validate output"""
    print("\n" + "=" * 60)
    print("TEST: Dry-run mode validation")
    print("=" * 60)
    
    result = run_prepare_workflow(test_dir, mock_env, dry_run=True)
    
    # In dry-run mode, check stdout for expected commands
    output = result.stdout + result.stderr
    
    checks = [
        ('--fasta', 'genome'),
        ('--gtf', ''),
        ('--aligner', 'star'),
    ]
    
    all_passed = True
    for param, expected_content in checks:
        if param in output:
            print(f"✓ {param} found in dry-run output")
            if expected_content and expected_content not in output:
                print(f"  Warning: Expected '{expected_content}' not found")
        else:
            print(f"? {param} not found in dry-run output (may be normal)")
    
    return all_passed


def test_full_generation(test_dir, mock_env):
    """Test full script generation and validate output file"""
    print("\n" + "=" * 60)
    print("TEST: Full script generation")
    print("=" * 60)
    
    workdir = Path(test_dir) / 'workdir'
    
    # Also copy contaminant files to reference dir for this test
    contam_file = mock_env['contam_dir'] / f"{mock_env['species']}_final_contamination.fasta"
    ref_contam = mock_env['ref_dir'] / f"{mock_env['species']}_final_contamination.fasta"
    if contam_file.exists() and not ref_contam.exists():
        shutil.copy(contam_file, ref_contam)
    
    result = run_prepare_workflow(workdir, mock_env, dry_run=False)
    
    # Check generated script
    script_path = workdir / 'process' / 'run_pipeline.sh'
    
    if not script_path.exists():
        print(f"✗ Script not generated at expected path: {script_path}")
        return False
    
    print(f"✓ Script generated: {script_path}")
    
    # Validate script content
    errors, warnings, passed = validate_generated_script(script_path, mock_env)
    
    print("\n--- Validation Results ---")
    for p in passed:
        print(p)
    for w in warnings:
        print(w)
    for e in errors:
        print(e)
    
    print(f"\nSummary: {len(passed)} passed, {len(warnings)} warnings, {len(errors)} errors")
    
    # Print actual script for debugging
    print("\n--- Generated Script Content ---")
    print(script_path.read_text())
    
    return len(errors) == 0


def main():
    parser = argparse.ArgumentParser(description="Test prepare_workflow.py parameter generation")
    parser.add_argument('--test-dir', default=None,
                        help='Directory for test files (default: temp directory)')
    parser.add_argument('--keep-temp', action='store_true',
                        help='Keep temporary test directory after completion')
    args = parser.parse_args()
    
    # Create test directory
    if args.test_dir:
        test_dir = Path(args.test_dir)
        test_dir.mkdir(parents=True, exist_ok=True)
        use_temp = False
    else:
        test_dir = Path(tempfile.mkdtemp(prefix='test_prepare_workflow_'))
        use_temp = True
    
    print(f"Test directory: {test_dir}")
    
    try:
        # Setup mock environment
        mock_env = create_mock_environment(test_dir)
        print(f"Mock environment created:")
        print(f"  Data: {mock_env['data_dir']}")
        print(f"  Reference: {mock_env['ref_dir']}")
        print(f"  Contamination: {mock_env['contam_dir']}")
        
        # Run tests
        all_passed = True
        
        # Test 1: Dry-run validation
        test_dry_run(test_dir, mock_env)
        
        # Test 2: Full generation
        if not test_full_generation(test_dir, mock_env):
            all_passed = False
        
        # Final result
        print("\n" + "=" * 60)
        if all_passed:
            print("ALL TESTS PASSED ✓")
        else:
            print("SOME TESTS FAILED ✗")
            sys.exit(1)
        print("=" * 60)
        
    finally:
        if use_temp and not args.keep_temp:
            print(f"\nCleaning up: {test_dir}")
            shutil.rmtree(test_dir)
        else:
            print(f"\nTest files kept at: {test_dir}")


if __name__ == "__main__":
    main()
