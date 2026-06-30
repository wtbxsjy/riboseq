#!/usr/bin/env python3
"""
Validate gencode converter output files

Checks:
1. BED file has unique ORF names in column 4
2. FASTA sequences match their declared length in the header
3. All sequences end with stop codon
4. BED and FASTA have matching ORF names
"""

import argparse
import sys
import re

def validate_bed_file(bed_file):
    """Validate BED file format"""
    issues = []
    orf_names = []
    line_count = 0
    
    with open(bed_file, 'r') as f:
        for line_num, line in enumerate(f, 1):
            if line.startswith('#'):
                continue
                
            line_count += 1
            parts = line.strip().split('\t')
            
            if len(parts) != 6:
                issues.append(f"Line {line_num}: Expected 6 columns, got {len(parts)}")
                continue
            
            chrom, start, end, orf_name, study_id, strand = parts
            
            # Check ORF name format
            if not re.match(r'^[A-Z0-9]+_\d+_\d+aa$', orf_name):
                issues.append(f"Line {line_num}: Invalid ORF name format '{orf_name}' (expected GENE_START_LENaa)")
            
            # Check for duplicate ORF names
            if orf_name in orf_names:
                issues.append(f"Line {line_num}: Duplicate ORF name '{orf_name}'")
            orf_names.append(orf_name)
            
            # Check coordinates
            try:
                start_int = int(start)
                end_int = int(end)
                if start_int >= end_int:
                    issues.append(f"Line {line_num}: Invalid coordinates start={start} >= end={end}")
            except ValueError:
                issues.append(f"Line {line_num}: Non-numeric coordinates start={start}, end={end}")
            
            # Check strand
            if strand not in ['+', '-']:
                issues.append(f"Line {line_num}: Invalid strand '{strand}' (expected + or -)")
    
    return {
        'file': bed_file,
        'total_lines': line_count,
        'unique_orfs': len(set(orf_names)),
        'issues': issues,
        'orf_names': orf_names
    }

def validate_fasta_file(fasta_file):
    """Validate FASTA file format"""
    issues = []
    orf_data = {}
    current_header = None
    current_seq = []
    line_count = 0
    
    with open(fasta_file, 'r') as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line:
                continue
                
            if line.startswith('>'):
                # Process previous entry
                if current_header:
                    seq = ''.join(current_seq)
                    orf_data[current_header] = seq
                    
                    # Validate sequence
                    if not seq.endswith('*'):
                        issues.append(f"ORF '{current_header}': Sequence doesn't end with stop codon")
                    
                    # Extract expected length from header
                    match = re.search(r'_(\d+)aa--', current_header)
                    if match:
                        expected_length = int(match.group(1))
                        actual_length = len(seq.rstrip('*'))
                        if actual_length != expected_length:
                            issues.append(f"ORF '{current_header}': Length mismatch (expected {expected_length}aa, got {actual_length}aa)")
                
                # Start new entry
                current_header = line[1:]  # Remove '>'
                current_seq = []
                line_count += 1
                
                # Validate header format
                if not re.match(r'^[A-Z0-9]+_\d+_\d+aa--\S+$', current_header):
                    issues.append(f"Line {line_num}: Invalid header format '{current_header}'")
            else:
                current_seq.append(line)
        
        # Process last entry
        if current_header:
            seq = ''.join(current_seq)
            orf_data[current_header] = seq
            
            if not seq.endswith('*'):
                issues.append(f"ORF '{current_header}': Sequence doesn't end with stop codon")
            
            match = re.search(r'_(\d+)aa--', current_header)
            if match:
                expected_length = int(match.group(1))
                actual_length = len(seq.rstrip('*'))
                if actual_length != expected_length:
                    issues.append(f"ORF '{current_header}': Length mismatch (expected {expected_length}aa, got {actual_length}aa)")
    
    return {
        'file': fasta_file,
        'total_orfs': line_count,
        'issues': issues,
        'orf_data': orf_data
    }

def cross_validate(bed_result, fasta_result):
    """Cross-validate BED and FASTA files"""
    issues = []
    
    # Extract ORF names from FASTA headers (before --)
    fasta_orfs = set()
    for header in fasta_result['orf_data'].keys():
        orf_name = header.split('--')[0]
        fasta_orfs.add(orf_name)
    
    bed_orfs = set(bed_result['orf_names'])
    
    # Check if sets match
    only_in_bed = bed_orfs - fasta_orfs
    only_in_fasta = fasta_orfs - bed_orfs
    
    if only_in_bed:
        issues.append(f"ORFs only in BED file: {len(only_in_bed)} ({list(only_in_bed)[:5]}...)")
    
    if only_in_fasta:
        issues.append(f"ORFs only in FASTA file: {len(only_in_fasta)} ({list(only_in_fasta)[:5]}...)")
    
    return issues

def main():
    parser = argparse.ArgumentParser(description='Validate gencode converter output')
    parser.add_argument('--bed', required=True, help='BED file to validate')
    parser.add_argument('--fasta', required=True, help='FASTA file to validate')
    
    args = parser.parse_args()
    
    print("=" * 70)
    print("GENCODE Output Validation")
    print("=" * 70)
    
    # Validate BED
    print(f"\n📄 Validating BED file: {args.bed}")
    bed_result = validate_bed_file(args.bed)
    print(f"   Total entries: {bed_result['total_lines']}")
    print(f"   Unique ORFs: {bed_result['unique_orfs']}")
    
    if bed_result['issues']:
        print(f"   ⚠️  Issues found: {len(bed_result['issues'])}")
        for issue in bed_result['issues'][:10]:  # Show first 10
            print(f"      - {issue}")
        if len(bed_result['issues']) > 10:
            print(f"      ... and {len(bed_result['issues']) - 10} more")
    else:
        print("   ✅ No issues found")
    
    # Validate FASTA
    print(f"\n📄 Validating FASTA file: {args.fasta}")
    fasta_result = validate_fasta_file(args.fasta)
    print(f"   Total ORFs: {fasta_result['total_orfs']}")
    
    if fasta_result['issues']:
        print(f"   ⚠️  Issues found: {len(fasta_result['issues'])}")
        for issue in fasta_result['issues'][:10]:  # Show first 10
            print(f"      - {issue}")
        if len(fasta_result['issues']) > 10:
            print(f"      ... and {len(fasta_result['issues']) - 10} more")
    else:
        print("   ✅ No issues found")
    
    # Cross-validate
    print(f"\n🔗 Cross-validating BED and FASTA")
    cross_issues = cross_validate(bed_result, fasta_result)
    
    if cross_issues:
        print(f"   ⚠️  Issues found: {len(cross_issues)}")
        for issue in cross_issues:
            print(f"      - {issue}")
    else:
        print("   ✅ BED and FASTA files are consistent")
    
    # Summary
    print("\n" + "=" * 70)
    total_issues = len(bed_result['issues']) + len(fasta_result['issues']) + len(cross_issues)
    
    if total_issues == 0:
        print("✅ VALIDATION PASSED - No issues found!")
        return 0
    else:
        print(f"⚠️  VALIDATION FAILED - {total_issues} issues found")
        return 1

if __name__ == '__main__':
    sys.exit(main())
