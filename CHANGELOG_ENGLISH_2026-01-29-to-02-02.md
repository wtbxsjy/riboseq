# Changelog - January 29 to February 2, 2026

## Overview

This document summarizes all changes to the RiboseQ Nextflow pipeline from January 29, 2026 through February 2, 2026. During this 4-day period, the project received significant enhancements including new features for unified ORF analysis, performance optimizations, bug fixes, and documentation improvements across approximately 60 commits.

---

## February 2, 2026

### Features

#### Input Validation and Error Handling
- **Commit**: f04f3d1
- Added input file existence checks for RiboseQC analysis
- Enhanced error handling to improve pipeline robustness and provide meaningful error messages to users
- Ensures validation of required input files before processing begins

#### Enhanced max_coverage Processing
- **Commit**: 2e5de46
- Improved max_coverage parameter handling to support both logical and numerical representations
- Added debug output for data structure inspection
- Enables better handling of complex coverage scenarios

---

## February 1, 2026

### Major Features

#### Unified ORF Prediction and Classification Integration
- **Commit**: 5f04842
- Implemented comprehensive unified ORF prediction functionality
- Integrated RiboseQC P-site statistics with ORF predictions
- Updated related documentation and scripts
- Provides consolidated ORF analysis pipeline combining multiple prediction tools

#### ORF Merging Capabilities
- **Commits**: 04d5782, 062f197, a82bdf7
- Added advanced ORF merging options supporting multiple configuration modes:
  - Tolerance-based merging
  - Overlapping group-based merging
  - Configurable merging parameters
- Updated minimum amino acid length parameter from default to 24 for improved ORF merging accuracy
- Provides flexibility in handling redundant ORF predictions

#### Unified ORF Container Support
- **Commit**: 0d90c93
- Added support for unified ORF container with improved examples
- Updated documentation to reflect container usage
- Enables containerized execution of ORF analysis workflows

#### Custom Container Support for ORF Analysis
- **Commit**: 9fc13db
- Implemented custom container support for ORF prediction and classification processes
- Optimized dependency management during container execution
- Allows users to specify custom containers for specialized workflows

#### Optimized ORF Prediction Process
- **Commit**: b9f4fb7
- Updated UNIFY_ORF_PREDICTIONS process to use biopython container
- Optimized dependency installation workflow
- Streamlined execution of ORF prediction analysis

#### Parallelization Enhancements
- **Commit**: 99012bd
- Parallelized unify_orf_prediction functionality
- Improved performance for large-scale ORF analysis
- Enables efficient processing of multiple samples

### Documentation

#### Statistical Models Documentation
- **Commit**: e2c95a3
- Added comprehensive statistical models documentation for Ribo-seq ORF detection tools
- Provides reference materials for understanding ORF detection methodologies
- Facilitates transparency in statistical approaches used by the pipeline

### Bug Fixes

#### ORFquant Output Filename Mismatch
- **Commit**: 59989d0
- Fixed issue with ORFquant output filenames not matching expected naming convention
- Ensures proper file identification and downstream processing

#### Dependency Installation Errors
- **Commit**: 1eaf0a0
- Fixed import error messages in unify_orf_predictions.py
- Improved error messaging to guide users in installing required dependencies
- Provides clear instructions when dependencies are missing

#### Input Connection Optimization
- **Commit**: 8b5cf17
- Optimized pre-filter and post-filter input connection method
- Fixed to ensure connections through sample IDs rather than other attributes
- Resolves data flow issues between analysis stages

### Infrastructure and Configuration Updates

#### Environment Configuration Cleanup
- **Commits**: d6f7ef9, 6baff64
- Removed extract_rl_cutoff environment configuration file
- Removed prepare_for_orfquant environment configuration file
- Updated main workflow to use RiboseQC environment configuration
- Updated main workflow to use ORFQUANT_RUN environment configuration
- Simplifies configuration management by consolidating environment definitions

#### Default Preparation Script Updates
- **Commit**: ed49afb
- Added skip_rpbp as default option in preparation scripts
- Reduces configuration overhead for standard workflows

#### Container Dependency Updates
- **Commit**: fbbc03a
- Added GenomicFeatures package to ORFquant Singularity container definition
- Ensures all required genomics analysis packages are available

### Workflow Preparation and Documentation
- **Commits**: f7d9bea, 9b918eb
- Updated prepare_workflow.py script with enhanced documentation
- Added detailed explanations for species and genome selection options
- Included multiple species usage examples
- Added thread option support for SRA conversion processes
- Optimized SRA to FASTQ conversion workflow

---

## January 31, 2026

### Features

#### Enhanced Workflow Preparation Tools
- Updated prepare_workflow.py with improved usability
- Added support for multiple species configurations
- Implemented threading options for SRA data conversion

### Bug Fixes and Improvements

#### ORFquant Processing
- **Commits**: 8e25bd3, f9b3600, 71dc017
- Fixed ORFquant parallelization issues
- Resolved ORFquant execution errors
- Corrected old variable name references

#### RiboseQC Data Flow
- **Commit**: f9406e6
- Fixed data flow issues in RiboseQC analysis
- Ensures proper data transmission between pipeline stages

#### Workflow Logic Refinement
- **Commit**: 45c2fd7
- Reorganized and clarified RiboseQ and ORFquant workflow logic
- Improves maintainability and clarity of pipeline structure

#### Workflow Automation
- **Commit**: f2c26c8
- Added new preparation workflow script
- **Commit**: b10d000
- Updated preparation workflow scripts
- Streamlines workflow initialization and configuration

---

## January 30, 2026

### Major Features

#### GENCODE ORF Annotation Workflow
- **Commit**: 6f8ac01
- Added complete GENCODE ORF annotation workflow
- Integrated unified ORF prediction capabilities
- Enables annotation of ORF predictions using GENCODE reference data

#### Pre-filtered QC Analysis
- **Commit**: 7708d3b
- Added pre-filter QC analysis option
- Supports comparative analysis between filtered and unfiltered BAM files
- Updated documentation to reflect feature
- Redesigned output file structure for comprehensive QC reporting

### Documentation and Output Structure

#### Output Directory Layout Redesign
- **Commit**: 8cbbe2d
- Redesigned complete output directory structure
- Updated documentation to reflect RiboseQC use of filtered BAM files
- Added descriptions for experimental output files for unified ORF prediction and classification
- Improves organization and accessibility of analysis results

### Bug Fixes

#### Null Value Handling
- **Commits**: a31e5b3, 1267eb3, 2ec0ac5
- Fixed multiple null value handling issues throughout pipeline
- **Commit**: f8d4437
- Updated default parameters to enable unified ORF prediction and classification
- Fixed related logic to properly handle null values
- Ensures stable execution with all parameter combinations

#### Lexical and Syntax Corrections

##### Bash Variable Escaping
- **Commit**: 0119c0b
- Escaped Bash variables in script blocks to prevent Groovy variable interpretation
- Fixed issues with variables like `$mode` being misinterpreted

##### AWK Field Separator Unification
- **Commit**: c8953a8
- Unified awk field separator (FS) and output field separator (OFS) to explicit tab characters
- Prevents misinterpretation as spaces or string errors
- Ensures consistent field parsing across scripts

##### Versions File Fixes
- **Commits**: 44f13d7, 70fbd90, 05f2b18
- Converted backticks in versions.yml to `$(...)` syntax
- Properly escaped all dynamic expressions
- Fixed Groovy lexical parsing errors in main.nf
- Ensures correct version tracking and reporting

#### File Path Handling
- **Commit**: d7c1bd8
- Corrected file paths to ensure workflow.projectDir references for scripts
- **Commit**: c01a3f1
- Optimized file path processing
- Removed unnecessary Channel.value wrapping
- Improves performance and clarity of file handling

#### Data Processing Issues
- **Commit**: bdfa3cc
- Fixed file object and filename list passing logic
- Ensures proper data structures throughout processing pipeline

#### String Parsing
- **Commit**: b6a0d21
- Resolved triple-quote string interpretation issues
- **Commit**: e453446
- Fixed pip installation failures

#### Dependency Management
- **Commit**: 94a0e7c
- Added patched packages to workdir
- Ensures custom package versions are available during execution

---

## Performance Optimization: ORFquant Series

### ORFquant v1.2.0 Baseline
- **Commit**: 9ba6649
- Established baseline version prior to optimization efforts
- Reference point for performance improvements

### ORFquant v1.3.0 - Performance Optimization
- **Commits**: 124312a, 7fd2fbb, 7909335
- **DPSS Caching Optimization**: 3-40x speedup achieved
- Complete performance optimization implementation
- Updated Singularity container definition
- Enables significantly faster ORF quantification for large datasets

### ORFquant v1.3.1 - Parallel Processing
- **Commits**: cb97685, c846ed0, e90b919
- Enabled parallel processing with FaFile improvements
- Incorporated FaFile parallel fix for concurrent file access
- Updated Singularity definition for parallel support
- Further enhances performance for multi-threaded execution

---

## Summary Statistics

| Metric | Value |
|--------|-------|
| Total Commits | ~60 |
| Date Range | 4 days (Jan 29 - Feb 2, 2026) |
| Major Features Added | 5+ |
| Bug Fixes | 20+ |
| Performance Optimizations | 2 (ORFquant v1.3.0 & v1.3.1) |
| Performance Improvement | 3-40x speedup |
| Documentation Updates | 3+ |
| Container Updates | 2+ |

---

## Key Highlights

### 🚀 Performance
- **ORFquant v1.3.0-1.3.1**: DPSS caching optimization achieving 3-40x speedup
- **Parallelization**: Unified ORF prediction parallelization for efficient multi-sample processing
- **FaFile Parallel Processing**: Enabled concurrent file access for improved throughput

### 🛡️ Robustness
- Comprehensive input file validation
- Enhanced error handling and reporting
- Improved null value handling throughout pipeline
- Better variable escaping and lexical parsing

### ✨ Functionality
- Complete unified ORF prediction pipeline integration
- GENCODE ORF annotation workflow
- Pre-filtered QC analysis capabilities
- Advanced ORF merging with multiple configuration options

### 🔧 Flexibility
- Custom container support for dependency management
- Unified container support for standardized execution
- Configurable ORF merging parameters
- Multiple species and genome support in preparation tools

### 📊 Quality
- Pre-filter QC analysis for comprehensive quality assessment
- Statistical models documentation for transparency
- Improved output directory structure for result organization
- Consolidated dependency management through environment configuration

---

## Breaking Changes

None identified in this release period.

---

## Migration Notes

- Users should update their ORFquant containers to v1.3.1 to benefit from performance improvements
- Output directory structure has been redesigned; scripts expecting old structure should be updated
- Some environment configuration files have been consolidated; verify custom configurations remain functional

---

## Known Issues

None reported at the time of this changelog.

---

## Next Steps

- Monitor performance of new unified ORF prediction pipeline in production environments
- Gather user feedback on new pre-filter QC analysis options
- Continue optimization of container execution workflows
- Expand statistical documentation as additional methodologies are refined

