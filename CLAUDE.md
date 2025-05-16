# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

hap.py is a toolkit for benchmarking variant calls against gold standard truth datasets. It performs genotype-level haplotype comparison for diploid samples, handles complex variant representations through graph-based analysis, provides variant preprocessing/normalization, and stratifies variant counts by type and genomic regions.

This repository is currently undergoing a Python 3 migration, focusing on core functionality using the vcfeval engine and stratified metrics.

## Build and Installation

### Python 3 Installation (recommended)

```bash
# Create and activate a Python 3 virtual environment
python3 -m venv venv_py3
source venv_py3/bin/activate

# Install Python dependencies
pip install -r happy.requirements.py3.txt

# Install using the Python 3 installer
python3 install_py3.py /path/to/install/dir
```

### Manual Build with CMake

```bash
# Configure and build the project
./configure.sh Release
make -j$(nproc)
make install

# Alternative simple build for Python 3
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release -DBUILD_PYTHON3=ON
make
make install
```

### Docker

```bash
docker pull pkrusche/hap.py
docker run -it -v `pwd`:/data pkrusche/hap.py /opt/hap.py/bin/hap.py [args]
```

## Running Tests

```bash
# Run core Python 3 functionality tests
./test_py3_core.sh

# Test Cython module loading
python3 test_cython_module_py3.py --build-dir /path/to/build/dir

# Run comprehensive Cython integration tests
./test_cython_integration_py3.sh

# Run all Python 3 tests
./run_all_py3_tests.sh
```

## Usage Examples

### Basic Variant Comparison

```bash
# Using Python 3
python3 /path/to/bin/hap.py.py3 truth.vcf query.vcf -r reference.fa -o output_prefix

# Using vcfeval engine
python3 /path/to/bin/hap.py.py3 truth.vcf query.vcf -r reference.fa -o output_prefix --engine vcfeval

# With confidence regions
python3 /path/to/bin/hap.py.py3 truth.vcf query.vcf -f confident.bed -r reference.fa -o output_prefix
```

### Preprocessing and Normalization

```bash
# Clean and normalize VCF files
python3 /path/to/bin/pre.py.py3 input.vcf -o output_prefix -r reference.fa
```

### Variant Quantification

```bash
# Run just the quantification step
python3 /path/to/bin/qfy.py.py3 truth.vcf query.vcf -r reference.fa -o output_prefix
```

## Python 3 Migration Status

The codebase is being migrated from Python 2 to Python 3. Key points:

1. Core functionality focus:
   - vcfeval comparison engine
   - Stratified performance metrics
   - Preprocessor and quantification steps

2. Non-core components moved to feature branches:
   - Somatic variant calling
   - bamstats
   - scmp/xcmp comparison engines

3. Migration progress:
   - Core utilities fully migrated
   - Tools module fully migrated
   - Haplo module partially migrated (~70%)

## Key Development Tasks for Python 3 Migration

### Checking Migration Issues

```bash
# Identify Python 3 compatibility issues
python3 check_py3_issues.py --dir src/python

# Fix truncated Python files
python3 fix_truncated_files.py
```

### Updating Cython Integration

```bash
# Update Cython modules for Python 3
python3 update_cython_modules_py3.py --src-dir src/python

# Update single module
python3 update_cython_for_py3.py src/python/path/to/module.pyx
```

## Architecture Overview

hap.py consists of several key components:

1. **Variant Reading and Preprocessing (`pre.py`)**:
   - Handles VCF loading, validation, and normalization
   - Applies variant decomposition (splitting multi-allelic variants)
   - Normalizes variants (left-shifting, redundant base removal)

2. **Comparison Engine**:
   - **vcfeval**: Primary engine, based on RTG's vcfeval
   - Performs best-matching between truth and query variants

3. **Quantification (`qfy.py`)**:
   - Produces precision/recall metrics
   - Stratifies results by variant type and regions
   - Generates output tables and VCFs

4. **C++ Core Library**:
   - Provides variant manipulation utilities
   - Implements graph-based comparison algorithms
   - Exposed to Python through Cython

The migration to Python 3 preserves this architecture while focusing on the core functionality.