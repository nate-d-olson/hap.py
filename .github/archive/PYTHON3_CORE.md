# Modernized hap.py: Core Functionality with Python 3

This branch contains a modernized version of hap.py that focuses on the core functionality:
running hap.py with vcfeval as the comparison engine and stratified performance metrics.

## Changes Made

1. **Python 3 Compatibility**:
   - Updated all core files to work with Python 3 (requires Python 3.6+)
   - Fixed string handling and encoding issues
   - Updated exception handling to use Python 3 syntax
   - Improved subprocess handling

2. **Streamlined Functionality**:
   - Removed non-core components (somatic variant calling, bamstats, xcmp and scmp engines)
   - Preserved removed functionality in separate feature branches
   - Focused on vcfeval as the primary comparison engine

3. **Build System Improvements**:
   - Updated CythonSupport.cmake for better NumPy and Cython detection
   - Created simplified requirements file for core functionality
   - Improved dependency handling in the install script

4. **Code Quality**:
   - Added type hints to improve code maintainability
   - Improved error handling and logging
   - Enhanced documentation with more descriptive docstrings

## Installation

### Requirements

- Python 3.6+
- C++ compiler with C++11 support
- CMake 3.5+
- Boost libraries
- (optional) RTG Tools for vcfeval

### Basic Installation

```bash
# Build hap.py with Python 3
python3 install.py /path/to/install/dir

# To run the core tests
./test_py3_core.sh
```

### Custom Installation

```bash
# Using a specific Python interpreter
python3 install.py /path/to/install/dir --python-interpreter=/path/to/python3

# Install with RTG Tools for vcfeval
python3 install.py /path/to/install/dir --with-rtgtools
```

## Usage

The basic usage remains unchanged, but now you should specify vcfeval as the engine:

```bash
hap.py truth.vcf query.vcf -r reference.fa -o output-prefix --engine vcfeval
```

For stratification by regions:

```bash
hap.py truth.vcf query.vcf -r reference.fa -o output-prefix --engine vcfeval \
    -f confident.bed --stratification strat.tsv
```

## Core Components

The modernized version maintains these key components:

1. **hap.py**: Main entry point for variant comparison
2. **vcfeval**: RTG's variant comparison engine
3. **Stratification**: Region-based performance stratification
4. **ROC Curves**: Score-based performance analysis

## Known Issues

- Python 3's handling of Unicode strings may cause issues with older VCF files
- Some features may behave differently due to changes in Python 3's behavior
- The preprocessing step can be slow for large VCF files due to file splitting

## Future Work

1. **Further Python 3 Modernization**:
   - Use pathlib for file operations
   - Implement more type hints for better IDE support
   - Add proper exception handling with context managers

2. **Performance Improvements**:
   - Optimize preprocessing for large VCF files
   - Improve parallelization for multi-sample VCFs

3. **Testing Framework**:
   - Convert shell script tests to pytest
   - Add more unit tests for core functionality
