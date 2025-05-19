---
mode: 'agent'
tools: ['githubRepo', 'codebase', 'use_tool', 'debuggingapproach', 'metacognitivemonitoring' ]
description: 'continue hap.py python3 migration'
---

# Continue Python 3 Migration: hap.py

## Migration Status Summary

Based on current progress (Updated May 17, 2025):

- **Overall Progress**: 85.7% of files fully migrated (42/49 files)
- **Remaining Issues**: 28 issues across 7 files
- **Priority Areas**:
  - Core functionality optimization (primary focus)
  - Testing and validation of Python 3 code
  - Documentation and deployment

### Recently Completed Tasks

The following major tasks have been successfully completed:

- **Core Functionality Streamlining**:
  - Focused on vcfeval comparison engine and stratified metrics
  - Removed non-core components (som.py, bamstats.py, scmp and xcmp engines)
  - Created feature branches to preserve removed functionality:
    - `feature/alternative-engines` for scmp and xcmp
    - `feature/somatic-support` for somatic variant calling

- **Python 3 Compatibility**:
  - Updated shebang lines from `python` to `python3`
  - Fixed string handling (bytes vs Unicode)
  - Fixed subprocess execution and error handling
  - Improved file I/O with proper encoding
  - Added improved type hints for better code clarity

- **Core Implementations**:
  - Enhanced `u_unhappy` and `v_vcfeval` functions in quantify.py
  - Added missing `run_quantify` function in quantify.py
  - Fixed Haplo module imports and dependencies

- **Build System Enhancements**:
  - Updated install.py to support Python 3
  - Added better dependency handling with fallback to different requirements files
  - Improved CythonSupport.cmake for NumPy and Cython detection

- **Testing and Documentation**:
  - Created test_py3_core.sh script for testing core functionality
  - Added update_shebangs.py script to update Python files
  - Created comprehensive documentation in PYTHON3_CORE.md
  - Added PYTHON3_MIGRATION_SUMMARY.md for quick reference

### Next Steps

1. **Testing and Validation (Priority)**:
   - Run comprehensive tests on the updated codebase
   - Create unit tests for core functionality
   - Verify output matches Python 2 version within acceptable tolerance
   - Add test cases for edge cases in string handling

   ```bash
   # Build test installation and run core tests
   ./test_py3_core.sh

   # Run focused tests on specific modules
   python3 -m pytest src/python/Haplo/tests/
   ```

2. **Performance Optimization**:
   - Profile and optimize critical paths
   - Improve memory usage for large VCF files
   - Enhance parallel processing capabilities
   - Optimize Python-C++ interface

   ```bash
   # Profile a typical run
   python3 -m cProfile -o happy_profile.prof src/python/hap.py [typical arguments]

   # Analyze profile
   python3 -m pstats happy_profile.prof
   ```

3. **Further Code Modernization**:
   - Complete type annotations throughout codebase
   - Use pathlib for file operations
   - Implement context managers for resource handling
   - Use f-strings consistently for string formatting

4. **CI/CD Setup**:
   - Create GitHub Actions workflow for automated testing
   - Implement code quality checks in the CI pipeline
   - Add test coverage reporting
   - Set up automatic documentation updates

5. **Documentation Updates**:
   - Expand user documentation with Python 3 specific information
   - Create migration guide for users
   - Update API documentation
   - Add examples for common use cases

## Migration Workflow

For continuing the migration work:

1. **Setup Development Environment**:

   ```bash
   # Clone the repository (if needed)
   git clone https://github.com/Illumina/hap.py.git
   cd hap.py

   # Create a Python 3 virtual environment
   python3 -m venv venv
   source venv/bin/activate

   # Install development dependencies
   pip install -r happy.core-requirements.py3.txt
   pip install pytest black ruff mypy
   ```

2. **Code Quality Tools**:

   ```bash
   # Format code with Black
   black src/python/

   # Check typing with mypy
   mypy src/python/

   # Use ruff for linting and auto-fixing
   ruff check --fix src/python/
   ```

3. **Testing Process**:

   ```bash
   # Build and install
   python3 install.py /tmp/happy-build

   # Run tests
   cd /tmp/happy-build
   src/sh/run_tests.sh
   ```

## Implementation Details

### Core Functionality Focus

The core functionality of hap.py has been streamlined to focus on:

1. **vcfeval Comparison Engine**:
   - Preserved and enhanced all vcfeval integration
   - Improved string handling in vcfeval result parsing
   - Fixed command execution and output processing

2. **Stratified Metrics**:
   - Maintained stratification functionality
   - Enhanced region handling and BED file processing
   - Improved ROC curve generation

3. **Python-C++ Integration**:
   - Updated Cython modules for Python 3
   - Added improved error handling
   - Created fallback implementations for testing

### Feature Branches

The following feature branches have been created to preserve removed functionality:

- `feature/alternative-engines`: Contains xcmp and scmp comparison engines
- `feature/somatic-support`: Contains somatic variant calling functionality
- `feature/bamstats-support`: Contains BAM statistics functionality

These can be used to restore functionality if needed or maintained as separate modules.
