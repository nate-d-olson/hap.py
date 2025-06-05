---
description: AI rules derived by SpecStory from the project AI interaction history
globs: *
---

# hap.py Copilot Instructions

## Project Overview
hap.py is a bioinformatics tool for benchmarking small variant calls.
The tool is a critical resource in the genomics community for evaluating the accuracy of variant callers.
The original codebase used Python 2, which is no longer supported, and had outdated dependencies.
This fork aims to modernize the codebase for continued use and development. When appropriate, refer to the initial implementation of [hap.py](https://github.com/Illumina/hap.py) when there is ambigutity in the functionality and logic.

## Repository Structure
- `src/`: Main source code directory
  - `hap_py/`: Core Python package (modernized from src/python)
  - `c++/`: Performance-critical algorithms implemented in C++ (Note: C++ code was largely replaced with Python for simplicity in the modernized version.)
  - `sh/`: Shell scripts for testing and utility functions
  - `data/`: Reference data files
- `external/`: External dependencies (e.g. htslib, rtg-tools)
- `example/`: Test data and usage examples
- `tests/`: Unit and integration test suite
- `scripts/`: Development and build scripts
- `.github/instructions/`: Instructions for project implementation and debugging
- `doc/`: Documentation files

## Current Project Status
✅ **Completed Modernization Tasks:**
- Python 3 conversion complete
- Modern Python package structure with pyproject.toml
- All unit and integration tests converted to pytest framework
- Type hints and Google-style docstrings added
- Error handling and logging implemented
- Pre-commit hooks and code quality tools configured
- Fixed `normalize_variant` method
- Fixed missing FILTER detection in VCF header checks
- Fixed RTG tool detection in tests
- Fixed test expectations in `findVCFEval` test
- Type annotation fix
- Fixed SDF Template Directory Creation (2025-05-27)
- Enhanced RTG Path Detection (2025-05-27)
- Fixed RTG Detection Warning (2025-05-27)
- Fixed Test Package Structure (2025-05-27)
- Verified Binary Files Exist (2025-05-27)
- Fixed MultiSampleQuantifier inconsistency by aliasing `register_sample` to `add_sample` to resolve test failures.

🔄 **In Progress:**
- C++ component modernization (Note: C++ code was largely replaced with Python for simplicity in the modernized version.)
- Performance optimization
- Documentation updates
- Implementation of the `quantify` package functionality
    - ✅ **Completed Analysis**: Comprehensive review of current quantify module implementation status
    - ✅ **Architecture Review**: Analyzed original C++ vs modernized Python implementation differences
    - ✅ **Gap Analysis**: Identified missing critical components including `_match_variants` method
    - ✅ **Development Plan**: Created comprehensive implementation roadmap (see `QUANTIFY_IMPLEMENTATION_PLAN.md`)
    - ✅ **Phase 1**: Core variant matching implementation (`_match_variants` method, benchmarking decision tracking)
        - ✅ Fixed critical test failures related to allele compatibility, variant classification, and performance.
        - ✅ Enhanced method compatibility to properly handle both pandas Series and dictionary inputs
        - ✅ Added robust tests for sophisticated variant matching algorithms
        - ✅ Established realistic performance expectations for current implementation
        - Phase 1 is now considered complete with core variant matching functionality working and tested.
    - ✅ **Phase 2**: Enhanced ROC analysis with confidence intervals and quality score stratification
        - Implementation includes methods: `_perform_roc_analysis()`, `_perform_quality_stratification()`, `_generate_roc_curve()`, `_calculate_bootstrap_confidence_intervals()`, `_perform_multi_threshold_analysis()`, and `_write_roc_results()`.
    - ✅ **Phase 3**: Superlocus analysis and region-based quantification
        - ✅ Fixed MultiSampleQuantifier inconsistency by aliasing `register_sample` to `add_sample` to resolve test failures.
        - Phase 3 is now considered complete with all core functionality working and tested.
    - ✅ **Phase 5**: GA4GH compliance and standards support
        - Implementation includes comprehensive GA4GH compliance classes (`GA4GHFormatter`, `GA4GHStratification`, `GA4GHMetrics`) and integration module
        - Added unit tests in `tests/unit/test_ga4gh_compliance.py` and integration tests in `tests/integration/test_ga4gh_integration.py`
        - Created validation script in `test_ga4gh_implementation.py`
        - Added detailed documentation in `PHASE5_IMPLEMENTATION_SUMMARY.md` and `GA4GH_IMPLEMENTATION_DETAILS.md`
        - Ensures compliance with the GA4GH benchmarking standards
        - Phase 5 is now considered complete with all GA4GH functionality implemented and documented.
    - 🚧 **Phase 4**: Performance optimization for large datasets (will be revisited after Phase 5)
    - Ensure the original functionality of the quantify module is maintained in the updated (modernized codebase).
    - Evaluate whether any of the original cython provides significant performance improvements compared to the in progress updated implementation that would justify the additional layer of complexity for package maintenance, development, and install.

## Development Environment Setup

### Prerequisites
- Python 3.8+ (recommended: 3.11)
- CMake 3.10+
- C++ compiler (GCC 7+ or Clang 10+)
- Git
- micromamba (recommended) or conda/mamba
- Standard bioinformatics tools: bcftools, samtools, tabix
- `pybedtools` (optional, for Phase 3 quantify implementation)

### Initial Setup

1. **Clone and navigate to the repository:**
```bash
git clone <repository-url>
cd hap.py
```

2. **Create and activate the development environment:**
```bash
# Using micromamba (RECOMMENDED for this project)
micromamba create -n happy-dev python=3.11
micromamba activate happy-dev

# Alternative: using venv
python -m venv .venv
source .venv/bin/activate  # On Windows: .venv\Scripts\activate
```

**IMPORTANT:** For all development work, use the `happy-dev` micromamba environment:
```bash
micromamba activate happy-dev
```

3. **Install development dependencies:**
```bash
# Install the package in development mode with all dependencies
pip install -e '.[dev,cpp]'

# Install pre-commit hooks for code quality
pre-commit install
```

4. **Build external dependencies:**
```bash
# External dependencies are now managed through the build system
# RTG tools and other dependencies are automatically configured
cmake -B build -S .
cmake --build build
```

5. **Configure environment variables:**
```bash
# Add to your shell profile (.bashrc, .zshrc, etc.)
export HGREF="/path/to/reference/genome.fa"  # Optional: for testing
```

### Code Quality Tools

The project uses modern Python development tools:

- **Black**: Code formatting (88 character line limit)
- **Ruff**: Fast linting and Python 3 compatibility checks
- **isort**: Import sorting compatible with Black
- **mypy**: Optional static type checking
- **pre-commit**: Automated checks on commit

### Running Code Quality Checks

```bash
# Format code with Black
black src/ tests/

# Check and fix linting issues
ruff check src/ tests/ --fix

# Sort imports
isort src/ tests/

# Type checking (optional)
mypy src/hap_py/

# Run all pre-commit hooks
pre-commit run --all-files
```

## Testing

### Test Structure
- `tests/unit/`: Unit tests for individual modules (✅ Converted to pytest)
- `tests/integration/`: End-to-end integration tests (✅ Converted to pytest)
- `tests/utils.py`: Shared test utilities
- `conftest.py`: Pytest configuration and fixtures

### Running Tests

**Always activate the environment first:**
```bash
micromamba activate happy-dev
```

**Unit Tests:**
```bash
# Run all unit tests
pytest tests/unit/ -v

# Run specific test file
pytest tests/unit/test_vcfeval.py -v

# Run with coverage
pytest tests/unit/ --cov=hap_py --cov-report=html
```

**Integration Tests:**
```bash
# Run all integration tests (requires external tools)
pytest tests/integration/ -v

# Run specific integration test
pytest tests/integration/test_performance.py -v

# Skip slow tests
pytest tests/integration/ -v -m "not slow"
```

**All Tests:**
```bash
# Run complete test suite
pytest tests/ -v

# Run tests in parallel (if pytest-xdist installed)
pytest tests/ -n auto
```

### Test Markers
- `@pytest.mark.integration`: Integration tests requiring external tools
- `@pytest.mark.slow`: Long-running tests
- `@pytest.mark.cpp`: Tests requiring C++ components (Note: These tests may need adjustment or skipping as C++ components were largely replaced with Python.)

### Common Test Issues During Modernization

**1. Path-based Import Errors:**
```python
# Old style (Python 2)
import haplo.vcfeval

# New style (Python 3)
import hap_py.haplo.vcfeval
```

**2. Module Not Found Errors:**
Check and fix:
- Package structure in `src/hap_py/`
- `__init__.py` files in all package directories
- Import statements using relative vs absolute imports
- `sys.path` modifications in test files

**3. String/Bytes Compatibility:**
```python
# Use utility functions for Python 3 compatibility
from hap_py.haplo.string_handling import ensure_str, ensure_bytes

# Handle both string and bytes inputs
def process_sequence(seq):
    seq = ensure_str(seq)  # Convert bytes to str if needed
    return seq.upper()
```

**4. File Path Issues:**
```python
# Use pathlib for cross-platform compatibility
from pathlib import Path

# Old way
import os
test_file = os.path.join(os.path.dirname(__file__), "data", "test.vcf")

# New way
test_file = Path(__file__).parent / "data" / "test.vcf"
```

**5. External Tool Dependencies:**
```python
# Check tool availability before running tests
import shutil
if not shutil.which("bcftools"):
    pytest.skip("bcftools not available")

# Use proper RTG tools path
rtg_path = Path("build/external/rtg-tools/rtg")
if not rtg_path.exists():
    pytest.skip("RTG tools not built")
```

### Debugging Test Failures

1. **Check Python path setup:**
```python
import sys
print("Python path:", sys.path)
print("Current working directory:", os.getcwd())
```

2. **Verify package installation:**
```bash
pip list | grep hap
python -c "import hap_py; print(hap_py.__file__)"
```

3. **Check external dependencies:**
```bash
# Verify build artifacts
ls build/external/

# Check other tools
bcftools --version
samtools --version
```

4. **Run tests with verbose output:**
```bash
pytest tests/unit/test_specific.py -v -s --tb=long
```

## Build and Installation

### Getting the Version from `pyproject.toml`

To have the install process use the version from `pyproject.toml`, configure the build system to read the version dynamically. The recommended approach is using `setuptools_scm` as it automatically manages versions based on Git tags.

**Steps:**

1.  **Update `pyproject.toml`:**

```toml
[build-system]
requires = ["setuptools>=61.0", "setuptools_scm[toml]>=6.2"]
build-backend = "setuptools.build_meta"

[project]
name = "hap-py"
dynamic = ["version"]
description = "Haplotype VCF comparison tools"
# ... other metadata

[tool.setuptools_scm]
write_to = "src/hap_py/_version.py"
fallback_version = "0.4.0"
version_scheme = "post-release"
local_scheme = "dirty-tag"

[tool.setuptools.packages.find]
where = ["src"]
```

2.  **Create `src/hap_py/_version.py`:**

This file will be automatically generated and updated by `setuptools_scm`.

3.  **Import the version in `src/hap_py/__init__.py`:**

```python
try:
    from ._version import version as __version__
except ImportError:
    # Fallback for development installs
    __version__ = "0.4.0"
```

4.  **Tag your current version:**

```bash
git tag v0.4.0
git push origin v0.4.0
```

5.  **Install in development mode:**

```bash
pip install -e .
```

6.  **Verify Installation:**

```bash
# Check the version
python -c "import hap_py; print(hap_py.__version__)"

# Or check from command line if you have a CLI
python -m hap_py --version
```

### Development Installation
```bash
# Install in development mode (changes reflected immediately)
pip install -e .[dev]
```

### Production Build
```bash
# Build source distribution and wheel
python -m build

# Install from wheel
pip install dist/hap.py-*.whl
```

### CMake Build (for C++ components)
```bash
cmake -B build -S .
cmake --build build --config Release
```

## Development Plan

### ✅ Completed: Modernization
- ✅ Add type hints to improve code reliability
- ✅ Update documentation with Google-style docstrings
- ✅ Implement proper package structure with pyproject.toml
- ✅ Convert shell script tests to pytest framework
- ✅ Implement proper error handling and logging
- ✅ Fixed `normalize_variant` method
- ✅ Fixed missing FILTER detection in VCF header checks
- ✅ Fixed RTG tool detection in tests
- ✅ Fixed test expectations in `findVCFEval` test
- ✅ Type annotation fix
- ✅ Fixed SDF Template Directory Creation (2025-05-27)
- ✅ Enhanced RTG Path Detection (2025-05-27)
- ✅ Fixed RTG Detection Warning (2025-05-27)
- ✅ Fixed Test Package Structure (2025-05-27)
- ✅ Verified Binary Files Exist (2025-05-27)
- ✅ Fixed MultiSampleQuantifier inconsistency by aliasing `register_sample` to `add_sample` to resolve test failures.

### ✅ Completed: Phase 3 Validation
- The Phase 3 implementation is now complete and ready for submission. The complete set of features has been implemented, tested, and documented according to the project requirements.

### ✅ Completed: Phase 5 GA4GH Compliance
- The Phase 5 implementation is now complete and documented. The GA4GH compliance functionality for the hap.py project has been successfully implemented. This implementation includes comprehensive GA4GH compliance classes, integration with the QuantifyEngine, unit and integration tests, and detailed documentation.

### 🔄 In Progress: C++ Modernization and Optimization
- Update C++ code to use modern standards (Note: C++ code was largely replaced with Python for simplicity in the modernized version.)
- Optimize memory usage for large genomic datasets
- Improve parallelization for performance
- 🔄 Implementation of the `quantify` package functionality
    - ✅ **Completed Analysis**: Comprehensive review of current quantify module implementation status
    - ✅ **Architecture Review**: Analyzed original C++ vs modernized Python implementation differences
    - ✅ **Gap Analysis**: Identified critical components including sophisticated variant matching algorithms
    - ✅ **Development Plan**: Created comprehensive implementation roadmap (see `QUANTIFY_IMPLEMENTATION_PLAN.md`)
    - ✅ **Phase 1**: Core variant matching implementation (`_match_variants` method, benchmarking decision tracking)
        - ✅ Fixed critical test failures related to allele compatibility, variant classification, and performance.
        - ✅ Enhanced method compatibility to properly handle both pandas Series and dictionary inputs
        - ✅ Added robust tests for sophisticated variant matching algorithms
        - ✅ Established realistic performance expectations for current implementation
        - Phase 1 is now considered complete with core variant matching functionality working and tested.
    - ✅ **Phase 2**: Enhanced ROC analysis with confidence intervals and quality score stratification
        - Implementation includes methods: `_perform_roc_analysis()`, `_perform_quality_stratification()`, `_generate_roc_curve()`, `_calculate_bootstrap_confidence_intervals()`, `_perform_multi_threshold_analysis()`, and `_write_roc_results()`.
    - ✅ **Phase 3**: Superlocus analysis and region-based quantification
        - Phase 3 is now considered complete with all core functionality working and tested.
    - ✅ **Phase 5**: GA4GH compliance and standards support
        - Implementation includes comprehensive GA4GH compliance classes (`GA4GHFormatter`, `GA4GHStratification`, `GA4GHMetrics`) and integration module
        - Added unit tests in `tests/unit/test_ga4gh_compliance.py` and integration tests in `tests/integration/test_ga4gh_integration.py`
        - Created validation script in `test_ga4gh_implementation.py`
        - Added detailed documentation in `PHASE5_IMPLEMENTATION_SUMMARY.md` and `GA4GH_IMPLEMENTATION_DETAILS.md`
        - Ensures compliance with the GA4GH benchmarking standards
        - Phase 5 is now considered complete with all GA4GH functionality implemented and documented.
    - 🚧 **Phase 4**: Performance optimization for large datasets (will be revisited after Phase 5)
    - Ensure the original functionality of the quantify module is maintained in the updated (modernized codebase).
    - Evaluate whether any of the original cython provides significant performance improvements compared to the in progress updated implementation that would justify the additional layer of complexity for package maintenance, development, and install.

### 📋 Future: Infrastructure and Documentation
- Add CI/CD pipeline for testing
- Create containerized deployment options
- Update user documentation

## Troubleshooting Common Issues

### Import Errors
1. Ensure virtual environment is activated
2. Check that package is installed in development mode: `pip install -e .`
3. Verify Python path includes the src directory
4. Check for missing `__init__.py` files

### External Tool Errors
1. Verify build completed: `ls build/external/`
2. Check PATH includes external tools
3. Ensure required bioinformatics tools are installed

### Test Data Issues
1. Check that example data exists: `ls example/`
2. Verify file permissions and formats
3. Ensure reference files are properly indexed

### Build Issues
1. Update CMake: minimum version 3.10 required
2. Check C++ compiler compatibility
3. Verify external dependencies built successfully
4. Clear build directory and rebuild: `rm -rf build && cmake -B build -S .`

## Workflow & Release Rules

### Environment Management

**CRITICAL:** Always start terminal sessions with:
```bash
micromamba activate happy-dev
```

This environment provides:
- Python 3.11.12
- All required dependencies for hap.py development
- Isolated environment preventing conflicts with other projects

**Environment Verification:**
```bash
# Verify correct environment is active
which python  # Should show: /Users/nolson/micromamba/envs/happy-dev/bin/python
python --version  # Should show: Python 3.11.12
```

**For all development tasks:**
- Running tests: `micromamba activate happy-dev && pytest tests/`
- Code formatting: `micromamba activate happy-dev && black src/`
- Building: `micromamba activate happy-dev && pip install -e .`

- When fixing failing tests, always clean up any temporary files created and used during the debugging process.
- Tests should accurately reflect tool functionality; failures should indicate errors in the source code, not the tests themselves.
- Use the micromamba environment `happy` (`micromamba activate happy`) for an isolated Python environment and installation. **NOTE:** Use `micromamba activate happy-dev` instead of `happy`.
- When debugging tests:
    - First run the unit tests to identify failures.
    - Then proceed with debugging the integration tests.
    - Save test output to a file for analysis.
    - Start by examining the first failing test.
- When addressing "command not found" errors during testing, ensure the correct path to the executable is used.
    - Example: If `rtg` executable is available at `/Users/nolson/hap.py-modern-claude4/hap.py/external/rtg-tools-3.12.1/rtg`, modify `vcfeval.py` to use this path.
- When running integration tests and encountering "rtg: command not found" errors, explicitly pass the `--engine-vcfeval-path` argument to the `hap.py` call within the test to ensure the correct `rtg` location is used.
- The original hap.py code base is available at [https://github.com/Illumina/hap.py](https://github.com/Illumina/hap.py). Note that the C++ code from this repository was largely replaced with Python for simplicity in the modernized version.
- The tests should use relative paths or environment variables to specify file locations, making the tests more portable.
- Employ the `pathlib` module for handling file paths in a cross-platform manner, ensuring that the tests work correctly on different operating systems.
- If a test fails due to a logic error in the test itself (rather than a bug in the tool), correct the test logic.
- When committing changes, clean up the codebase by removing extra files generated during coding sessions, renaming files as it makes sense, or restoring deleted files.
- When committing changes, clean up the codebase a bit. removing extra files generated during the coding sessions, renaming files as it makes sense, or restoring deleted files.
- When addressing integration test errors, analyze the integration test log file to understand the failures and identify common error patterns before developing a comprehensive fix strategy.
- When integration tests are failing with "rtg: command not found" errors, ensure the `--engine-vcfeval-path` argument is correctly passed to the `hap.py` call within the test.
- If the `--engine-vcfeval-path` argument is being used, and the tests are still failing, examine the `src/hap_py/haplo/tool_runner.py` and `src/hap_py/haplo/vcfeval.py` files to ensure that the path is correctly utilized by the `ToolRunner` class when finding the `rtg` executable.
- If the tests are still failing and `rtg vcfeval` commands are not executing correctly, check how `hap.py` parses the `--engine-vcfeval-path` argument and how it's passed to the `runVCFEval` function in `vcfeval.py`.
- **When debugging integration tests, review `integration_test_output.txt` to identify failing tests and their causes.**
- **Examine the corresponding test files in `tests/integration/` and their fixtures/utilities.**
- **Check the example/reference data in `example/integration/` for correctness and completeness.**
- **Trace failures to the source code in `src/hap_py/` and fix bugs or update tests as needed.**
- **Use the `happy` micromamba environment and ensure all dependencies (like RTG) are available.** **NOTE:** Use `micromamba activate happy-dev` instead of `happy`.
- If the test is hanging, consider that it might be related to temporary directory permissions or cross-filesystem issues. Try using a temporary directory within the repo as a test.
- When running integration tests, capture the actual error output to see what's going wrong.
- If tests are still failing, check if you need to restart Python or reinstall the package.
- To ensure RTG tools are accessible during testing, set up the RTG path by creating a symbolic link (e.g., `ln -s external/rtg-tools-3.12.1/rtg .`) and including it in the PATH environment variable (e.g., `export PATH=$PATH:$(pwd)/`).
- When running tests, ensure RTG is in the PATH: `export PATH=$PATH:$(pwd)/`.
- When running tests, make sure to include rtg in path e.g. `export PATH=$PATH:$(pwd)/` after creating a symbolic link for rtg e.g. ` ln -s external/rtg-tools-3.12.1/rtg .`
- Integration tests should only fail when the output files don't match the equivalent expected files, and not when there are additional output files than the expected data files.
- When analyzing integration test failures, remember that tests can fail because the `hap.py` command itself is returning non-zero exit status, not always because of file comparison issues. Examine the error messages to understand the root cause.
- **Integration tests MUST NOT fail when there are additional output files than the expected data files. Tests should only fail when the output files don't match the equivalent expected files.**
- Before committing any changes to the repository, the Phase 3 validation testing MUST be finalized and debugged, and the relevant documentation updated.
- All changes to the repository while implementing phase 3 must be committed and pushed to github.
- Phase 4 (performance optimization) will be revisited after we have a functioning codebase and are able to run a whole genome callset end ot end without error.
- The changes for Phase 3 should be committed and pushed to GitHub.
- The `PHASE5_IMPLEMENTATION_PLAN.md` file should be added to git and finalized.
- Changes to the `test_phase3_implementation.py` and `validate_phase3_complete.py` files should be committed.
- All changes to the repository while implementing phase 3 MUST be committed and pushed to github.
- Use the `TEST_STATUS.md` file to track the current status of the tests.
- Implement fixes according to the `TEST_ERROR_ANALYSIS_AND_FIXING_PLAN.md` file.

### Detailed Plan for Debugging Failing Integration Tests

This outlines a strategic approach to debug and resolve failing integration tests, incorporating specific steps and considerations based on previous findings and project guidelines.

#### 1. Preparation and Environment Setup

- **Activate the `happy-dev` environment:**
  ```bash
  micromamba activate happy-dev
  ```
  This ensures all necessary dependencies are available and isolated, preventing conflicts with other projects or system-level packages.

- **Verify RTG tools availability:**
  Confirm that the RTG tools are accessible at the expected path: `/Users/nolson/hap.py-modern-claude4/hap.py/external/rtg-tools-3.12.1/rtg`.

#### 2. Identification of Failing Tests

- **Review `integration_test_output.txt`:**
  This file contains the output and error messages from the most recent integration test run. It lists which tests failed, the error types, and any stack traces or assertion errors. This is the primary resource for diagnosing the failures.

- **Focus on key integration tests:**
  Pay close attention to tests in the `tests/integration/` directory, such as `test_integration.py`, `test_happy_pg.py`, and `test_giab.py`, which cover core functionalities.

#### 3. Common Failure Causes and Resolutions

##### a. RTG/VCFEval Not Found Errors

- **Explicitly specify RTG path:**
  Ensure that all test calls to `hap.py` that use the `--engine=vcfeval` option also include the `--engine-vcfeval-path` argument, pointing to the correct RTG executable:
  ```
  --engine-vcfeval-path /Users/nolson/hap.py-modern-claude4/hap.py/external/rtg-tools-3.12.1/rtg
  ```
  Update the relevant integration tests to pass this argument explicitly.

##### b. Temporary File Cleanup

- **Implement cleanup logic:**
  Verify that all tests properly clean up any temporary files created during execution. This is crucial to prevent interference between tests and to avoid filling up disk space. Use `tmp_path` fixtures provided by pytest for managing temporary directories.

- **Add teardown/cleanup functions:**
  If cleanup logic is missing in any test, add it to ensure that temporary files are removed after the test completes.

##### c. Reference Data and File Paths

- **Ensure data accessibility:**
  Confirm that all test data files exist and are accessible to the tests. Use relative paths or environment variables to specify file locations, making the tests more portable.

- **Use `pathlib` for cross-platform compatibility:**
  Employ the `pathlib` module for handling file paths in a cross-platform manner, ensuring that the tests work correctly on different operating systems.

##### d. Test Logic vs. Source Code Issues

- **Distinguish between test errors and source code bugs:**
  If a test fails due to a logic error in the test itself (rather than a bug in the tool), correct the test logic. If the failure is caused by a bug in the tool's source code, fix the underlying issue in the relevant module.

#### 4. Debugging Workflow

- **Run unit tests first:**
  Execute all unit tests using `pytest tests/unit/ -v` to identify any low-level issues before running the integration tests.

- **Run integration tests with output capture:**
  Execute the integration tests using `pytest tests/integration/ -v | tee integration_test_output.txt`. The `tee` command captures the output to both the console and the `integration_test_output.txt` file, allowing for detailed analysis.

- **Systematic debugging:**
    - Start with the first failing integration test.
    - Carefully read the error message and traceback to understand the cause of the failure.
    - If the error indicates a missing command or file, verify the path and dependency setup.
    - If the error indicates a logic or output mismatch, compare the actual output with the expected output to identify discrepancies.
    - Clean up any temporary files created during the debugging process.

#### 5. Implementation of Fixes

- **Update test files:**
  Modify the test files to pass the correct `--engine-vcfeval-path` argument where necessary.

- **Add/verify cleanup logic:**
  Implement or verify cleanup logic for temporary files in each test.

- **Correct test logic:**
  Fix any errors in the test logic that are causing the tests to fail.

- **Fix source code bugs:**
  If a bug is found in the source code, correct it in the relevant module.

#### 6. Code Quality and Commit Process

- **Run pre-commit hooks:**
  Execute the pre-commit hooks using `pre-commit run --all-files` to ensure that the code is properly formatted and linted before committing.

- **Document changes:**
  Document any new test cases or fixes that were implemented.

#### 7. Re-run and Verification

- **Re-run the integration test suite:**
  After implementing the fixes, re-run the full integration test suite to ensure that all tests pass and that no temporary files are left behind.

- **Confirm all tests pass:**
  Verify that all tests pass and that there are no unexpected errors or warnings.

- **Clean environment:**
  Ensure the environment is clean and ready for further development or testing.

### Patching RTG Template Directory Creation

When running `rtg format`, the command fails if the SDF template directory already exists. The `tempfile.NamedTemporaryFile` creates the directory as a file, not a directory, and RTG refuses to overwrite it.

**Solution:**
- Replace the use of `tempfile.NamedTemporaryFile` for the SDF template with `tempfile.mkdtemp`, which creates a unique directory and can be safely removed and recreated as needed.

**Implementation:**
1. **Locate this block in `runVCFEval` in `vcfeval.py`:**
   ```python
   try:
       with tempfile.NamedTemporaryFile(
           dir=args.scratch_prefix, prefix="vcfeval.sdf", suffix=".dir"
       ) as stf:
           template_dir = stf.name

       # Remove template dir if it already exists (RTG format will fail otherwise)
       if os.path.exists(template_dir):
           logging.warning(f"SDF template directory {template_dir} already exists. Removing it before running rtg format.")
           shutil.rmtree(template_dir)
       os.makedirs(template_dir, exist_ok=True)
   ```

2. **Replace it with:**
   ```python
   try:
       # Use mkdtemp to create a unique directory for the SDF template
       template_dir = tempfile.mkdtemp(dir=args.scratch_prefix, prefix="vcfeval.sdf.")
       # Remove template dir if it already exists (should not happen, but for safety)
       if os.path.exists(template_dir):
           shutil.rmtree(template_dir)
       os.makedirs(template_dir, exist_ok=True)
   ```

### Test Failure Analysis & Fix Plan

#### Issues Identified:

1.  **Unit Test Failures:**
    *   `test_normalize_variant`: Expected position 101 but got 102 (off-by-one error in variant normalization)
    *   Two `vcfeval.py` tests failing due to RTG executable path issues

2.  **Integration Test Issues:**
    *   RTG executable path not properly resolved
    *   Tests are hardcoded with the specific RTG path that may not be consistent
    *   VCF header validation errors: Missing FILTER field and duplicate FORMAT entries.
    *   `roc.tsv` output file is missing.

3.  **Path Configuration:**
    *   RTG tools are available at `rtg` but tests need proper path handling

### Updating the `init()` function

- The `init()` function in `__init__.py` should be updated to also check our custom RTG location, instead of only the PATH.
- Modify the `init()` function to properly detect our included RTG tools.

### Summary of Fixes Made (2025-05-27)

### 1. **Fixed SDF Template Directory Creation** ✅
- **Issue**: RTG `format` command was failing with "directory already exists" errors
- **Fix**: Simplified the `mkdtemp` logic in `vcfeval.py` to avoid problematic directory existence checks
- **File**: `vcfeval.py`

### 2. **Enhanced RTG Path Detection** ✅
- **Issue**: `findVCFEval()` wasn't properly locating RTG tools in the modernized directory structure
- **Fix**: Updated the function to check for RTG tools in `rtg` path
- **File**: `vcfeval.py`

### 3. **Fixed RTG Detection Warning** ✅
- **Issue**: "WARNING:root:Executable for rtg not found" messages appearing despite RTG being available
- **Fix**: Updated the `init()` function in `__init__.py` to check for our included RTG tools before issuing warnings
- **File**: `__init__.py`

### 4. **Fixed Test Package Structure** ✅
- **Issue**: Integration tests failing with `ModuleNotFoundError: No module named 'tests.utils'`
- **Fix**: Created missing `__init__.py` files to make tests a proper Python package
- **Files**:
  - `tests/__init__.py`
  - `tests/integration/__init__.py`

### 5. **Verified Binary Files Exist** ✅
- **Status**: Confirmed that all required binary files exist in the `build/bin` directory:
  - `hap.py` ✅
  - `hapenum` ✅
  - `hapcmp` ✅
  - `multimerge` ✅ (placeholder script)
  - `qfy.py
