---
description: AI rules derived by SpecStory from the project AI interaction history
globs: *
---

# hap.py CODEX Development Guide

## Project Overview
hap.py is a bioinformatics tool for benchmarking small variant calls, widely used for evaluating the accuracy of variant callers in the genomics community. The original codebase used Python 2 (now end-of-life) and had outdated dependencies. This fork aims to modernize the codebase for continued use and development.

## Repository Structure
- `src/`: Main source code directory
  - `hap_py/`: Core Python package (modernized from the original `src/python` code)
- `example/`: Example usage data and test data sets
- `tests/`: Unit and integration tests
- `scripts/`: Development and build scripts

## Project Status and Roadmap

- [ ] **CI/CD pipeline** – Set up continuous integration (automated testing, linting) and continuous deployment for the project.
- [ ] **Containerization** – Provide Docker or Conda environments for easier deployment and reproducibility of hap.py in different systems.
- [ ] **Documentation** – Expand user documentation and tutorials (e.g. README updates, example usage guides) once the codebase changes stabilize.

## Development Environment Setup

- **Environment management:** micromamba (recommended) or conda/mamba, dependencies in `environment.yml`
- **Bioinformatics tools:** Ensure `bcftools` and `samtools` are installed for certain tests. `bgzip`/`tabix` are optional as the code now defaults to the `pysam` Python implementation.

### Initial Setup

1. **Clone the repository and navigate to it:**
   ```bash
   git clone <repository-url>
   cd hap.py
   ```

2. **Create and activate the development environment:**
   - Using **micromamba** (recommended):
     ```bash
     micromamba create -n happy-dev python=3.11
     micromamba activate happy-dev
     ```
   - Using **venv** (alternative):
     ```bash
     python3 -m venv .venv
     source .venv/bin/activate  # (Windows: .venv\Scripts\activate)
     ```

3. **Install project in development mode with dependencies:**
   ```bash
   pip install -e ".[dev]"
   # Install pre-commit hooks for code quality
   pre-commit install
   ```

### Code Quality Tools

The project uses several tools to maintain code quality and style:

- **Black** – code formatter (with 88-character line limit)
- **Ruff** – linter for Python (fast, includes flake8/pyflakes checks)
- **isort** – import statement sorter (configurations compatible with Black)
- **mypy** – static type checker for Python
- **pre-commit** – framework for running linters/formatters on each commit

### Running Code Quality Checks

Use the following commands to format, lint, and type-check the code:
```bash
# Format code
black src/ tests/

# Lint and auto-fix issues
ruff check src/ tests/ --fix

# Sort imports
isort src/ tests/

# Static type checking (optional, if mypy is configured)
mypy src/hap_py/

# Run all pre-commit checks on the entire codebase
pre-commit run --all-files
```

## Testing

### Test Structure

- **Unit tests:** `tests/unit/` cover individual modules and functions (fast, isolated tests).
- **Integration tests:** `tests/integration/` cover end-to-end scenarios and require external tools or data.
- **Shared test utilities:** `tests/utils.py` provides common helper functions for tests, with configuration in `conftest.py`.

### Running Tests

**Note:** Always activate the `happy-dev` environment before running tests:

```bash
micromamba activate happy-dev
```

- **Run all unit tests:**
  ```bash
  pytest tests/unit/ -v
  ```
- **Run a specific unit test file:**
  ```bash
  pytest tests/unit/test_vcfeval.py -v
  ```
- **Run with coverage:**
  ```bash
  pytest tests/unit/ --cov=hap_py --cov-report=html
  ```
- **Run all integration tests** (requires external tools like RTG and access to example data):
  ```bash
  pytest tests/integration/ -v
  ```
- **Run a specific integration test:**
  ```bash
  pytest tests/integration/test_performance.py -v
  ```
- **Skip slow tests** (skip tests marked as slow):
  ```bash
  pytest tests/integration/ -v -m "not slow"
  ```
- **Run full test suite (all tests):**
  ```bash
  pytest tests/ -v
  ```
  (Use `-n auto` with `pytest-xdist` to run tests in parallel, if installed.)

### Test Markers
- `@pytest.mark.integration` – marks tests that require external tools or large data.
- `@pytest.mark.slow` – marks long-running tests.

### Common Test Issues (from Modernization)
During the transition from Python 2 to Python 3 and restructuring of the project, a few common issues were addressed:
1. **Import paths** – The package name changed. For example:
   ```python
   # Old import style (pre-modernization)
   import haplo.vcfeval

   # New import style (post-modernization)
   import hap_py.haplo.vcfeval
   ```
   Ensure tests use the updated `hap_py` package imports. Also verify that each package directory contains an `__init__.py` (so Python recognizes the package).
2. **Module not found errors** – Double-check that the `src/hap_py/` directory is on the Python path during testing. Installing in dev mode (`pip install -e .`) or using `pytest` from the repo root helps set this up. Missing `__init__.py` files in `tests/` subdirectories can also cause import errors (make sure `tests/` and its subfolders have `__init__.py`).
3. **String vs. bytes** – Use utility functions to handle byte strings vs Unicode strings. For instance, the project provides `ensure_str()` and `ensure_bytes()` in `hap_py.haplo.string_handling`. Use these when reading outputs from subprocesses to avoid type mismatches between Python 3 (which uses Unicode `str`) and older code expecting bytes.
4. **File path differences** – Use `pathlib.Path` for file paths to ensure cross-platform compatibility:
   ```python
   from pathlib import Path

   # Old way:
   test_file = os.path.join(os.path.dirname(__file__), "data", "test.vcf")

   # New way:
   test_file = Path(__file__).parent / "data" / "test.vcf"
   ```
   This makes path manipulations clearer and OS-agnostic.
5. **External tool availability** – Tests relying on external tools (e.g., `bcftools`, `rtg`) should check for tool presence and skip if not available:
   ```python
   import shutil, pytest
   if not shutil.which("bcftools"):
       pytest.skip("bcftools not available")
   ```
   Also ensure the RTG tools are built (expected at `build/external/rtg-tools/rtg`). If not present, tests should be skipped or the build instructions should be followed.

### Debugging Test Failures
When a test fails, consider the following steps to diagnose the issue:
1. **Check Python path and installation** – Confirm you are running tests in the correct environment and that `hap_py` is installed. For example:
   ```python
   import sys, hap_py
   print("Python path:", sys.path)
   print("hap_py location:", hap_py.__file__)
   ```
   Running `pip list | grep hap-py` can also verify that the package is installed in the environment.
2. **Verify external dependencies** – Ensure that required external tools and data are available:
   ```bash
   # List expected build outputs
   ls build/external/

   # Check that external tools can run
   bcftools --version
   samtools --version
   ```
3. **Run tests with verbose output** – Use `-s` (do not capture output) and `--tb=long` (full traceback) to get more insight into test failures:
   ```bash
   pytest tests/integration/test_some_failure.py -v -s --tb=long
   ```
   This can reveal detailed error messages from subprocesses or assertion failures.

## Build and Installation
### Versioning (from `pyproject.toml`)
The project uses **setuptools_scm** to manage the version number from Git tags. The `pyproject.toml` is configured so that no hard-coded version is needed; instead, the version is derived from the latest Git tag.

**Key configurations in `pyproject.toml`:**
```toml
[project]
name = "hap-py"
dynamic = ["version"]
# ...

[tool.setuptools_scm]
write_to = "src/hap_py/_version.py"
fallback_version = "0.4.0"
```
The version is written to `src/hap_py/_version.py` on build or install.

**Development install version:** If you do a development install (`pip install -e .`), and no git tag is present, it will fall back to `"0.4.0"` (as specified above) until a tag is added.

To properly set a version for release:
1. Tag the repository (e.g., `git tag v0.4.0` and push the tag).
2. Reinstall or build, and setuptools_scm will pick up the tag as the version.

### Development Installation
For an editable install (for development):
```bash
pip install -e .[dev]
```
This installs hap.py in development mode along with all development dependencies (`[dev]` extra in `pyproject.toml`).

### Production Build
To build a distributable package (wheel and source distribution):
```bash
python -m build  # requires build tool (pip install build)
# After building, install to verify:
pip install dist/hap.py-*.whl
```

## Troubleshooting Common Issues
### Import Errors
1. **Virtual environment not activated** – Always ensure you have `happy-dev` activated when working. If you encounter module import errors, double-check that you're using the intended Python environment.
2. **Package not installed** – If `import hap_py` fails, make sure you ran `pip install -e .` in your environment. Without the package installed (or the `src/` directory added to `PYTHONPATH`), tests will not find the modules.
3. **Missing `__init__.py`** – If tests cannot find modules under `tests/` (e.g., `ModuleNotFoundError: No module named 'tests.utils'`), add missing `__init__.py` files to those directories to make them part of a package.
4. **Conflicting module names** – Ensure there's no naming conflict (e.g., a script named `hap.py` in the working directory can shadow the package). Running tests from the repository root (so that `hap_py` package is found first) can help.

### External Tool Errors
1. **Build artifacts missing** – If tests fail due to missing tools (e.g., `rtg` not found), ensure external dependencies have been built. Check `build/external/` for expected directories (like `rtg-tools`).
2. **PATH not configured** – The `rtg` tool (and others like bcftools) need to be in your PATH for certain tests. You may need to update `PATH` to include the `external/` tools or specify their location via environment variables or test arguments.
3. **Tool installation** – Confirm that external bioinformatics tools are installed and accessible. If a tool is not installed, either install it or skip the tests requiring it (pytest will skip tests marked accordingly if the tool is missing, as configured).

### Test Data Issues
1. **Missing example data** – Ensure the `example/` directory with test datasets is present. Some integration tests expect data files in known locations; if they've been moved or not downloaded, tests will fail.
2. **File format or indexing** – Make sure reference files (like `.fa` reference genomes) are indexed properly (e.g., have `.fai` index if needed). If not, indexing tools should be run (e.g., `samtools faidx` for FASTA).
3. **File permissions** – Sometimes tests fail because of read/write permission issues, especially when writing temp files or logs. Ensure you have permission to write to the repository directories or use temp directories.

### Build Issues
1. **CMake version** – If CMake scripts fail, ensure your CMake is up-to-date (>= 3.10 as required). Older versions might not understand some CMake commands in the project.
2. **Compiler support** – Use a modern C++ compiler that supports C++11 or later, as needed. If the project uses any advanced C++17 features (unlikely given most code is Python now), your compiler version should support it.
3. **External dependency build** – If the `cmake --build` step fails, inspect the output. You might be missing a library (e.g., zlib for htslib) or have an incompatible toolchain. Installing required system packages or adjusting CMake options may be necessary.
4. **Clean rebuild** – If you encounter strange build errors, try deleting the `build/` directory and building from scratch:
   ```bash
   rm -rf build
   cmake -B build -S .
   cmake --build build
   ```

## Workflow & Development Guidelines
### Environment Management
- **Always activate the development environment** (`micromamba activate happy-dev`) at the start of each session. This ensures the correct Python version and all dependencies are in use. You can verify by checking `which python` (it should point to your `happy-dev` env path) and `python --version` (should show the expected Python 3.11.x).
- All commands for testing, building, etc., assume the `happy-dev` environment is active. If something isn't working, double-check that you're in the right environment.

### Development Workflow Tips
- **Run unit tests first** – Before tackling integration tests, ensure all unit tests pass. This catches low-level issues early.
- **Iterate on failing tests** – When debugging, start with the first failing test. Use `pytest -v -s` on that specific test to see detailed output and tracebacks. Once resolved, move to the next failure. This methodical approach prevents confusion.
- **Clean up after tests** – Remove or reset any temporary files or environment changes created during debugging. This avoids side effects that could cause subsequent tests to fail or pass incorrectly.
- **Ensure consistent paths** – Tests and code should use consistent paths. Use relative paths (or configurable paths) for test data so tests run on any system. We standardized many paths (using `pathlib.Path`) and environment variables for this reason.
- **External tools in tests** – If a test complains that `rtg` (or another tool) is not found, make sure:
  - You've built the external tools (the `rtg` binary should be in `build/external/rtg-tools/`).
  - The test is supplying the `--engine-vcfeval-path` argument to `hap.py` to point to the `rtg` binary. If not, update the test or the test configuration (e.g., fixture in `conftest.py`) to include the correct path.
  - As a quick workaround for development, you can create a symbolic link to `rtg` in the project root and add it to PATH:
    ```bash
    ln -s build/external/rtg-tools/rtg ./rtg  # link rtg to current directory
    export PATH="$PATH:$(pwd)"
    ```
    This ensures that when tests call `rtg`, it can be found. (Long-term, tests should handle this via configurations.)
- **Be mindful of test logic vs code issues** – If a test is failing, discern whether the test might be expecting outdated behavior (and thus needs updating) or if it revealed a bug in the code that needs fixing. We encountered both scenarios during modernization.
- **Keep commits clean** – When committing changes:
  - Remove any debug prints or temporary files used during development.
  - Restore any files that were accidentally deleted or renamed incorrectly during coding sessions.
  - Summarize the changes clearly in the commit message (e.g., "Fix off-by-one error in normalize_variant" or "Improve RTG tool detection in tests").
- **Reference original project for context** – If unsure about a behavior, consult the original hap.py repository ([Illumina/hap.py](https://github.com/Illumina/hap.py)). However, remember our fork has diverged (with C++ mostly replaced by Python), so not everything will directly apply.
- **Use consistent style** – Follow the established code style (PEP8, Black formatting) and structure in new code. This makes it easier for both humans and AI assistants to navigate the project.

## Debugging Integration Tests
Integration tests can be complex due to external dependencies and large data. Use this structured approach to diagnose and fix issues:
1. **Prepare the environment** – Ensure `happy-dev` is activated and all needed tools are built. For instance, verify the RTG binary exists and is executable at the expected location (`build/external/rtg-tools/rtg`). If not, rebuild or adjust the path.
2. **Identify failing tests** – Run the integration tests suite and note which tests fail. The output (or a log file if using `tee`, e.g. `pytest tests/integration/ -v | tee integration_test_output.txt`) will list failures. Focus on one test at a time.
3. **Examine error messages** – Open the detailed output for a failing test. Look at Python tracebacks and any tool output. Common causes include:
   - *Tool not found errors*: e.g., `"rtg: command not found"`. This indicates the test couldn't find the RTG binary. Solution: ensure `--engine-vcfeval-path` is provided and pointing correctly, or that `rtg` is in the PATH during the test.
   - *Assertion failures in outputs*: e.g., differences in expected vs actual output files. Check if the code produced all expected output (if an output like `roc.tsv` is missing, the code might not have created it due to an earlier error). Also verify that tests are not overly strict (additional output files should be ignored if they don't affect correctness).
   - *Hanging tests*: If a test stalls, it might be waiting on a subprocess or stuck writing to a temp file. Check if there's a prompt or if an external tool is waiting for input. Also consider if a temp directory issue is causing a deadlock (for example, trying to write large data to a full or unwritable location).
   - *VCF header or data issues*: For example, a test checking VCF headers might fail if the code doesn't include an expected FILTER field. Determine if the test expectation is correct or if the code needs to be adjusted (we fixed `_check_header` to handle optional FILTER fields).
4. **Apply fixes iteratively** – For each issue identified:
   - Update the test or code accordingly. If the test was missing a parameter (like the RTG path), add it. If the code had a bug (like off-by-one in normalization), fix the code.
   - Add any necessary cleanup in tests to remove temporary files (use `pytest` fixtures like `tmp_path` where possible to manage temp dirs automatically).
   - Rerun that test (or the subset of tests) to see if the issue is resolved.
5. **Re-run full test suite** – After fixing known issues, run the entire test suite (`pytest tests/ -v`) again. This ensures your fixes didn't break something else and that there are no new failures.
6. **Document and commit** – Once tests are all passing, make a commit detailing the fixes. Include references to test names or issue descriptions so others (or future you) can understand what was addressed.

## Recent Fixes and Notes (May 2025)
The following notable changes and fixes have been applied during the latest development cycle:
- **RTG template directory creation** – Fixed an issue with the `rtg format` command failing when a temporary SDF directory already existed. We replaced the use of `tempfile.NamedTemporaryFile` with `tempfile.mkdtemp` in `hap_py.haplo.vcfeval.runVCFEval`, ensuring a unique directory is created for RTG and preventing collisions.
- **RTG path detection** – Improved how the code locates the RTG toolkit. The function `findVCFEval()` in `vcfeval.py` now checks the project's bundled `rtg` path (in `build/external/rtg-tools`) in addition to checking the system PATH. This prevents false "executable not found" warnings when RTG is actually available in the expected location.
- **Suppressing false warnings** – Updated the package initialization (`hap_py/__init__.py`) to only warn about missing RTG tools if neither the PATH nor the bundled location has the executables. This removed redundant warnings when running tests.
- **Test package structure** – Added missing `__init__.py` files in the `tests/` directories (such as `tests/` and `tests/integration/`). This resolved import errors like `ModuleNotFoundError: No module named 'tests.utils'` by properly treating test directories as Python packages.
- **Binary output verification** – Confirmed that all expected binary scripts are present in `build/bin/` after building. This includes `hap.py`, `multimerge`, and the `qfy.py` wrapper.
- **Variant normalization fix** – Fixed an off-by-one error in the `normalize_variant` function (a test was expecting position 101 but the code returned 102). The logic was adjusted so that variant normalization now matches the expected behavior in tests.
- **VCF header check** – Modified the `_check_header` method in the VCF comparison module to handle scenarios where the VCF `FILTER` field may be missing. Tests expecting a strict check on FILTER were updated to either include the field or the code was made more flexible, in line with real-world VCFs.
- **Integration test enhancements** – Updated integration tests to consistently pass the `--engine-vcfeval-path` when using the RTG engine, ensuring the hap.py CLI knows where to find the `rtg` binary. The `get_rtg_path()` helper in `conftest.py` was also improved: it now respects an `RTG_PATH` environment variable and uses `shutil.which` to locate the RTG executable if not explicitly set, providing more robust test configuration.
- **Mocking and patching in tests** – Resolved an issue in a test (`test_runVCFEval_missing_output`) by correcting the order of `@patch` decorators. We noted that when multiple patches are applied, the order of arguments in the test function is the reverse of the patch application order. This detail was important to get the test working correctly and has been documented to avoid confusion in the future.

- **Deprecated tools removed** – Legacy binaries `hapcmp` and `hapenum` along with their tests have been removed. Documentation now notes their exclusion.
- **GA4GH annotation support** – Quantification writes GA4GH fields (`BD`, `BK`, `BI`, `BVT`, `BLT`, `QQ`) to output VCFs using a Python implementation.
- **Header comparisons** – Integration tests ignore VCF header lines when comparing outputs, matching the behavior of the original shell scripts.
- **Default compression with pysam** – Compression and indexing of VCF files use `pysam` by default; `bgzip` and `tabix` are optional fallbacks.

## Unit Test Coverage Plan

The goal is to rely on unit tests for most functionality. Add focused tests for
new Python modules and edge cases. Integration tests should eventually be used
only for a lightweight verification of the `hap.py` command line interface.

To increase unit test coverage:

1. Port complex logic from integration tests into unit tests. Target modules
   like `python_preprocess`, `vcfeval`, and CLI wrappers.
2. Use fixtures in `tests/conftest.py` to supply small reference files and VCF
   examples rather than large example data.
3. Mock external dependencies (pysam, subprocess) to isolate code paths.
4. Keep a minimal set of integration tests for the overall command line tools
   as sanity checks.
5. Add unit tests for helper modules in `hap_py.tools` such as `bcftools`
   and `bedintervaltree` to exercise parsing and interval logic.
6. Incrementally convert integration test logic to smaller unit tests that
   focus on specific functions. Use the integration suite only to verify
   the `hap.py` CLI behaves correctly end-to-end.
