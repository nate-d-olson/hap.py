# hap.py - Bioinformatics Variant Calling Benchmarking Tool

hap.py is a bioinformatics tool for benchmarking small variant calls. It analyzes and compares variant calls from different variant callers against a truth set to measure performance metrics like precision and recall. The tool is widely used in the genomics community despite its legacy codebase.

## Project Goal

We are modernizing this codebase by:

1. Converting from Python 2 to Python 3
2. Updating outdated dependencies
3. Implementing modern Python packaging standards
4. Improving installation process
5. Adding type hints and better documentation
6. Replacing ad-hoc testing with proper frameworks
7. Optimizing performance for genomic datasets

## Project Overview
### Tech Stack

- use `.codex/run_in_container.sh` for sandboxed code execution
- Python 3.6+ (migrating from Python 2)
- C++ (core algorithms and performance-critical components)
- CMake (build system)
- Cython (Python/C++ interface)
- Bcftools/htslib (for VCF/BAM handling)
- Boost C++ libraries (subset)
- Pandas, NumPy (data analysis)

### Key Considerations

- Bioinformatics data tends to be large - consider memory usage
- Maintain backward compatibility with existing workflows
- Preserve core algorithms while modernizing the code around them
- Focus on reliability and correctness - this is scientific software
- Make installation process more user-friendly

### Development Environment

- Use `install.py` to build the tool in a temporary directory
- Always run the full test suite after changes
- Implement proper version control practices
- Follow the coding standards in the instructions files


## Development Workflow Guidelines

### Build, install, and test

Run the installation script with a temporary build directory:
```bash
python install.py /tmp/happy-build
```

Verify the installation with the test suite:
```bash
cd /tmp/happy-build
src/sh/run_tests.sh
```

### Handling Build failures
If the build fails in the external C++ libraries:

1. Diagnose the issue:
```
python install.py /tmp/happy-build VERBOSE=1 2>&1 | tee build_log.txt
grep -A 10 "error:" build_log.txt
```

2. Incremental builds:
   * Attempt to build external dependencies only:
```
python install.py --build-externals-only /tmp/happy-ext
```

   * Check system dependencies are installed:

```
# For Ubuntu/Debian
sudo apt install build-essential cmake libz-dev libbz2-dev liblzma-dev libcurl4-openssl-dev

# For macOS
brew install cmake zlib bzip2 xz curl
```

3. Fix build issues following the guidance in build-system-fixes.instructions.md

## Development Workflow

### Making Code Changes

1. Testing workflow:

   * Ensure code is python3 compatible and follows best practices for formatting, documentation, type hints, and packaging
   * Test Python 3 updates with mock C++ interfaces when needed
   * Integrate both changes when each part is stable

2. Code quality workflow:

  * Install pre-commit hooks: `pip install pre-commit && pre-commit install`
  * Use pre-commit for automated checks: `pre-commit run --files <changed_files>`
  * Fix any issues reported by pre-commit
  * Run specific hooks as needed:

```bash
# Format specific files with Black
pre-commit run black --files src/python/path/to/file.py

# Run specific Python 3 compatibility checks
pre-commit run pyupgrade --files src/python/path/to/file.py

# Run linting with auto-fixes
pre-commit run ruff --files src/python/path/to/file.py
```

### Testing Process

1. Run build verification tests after any CMake changes
2. Run unit tests for components you modified
3. Run integration tests to ensure components work together
4. Capture and analyze test failures to identify root causes
5. Document any new test cases you add

### Code Review

1. Self-review your changes before submission
2. Address all comments from reviewers
3. Verify changes with appropriate tests
4. Update documentation to reflect changes

### Debugging

1. Use logging to trace execution flow
2. For build issues, use CMake's verbose output and message() commands
3. Test fixes thoroughly before committing
4. Document root causes of significant bugs for future reference

## Release Process

1. Update version numbers in relevant files
2. Update the RELEASES.md file with changes
3. Ensure all tests pass
4. Build and verify the release package
5. Create a tagged release in the repository

## Documentation

1. Update API documentation for any changed functions
2. Keep the main README up-to-date
3. Document any new features or significant changes
4. Ensure installation instructions remain accurate
