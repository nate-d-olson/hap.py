# Next Session Prompt for hap.py Python 3 Migration

## Progress So Far

### CLI Tools Update (Completed)

- Updated main CLI scripts (`hap.py`, `qfy.py`, `pre.py`) to be Python 3 compatible
- Removed auxiliary scripts (`ovc.py`, `cnx.py`) that weren't part of core functionality
- Added proper error handling and return status codes
- Created mock implementations for missing modules
- Added comprehensive tests for CLI tools
- Created documentation for CLI updates
- Fixed sys.exit() usage to propagate status codes properly
- Added proper type hints

### Codebase Cleanup (Completed)

- Removed unused C++ modules (scmp directory and XCmpQuantify files)
- Removed unnecessary Tools module (bamstats.py)
- Removed somatic and xcmp engine code from hap.py
- Updated tests to work without removed components
- Simplified engine selection to only support vcfeval

### Documentation Updates

- Created `cli_updates.md` describing CLI tool improvements
- Created `cli_updates_status.md` tracking update status
- Updated `python3_migration.md` with CLI upgrade information

### Tests

- Created test suite for CLI tools
- Added integration tests
- Implemented test utilities for debugging

## Next Steps

### Dependency Management (Priority: High)

- Update package dependencies in pyproject.toml
- Ensure all dependencies are Python 3 compatible
- Address any dependency conflicts or version issues

### Type Annotation Improvements (Priority: High)

- Fix mypy type errors identified during the commit phase
- Add complete type annotations to core modules
- Resolve Union[type] imports and usage

### Code Linting and Style (Priority: Medium)

- Fix ruff/flake8 linting errors:
  - E402 Module level import not at top of file
  - E741 Ambiguous variable names
  - F401 Unused imports
  - F821 Undefined names

### Module Reorganization (Priority: Medium)

- Fix import ordering issues
- Resolve circular imports
- Improve package structure

### Testing (Priority: High)

- Expand test coverage
- Fix failing tests
- Create additional integration tests

### Documentation (Priority: Medium)

- Update user documentation
- Improve code documentation with Google-style docstrings
- Document known issues and workarounds

## Current Issues to Address

### Type Issues

- Multiple mypy errors in various modules
- Missing type annotations in key functions
- Incorrect type annotations

### Import Issues

- Module level imports not at top of files
- Circular import dependencies

### Code Style Issues

- Ambiguous variable names (using single letter variables)
- Undefined names in exception handling (traceback)

### Missing Implementations

- Some mock implementations may need to be replaced with actual implementations

## Next Milestone

The next milestone is to address the type annotation and linting issues to ensure the code passes basic quality checks, followed by expanding test coverage to validate functionality.
