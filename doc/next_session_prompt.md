# Next Session Prompt for hap.py Python 3 Migration

## Current Status

The Python 3 migration of hap.py has made significant progress. Key milestones achieved and active work areas are outlined below.

## Completed Tasks

### CLI Tools Update

- ✅ Updated main CLI scripts (`hap.py`, `qfy.py`, `pre.py`) to be Python 3 compatible
- ✅ Removed auxiliary scripts (`ovc.py`, `cnx.py`) that weren't part of core functionality
- ✅ Added proper error handling and return status codes
- ✅ Created mock implementations for missing modules
- ✅ Added comprehensive tests for CLI tools
- ✅ Fixed sys.exit() usage to propagate status codes properly
- ✅ Added proper type hints to main functions

### Codebase Cleanup

- ✅ Removed unused C++ modules (scmp directory and XCmpQuantify files)
- ✅ Removed unnecessary Tools module (bamstats.py)
- ✅ Removed somatic and xcmp engine code from hap.py
- ✅ Updated tests to work without removed components
- ✅ Simplified engine selection to only support vcfeval

### Dependency Management

- ✅ Updated package dependencies in pyproject.toml to match happy.requirements.py3.txt
- ✅ Ensured all dependencies are Python 3 compatible
- ✅ Configured tool options (black, ruff, isort, mypy) in pyproject.toml
- ✅ Set specific versions for development tools to ensure consistency
- ✅ Added proper dependency organization with optional extras

### Documentation Updates

- ✅ Created `cli_updates.md` describing CLI tool improvements
- ✅ Created `cli_updates_status.md` tracking update status
- ✅ Updated `python3_migration.md` with CLI upgrade information
- ✅ Enhanced `dependency_management.md` with detailed information

## Active Work Areas

### Type Annotation Improvements (Priority: High)

- 🔄 Fixed critical type errors in string_handling.py using proper typing
- 🔄 Addressing mypy type errors in core modules
- ⬜ Complete type annotations for Haplo.cython.mock_internal
- ⬜ Fix remaining string handling issues at Python/C++ boundaries
- ⬜ Add proper type annotations to key functions in Tools modules

### Code Linting and Style (Priority: Medium)

- ⬜ Fix ruff/flake8 linting errors:
  - E402 Module level import not at top of file
  - E741 Ambiguous variable names
  - F401 Unused imports
  - F821 Undefined names

### Module Reorganization (Priority: Medium)

- ⬜ Fix import ordering issues
- ⬜ Resolve circular import dependencies
- ⬜ Improve package structure

### Testing (Priority: High)

- ⬜ Expand test coverage
- ⬜ Fix failing tests
- ⬜ Add additional integration tests
- ⬜ Migrate from nose to pytest

## Current Issues to Address

### Type Issues

- Incorrect return types in string_handling.py (fixed)
- Type errors in mock_internal.py (need to fix None initialization issues)
- String format type compatibility issues in various modules
- Missing type annotations in key functions

### Import Issues

- Module level imports not at top of files
- Circular import dependencies

### Code Style Issues

- Ambiguous variable names (using single letter variables)
- Undefined names in exception handling (traceback)

### Missing Implementations

- Some mock implementations may need to be replaced with actual implementations

## Next Steps

1. Continue fixing type annotation issues in remaining modules:
   - Focus on Haplo.cython.mock_internal.py first
   - Then address Tools.vcfextract.py

2. Update string formatting in modules to handle bytes/str conversion properly:
   - Tools.vcfcallerinfo.py

3. Fix import ordering and dependency issues in Tools modules

4. Address linting errors with ruff/black

## Working Environment Setup

For development, set up the workspace:

```bash
# Clone the repository (if not already done)
git clone https://github.com/Illumina/hap.py.git
cd hap.py

# Install development dependencies
pip install -e .[dev,cpp]

# Set up pre-commit hooks
pre-commit install

# Run type checking
mypy src/python/

# Run linting
ruff check src/python/
black --check src/python/
```
