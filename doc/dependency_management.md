# Dependency Management Updates

## Overview

This document describes the updates made to the dependency management system in the hap.py codebase as part of the Python 3 migration effort.

## Changes Made

### pyproject.toml Updates

1. **Core Dependencies**
   - Updated all dependency versions to match happy.requirements.py3.txt
   - Added setuptools and wheel as runtime dependencies
   - Fixed Cython version to 0.29.21 for compatibility
   - Specified exact versions for linting tools (black, ruff, isort, mypy)

2. **Optional Dependencies**
   - Organized dependencies into logical groups:
     - `cpp`: C++/Cython extension dependencies
     - `dev`: Development and testing tools
     - `viz`: Visualization libraries
     - `docs`: Documentation generation

3. **Tool Configurations**
   - Added target Python versions for black (py37-py310)
   - Configured ruff with basic rule sets (E, F, I)
   - Added file-specific ignores for common patterns
   - Set line length to 88 characters for all tools

## Dependency Versions

All dependencies have been updated to ensure compatibility with Python 3.7+ and to match the versions in happy.requirements.py3.txt:

| Package | Min Version | Notes |
|---------|-------------|-------|
| numpy | 1.19.0 | Required for array handling |
| pandas | 1.0.5 | For data manipulation |
| pysam | 0.16.0 | For SAM/BAM file access |
| scipy | 1.5.0 | For scientific computing |
| bx-python | 0.8.9 | For interval operations |
| Cython | 0.29.21 | For C++ extensions |
| setuptools | 50.3.2 | For packaging |
| wheel | 0.35.1 | For wheel distributions |
| pytest | 6.0.1 | For modern testing |
| matplotlib | 3.3.0 | For visualizations |
| seaborn | 0.11.0 | For statistical visualizations |

## Migration Progress

The dependency management system has been modernized to follow current Python packaging practices:

- Using pyproject.toml as the primary configuration file
- Proper dependency specification with minimum versions
- Optional dependencies organized by functionality
- Clear specification of Python version requirements
- Standardized code style tools with consistent configuration

## Future Work

1. **Type Annotation Improvements**
   - Fix remaining mypy type errors
   - Gradually enable stricter type checking options
   - Add complete type annotations to core modules

2. **Test Coverage Improvements**
   - Ensure tests cover all dependency combinations
   - Validate functionality with and without optional dependencies

3. **Documentation**
   - Update installation documentation
   - Create development environment setup instructions

## Usage Instructions

For a standard installation with all features:

```bash
pip install .[cpp,dev,viz]
```

For minimal installation (core functionality only):

```bash
pip install .
```

For specific feature sets:

```bash
pip install .[cpp]       # With C++/Cython acceleration
pip install .[viz]       # With visualization support
pip install .[dev]       # With development tools
```

## Using Pre-commit Hooks

After installing the development dependencies, set up pre-commit:

```bash
pip install .[dev]
pre-commit install
```

Run pre-commit on all files:

```bash
pre-commit run --all-files
```

Or on specific files:

```bash
pre-commit run --files src/python/Haplo/string_handling.py
```
