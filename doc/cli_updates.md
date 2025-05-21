# Command-Line Interface Updates

This document summarizes the updates made to the hap.py command-line tools during the Python 3 migration.

## Overview

The hap.py project includes these core command-line utilities for variant calling evaluation:

1. `hap.py` - Main tool for comparing variant calls
2. `qfy.py` - Quantification of variant calling performance
3. `pre.py` - VCF file preprocessing

**Note**: The auxiliary scripts `ovc.py` and `cnx.py` have been removed as part of the Python 3 migration to focus on core functionality.

## Key Updates

### Entry Point Structure

The entry points for the command-line tools are defined in `pyproject.toml`:

```toml
[project.scripts]
hap = "hap.py:main"
qfy = "qfy:main"
pre = "pre:main"
```

Each script has been updated with:

Each script has been updated with:

1. **Proper Return Types**: All `main()` functions now return an integer status code (0 for success, non-zero for failure)
2. **Type Hints**: Added Python type hints to improve code reliability
3. **Error Handling**: Consistent error handling with proper exit codes
4. **Python 3 Path Handling**: Updated import paths to use Python 3 library locations
5. **Docstrings**: Added descriptive docstrings to main functions

### Python 3 Compatibility

- Updated string formatting to use f-strings
- Fixed path handling with the `pathlib` module where appropriate
- Updated error handling to use context managers
- Improved command-line argument processing

### Individual Script Updates

#### hap.py

- Added return type annotation to `main()` function
- Updated error handling to use `sys.exit()` with proper return codes
- Ensured consistent path handling
- Removed support for xcmp and somatic comparison engines
- Simplified to only support vcfeval as the comparison engine

#### qfy.py

- Added return type annotation to `main()` function
- Updated error handling for Python 3 compatibility
- Ensured proper exit code propagation

#### pre.py

- Added return type annotation to `main()` function
- Added proper error handling with try/except
- Updated exit code handling to follow Python standards

## Testing

The CLI tools have been tested using integration tests to verify:

1. Basic functionality works as expected
2. Error handling works correctly
3. Exit codes are properly propagated

## Future Improvements

1. **Error Messages**: Further refinement of error messages for better user experience
2. **Help Text**: Update help text to be more consistent across tools
3. **Type Hints**: Add more comprehensive type hints throughout the codebase
4. **Logging**: Standardize logging approach across all tools

## Usage Examples

Basic usage patterns remain unchanged from previous versions:

```bash
# Compare variant calls
hap.py truth.vcf query.vcf -r reference.fa -o output

# Quantify variant calls
qfy.py variants.vcf -r reference.fa -o output

# Preprocess a VCF file
pre.py input.vcf output.vcf -r reference.fa
```
