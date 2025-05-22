# Add Python implementation of preprocess component

## Overview

This PR adds a pure Python implementation of the preprocess component, which completes the migration of all high-priority C++ components to Python. The preprocess component is responsible for normalizing, decomposing, and left-shifting variants in VCF files.

## Changes

- Added `python_preprocess.py` with a pure Python implementation using pysam
- Added comprehensive test suite in `test_python_preprocess.py`
- Updated the modernization plan to reflect progress
- Added preprocess demo to `demo_modernized.py`
- Created a simplified standalone demo in `simple_preprocess.py`
- Updated documentation to include the new component
- Created MODERNIZATION_SUMMARY.md to track overall progress

## Features Implemented

The Python preprocess implementation includes the following features:

- Variant decomposition (splitting multi-allelic variants)
- Left-alignment of variants using reference genome
- Normalization (trimming common prefixes/suffixes)
- Special handling for haploid regions (e.g., chrX in males)
- Region filtering
- PASS-only filtering
- BCF output support

## Testing

The implementation has been tested with:

- Unit tests for all major functions
- Integration tests with real-world VCF files
- Validation against the C++ implementation
- Standalone testing script for quick verification

## Documentation

- Updated `MODERNIZATION_PLAN.md` to mark preprocess as completed
- Created `MODERNIZATION_SUMMARY.md` to track overall modernization progress
- Updated `scripts/README.md` to include the new component
- Added a detailed migration guide in `doc/python3_migration.md`
- Added comments and docstrings throughout the code

## Known Issues

- Some edge cases in complex variant normalization need additional testing
- Left-shifting of complex variants may behave differently from the C++ implementation in rare cases
- Performance optimizations are still needed for large VCF files

## Next Steps

With all high-priority C++ components now migrated to Python, the next steps are:

1. Complete type hint coverage for the codebase
2. Migrate remaining shell script tests to pytest
3. Update documentation for Python 3 compatibility
4. Create CI/CD pipeline for automated testing
5. Begin migration of medium-priority components (xcmp, scmp, hapcmp)
