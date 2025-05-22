# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**hap.py** is a variant calling benchmarking toolkit that performs genotype-level haplotype comparison against gold standard truth datasets. The project uses graph-based representations of VCF alleles and compares haplotype sequences through alignment/exact matching.

## Build & Development Commands

### Installation & Setup
```bash
# Install with all dependencies
pip install -e .[cpp,dev]

# Basic install (Python-only components)
pip install -e .
```

### Testing
```bash
# Run all tests
pytest tests/

# Unit tests only (fast)
pytest tests/unit -v

# Integration tests (requires C++ components)
pytest tests/integration -v

# Skip C++ dependent tests during Python-only development
pytest -m "not cpp" -v

# Run specific test file
pytest tests/unit/test_python_blocksplit.py -v
```

### Code Quality
```bash
# Format code
black src/

# Lint and fix
ruff check src/ --fix

# Sort imports
isort src/

# Type checking
mypy src/

# All quality checks (via tox)
tox -e lint
```

### Modernization Development
```bash
# Track modernization progress
python scripts/track_progress.py

# Generate progress report
python scripts/generate_report.py

# Validate Python replacements against C++ originals
python scripts/validate_replacements.py
```

## Architecture Overview

### Core CLI Tools
- **hap.py** - Main haplotype comparison engine
- **pre.py** - VCF preprocessing and normalization
- **qfy.py** - Quantification and metrics generation

### Python 3 Modernization Status
The project is actively migrating from Python 2/C++ to modern Python 3:
- **62.5% Complete** (5/8 components migrated)
- **Completed**: blocksplit, quantify, vcfcheck, preprocess, hapcmp
- **Remaining**: xcmp, scmp, multimerge

### Hybrid Architecture Pattern
The codebase uses a hybrid C++/Python approach:
- **C++ components** (`src/c++/`) - Performance-critical algorithms
- **Python replacements** (`src/hap_py/haplo/python_*.py`) - Modern implementations using pysam
- **Cython extensions** (`src/hap_py/haplo/cython/`) - Performance bridges
- **Fallback mechanisms** - Graceful degradation when C++ unavailable

### Key Dependencies
- **pysam** - VCF/BCF parsing (replacing custom C++ htslib usage)
- **pandas** - Data manipulation and metrics
- **numpy/scipy** - Statistical calculations
- **BioPython** - Sequence manipulation

## Testing Architecture

### Test Organization
- `tests/unit/` - Component-specific unit tests
- `tests/integration/` - End-to-end workflow tests
- `tests/conftest.py` - Shared fixtures and utilities

### Test Markers
- `@pytest.mark.integration` - Slow, full workflow tests
- `@pytest.mark.cpp` - Tests requiring C++ components
- `@pytest.mark.slow` - Long-running tests

### Test Data
Test datasets located in `src/data/` and `example/` directories with expected outputs for validation.

## Modernization Development Guidelines

### When Adding Python Replacements
1. Create `python_<component>.py` in `src/hap_py/haplo/`
2. Use established packages (pysam, pandas) over custom implementations
3. Add comprehensive unit tests in `tests/unit/test_python_<component>.py`
4. Implement compatibility wrapper for gradual migration
5. Update progress tracking in modernization scripts

### Code Standards
- Python 3.8+ syntax and features
- Type hints required for new code
- 88-character line length (Black formatting)
- All new Python code must pass ruff linting

### Performance Considerations
- Maintain C++ fallbacks during transition period
- Benchmark critical paths when replacing C++ components
- Use Cython for performance-critical Python code when needed

## Key File Locations

### Entry Points
- `src/hap_py/hap.py` - Main comparison tool
- `src/hap_py/pre.py` - Preprocessing
- `src/hap_py/qfy.py` - Quantification

### Core Modules
- `src/hap_py/haplo/` - Comparison algorithms and Python replacements
- `src/hap_py/tools/` - Utility modules (vcfextract, bcftools interface, metrics)

### Build Configuration
- `pyproject.toml` - Modern build configuration
- `setup.py` - Legacy build (being phased out)
- `pytest.ini` - Test configuration

### Documentation
- `MODERNIZATION_PLAN.md` - Detailed migration strategy
- `doc/` - User documentation and examples