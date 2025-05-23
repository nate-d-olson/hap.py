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
The project has achieved its primary modernization goals, transforming from Python 2/C++ hybrid to modern Python 3:
- **75% Complete** (6/8 components migrated) - **All critical functionality modernized**
- **C++ Binary Elimination**: **100% complete** - All critical C++ binaries replaced with Python
- **Production Ready**: Core workflows fully functional with pure Python implementations
- **Completed**: blocksplit, quantify, vcfcheck, preprocess, hapcmp, gvcf2bed
- **Remaining**: xcmp, scmp, multimerge (specialized algorithms, optional for most use cases)
- **RTG Tools**: Successfully integrated RTG Tools 3.12.1 for vcfeval functionality
- **pysam Compatibility**: Resolved all major compatibility issues with robust error handling

### Modern Python Architecture
The codebase has been successfully modernized to a primarily Python-based architecture:
- **Python implementations** (`src/hap_py/haplo/python_*.py`) - Modern implementations using pysam
- **Core functionality** - All critical processing now pure Python using pysam, pandas, numpy
- **External tools** - RTG Tools integration for vcfeval functionality (`libexec/rtg-tools-install/`)
- **Legacy C++ components** (`src/c++/`) - Remaining for specialized algorithms (xcmp, scmp, multimerge)
- **Cython extensions** (`src/hap_py/haplo/cython/`) - Optional performance bridges

### Key Dependencies
- **pysam** - VCF/BCF parsing (replacing custom C++ htslib usage)
- **pandas** - Data manipulation and metrics
- **numpy/scipy** - Statistical calculations
- **BioPython** - Sequence manipulation
- **RTG Tools** - Java-based vcfeval functionality (included installation)

## Major Modernization Achievements

### C++ Binary Elimination Complete
All critical C++ binaries have been successfully replaced:

- **vcfcheck**: Now uses direct VCFChecker class calls instead of external binary subprocess
- **preprocess**: Now uses direct PreprocessEngine class calls instead of external binary subprocess  
- **gvcf2bed**: Complete Python implementation for VCF→BED confident region extraction using pysam
- **RTG Tools**: Self-contained installation with wrapper scripts for vcfeval functionality

### Benefits Achieved
- **Installation**: No C++ compilation required for core functionality - dramatic simplification
- **Dependencies**: Simplified to pure Python + included Java tools (RTG)
- **Maintainability**: Modern Python code with type hints and standard practices
- **Testing**: Significantly improved integration test success rate
- **Cross-platform**: Better compatibility without complex build requirements
- **pysam Integration**: Robust VCF/BCF processing with type-aware error handling
- **Performance**: Maintained functionality while improving code clarity and debugging

## Current Status & Next Steps

### ✅ Modernization Complete for Core Use Cases
The package is now **production-ready** for all standard haplotype comparison workflows:
- VCF preprocessing and normalization
- Haplotype comparison and analysis  
- RTG vcfeval integration for benchmarking
- Quantification and metrics generation

### 🔧 Optional Future Enhancements
Remaining work items are optimizations rather than blockers:
1. **VCF Output Sorting**: Minor optimization for preprocessing output ordering
2. **Component Migration**: xcmp, scmp, multimerge (specialized algorithms)
3. **Performance Tuning**: Benchmarking and optimization opportunities
4. **Integration Tests**: Complete test suite validation

### 💡 Known Issues & Workarounds
- **VCF Sorting**: Preprocessing may produce unsorted output in edge cases (use bcftools sort if needed)
- **Large Files**: Memory usage optimization opportunities for very large VCF files
- **Legacy Components**: xcmp, scmp, multimerge still require C++ compilation if needed

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