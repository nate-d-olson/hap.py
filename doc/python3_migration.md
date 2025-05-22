# Python 3 Migration Guide

## Overview

This document provides information about the migration of hap.py from Python 2 to Python 3 and the replacement of C++ components with pure Python implementations.

## Python 3 Compatibility

hap.py has been updated to work with Python 3.8 and newer. The following changes were made to ensure compatibility:

- Updated string handling to use unicode strings
- Replaced Python 2-specific idioms with Python 3 compatible code
- Updated imports for libraries that have changed in Python 3
- Added type hints to improve code clarity and enable static type checking
- Updated division operators to use true division by default

## Modern Python Features

The codebase now takes advantage of modern Python features:

- Type hints throughout the codebase
- F-strings for string formatting
- Pathlib for file path handling
- Modern packaging with pyproject.toml (PEP 517/518)
- Pre-commit hooks for code quality (black, ruff, isort, pyupgrade)

## Python Replacements for C++ Components

To simplify the codebase and make it more maintainable, the following C++ components have been replaced with pure Python implementations:

1. **blocksplit**: Python implementation using pysam
   - Functionality: Splits VCF files into blocks for parallel processing
   - Implementation: `Haplo.python_blocksplit.BlockSplitter`

2. **vcfcheck**: Python implementation using pysam
   - Functionality: Validates VCF files and reports issues
   - Implementation: `Haplo.python_vcfcheck.VCFChecker`

3. **quantify**: Python implementation using pandas and pysam
   - Functionality: Quantifies variant calls against truth datasets
   - Implementation: `Haplo.python_quantify.QuantifyEngine`

4. **preprocess**: Python implementation using pysam
   - Functionality: Preprocesses VCF files (normalization, decomposition, left-shifting)
   - Implementation: `Haplo.python_preprocess.PreprocessEngine`

5. **Sequence utilities**: Python implementation using BioPython
   - Functionality: FASTA handling, sequence manipulation
   - Implementation: `Haplo.sequence_utils.FastaReader`, `Haplo.sequence_utils.SequenceUtils`

## Installation

### Using pip

The recommended way to install hap.py is via pip:

```bash
pip install hap.py
```

For development purposes, you can install from the source with development dependencies:

```bash
pip install -e .[dev]
```

### Dependencies

The following dependencies are required:

- Python 3.8+
- pysam
- pandas
- biopython
- numpy
- pytest (for testing)
- black, ruff, isort (for code formatting)

## Testing

The codebase includes comprehensive test coverage for all the Python implementations:

```bash
# Run all tests
pytest

# Run specific component tests
pytest tests/unit/test_python_blocksplit.py
pytest tests/unit/test_vcfcheck.py
pytest tests/unit/test_quantify.py
pytest tests/unit/test_python_preprocess.py
pytest tests/unit/test_sequence_utils.py
```

## Known Issues and Limitations

- The Python implementations might have different performance characteristics compared to the C++ versions
- Some advanced features might not be fully implemented in the Python versions yet
- Edge cases that were handled by the C++ code might behave differently

## Migration Timeline

- **Phase 1 (Completed)**: Basic Python 3 compatibility
- **Phase 2 (Completed)**: Modern packaging and code quality tools
- **Phase 3 (Completed)**: Python implementations for core C++ components
  - blocksplit, vcfcheck, sequence utilities, quantify, preprocess
- **Phase 4 (In Progress)**: Comprehensive testing and documentation
- **Phase 5 (Planned)**: Performance optimization and edge case handling

## Demo Script

A demonstration script is available to showcase the Python implementations:

```bash
python scripts/demo_modernized.py --help
```

This script can be used to demonstrate each component:

```bash
# BlockSplitter demo
python scripts/demo_modernized.py --blocksplit --vcf example.vcf

# VCFChecker demo
python scripts/demo_modernized.py --vcfcheck --vcf example.vcf --reference ref.fa

# Quantify demo
python scripts/demo_modernized.py --quantify --truth-vcf truth.vcf --query-vcf query.vcf --reference ref.fa

# Preprocess demo
python scripts/demo_modernized.py --preprocess --vcf example.vcf --reference ref.fa

# Sequence utilities demo
python scripts/demo_modernized.py --sequence --reference ref.fa
```

## Feedback and Contributions

We welcome feedback and contributions to improve the Python implementations. Please open issues or pull requests on the GitHub repository.
