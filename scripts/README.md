# hap.py Modernization

This directory contains scripts for modernizing the hap.py codebase from Python 2 to Python 3 and replacing C++ components with pure Python implementations.

## Scripts Overview

- `run_modernization.py`: Main script to run all modernization steps
- `setup_project_config.py`: Set up modern Python project configuration
- `fix_python2_artifacts.py`: Fix common Python 2 artifacts
- `track_progress.py`: Track modernization progress
- `generate_report.py`: Generate a comprehensive progress report
- `generate_dashboard.py`: Create a visual dashboard of progress
- `validate_replacements.py`: Validate Python replacements for C++ components

## Python Replacements

The following C++ components have been replaced with Python implementations:

- `blocksplit`: Python implementation using pysam (see `python_blocksplit.py`)
- `vcfcheck`: Python implementation using pysam (see `python_vcfcheck.py`)
- `quantify`: Python implementation using pandas and pysam (see `python_quantify.py`)
- `preprocess`: Python implementation using pysam (see `python_preprocess.py`)
  - Provides variant decomposition, normalization, and left-shifting
  - Handles haploid regions and region filtering
  - Includes standalone version in `simple_preprocess.py`
- Sequence utilities: Python implementation using BioPython (see `sequence_utils.py`)

## Demonstration and Testing Scripts

- `demo_modernized.py`: Demonstrate the Python replacements for C++ components
- `simple_preprocess.py`: Standalone script for VCF preprocessing
- `replace_preprocess.py`: Helper script to replace C++ preprocess calls

## Getting Started

### Setup and Installation

1. Ensure you have Python 3.8+ installed
2. Install dependencies:

   ```bash
   pip install -e .[dev]
   ```

3. Install pre-commit hooks:

   ```bash
   pre-commit install
   ```

### Running the Modernization Process

To run all modernization steps:

```bash
python scripts/run_modernization.py --all
```

Or run specific steps:

```bash
# Just fix Python 2 artifacts
python scripts/run_modernization.py --fix-py2

# Just generate a progress report
python scripts/run_modernization.py --report

# Run pre-commit hooks on all files
python scripts/run_modernization.py --pre-commit
```

### Running Tests

To run tests for the Python replacements:

```bash
pytest tests/unit/test_python_blocksplit.py tests/unit/test_sequence_utils.py tests/unit/test_vcfcheck.py -v
```

## Progress Tracking

To track modernization progress:

```bash
python scripts/track_progress.py
```

To generate a visual dashboard:

```bash
python scripts/generate_dashboard.py --open
```

## Development Guidelines

1. Use Python 3.8+ syntax and features
2. Add type hints to all new code
3. Write comprehensive unit tests
4. Run pre-commit hooks before committing
5. Follow the modernization plan in `MODERNIZATION_PLAN.md`
6. Use pysam for VCF/BCF handling
7. Use BioPython for sequence manipulation
8. Validate Python replacements against C++ originals

## Contributing

1. Create a feature branch for your changes
2. Run pre-commit hooks on your changes
3. Add tests for new functionality
4. Ensure existing tests pass
5. Update documentation as needed
6. Submit a pull request

## Useful Resources

- [pysam Documentation](https://pysam.readthedocs.io/)
- [BioPython Documentation](https://biopython.org/docs/latest/api/)
- [Python Type Hints Documentation](https://docs.python.org/3/library/typing.html)
- [Pre-commit Documentation](https://pre-commit.com/)
- [Black Documentation](https://black.readthedocs.io/)
- [Ruff Documentation](https://beta.ruff.rs/docs/)
