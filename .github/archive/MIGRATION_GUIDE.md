# Migration Guide: Moving to Modern Python Installation

This guide explains how to migrate from the legacy installation method using `install.py`
to the new method using standard Python packaging tools.

## Legacy Installation (Deprecated)

Previously, hap.py was installed using the `install.py` script:

```bash
python install.py /path/to/installation
```

This method is now deprecated and will be removed in version 1.0.0.

## New Installation Method

The recommended way to install hap.py is now using pip:

```bash
# Install from source
git clone https://github.com/Illumina/hap.py.git
cd hap.py
pip install .

# Or install from PyPI (when available)
pip install hap.py
```

## What Changed?

- **Standard Python Packaging**: We now use a standard Python package structure with `pyproject.toml`
- **Entry Points**: All tools are available as command-line entry points after installation
- **Dependencies**: Dependencies are managed through pyproject.toml
- **C++/Cython Integration**: C++ components are automatically built as part of the installation

## Command Line Usage

Instead of calling scripts with full paths:

```bash
# Legacy method
/path/to/install/bin/hap.py truth.vcf query.vcf -r ref.fa -o output_prefix
```

You can now use the installed entry points directly:

```bash
# New method
hap truth.vcf query.vcf -r ref.fa -o output_prefix
```

## Available Entry Points

The following command-line tools are available after installation:

- `hap`: Main tool for haplotype comparison (was hap.py)
- `pre`: Preprocessing tool for VCF files
- `qfy`: Quantification tool for variant counts
- `som`: Somatic comparison tool (was som.py)
- `ftx`: Feature extraction tool
- `cnx`: Conversion tool
- `ovc`: Overlap checking tool

## Common Issues

### System Dependencies

You still need to install system dependencies:

```bash
# For Ubuntu/Debian
sudo apt install build-essential cmake libz-dev libbz2-dev liblzma-dev libcurl4-openssl-dev

# For macOS
brew install cmake zlib bzip2 xz curl
```

### Virtual Environments

Using a virtual environment is recommended:

```bash
python -m venv venv
source venv/bin/activate
pip install .
```

### Development Installation

For development, use an editable installation:

```bash
pip install -e .[dev]
```

## Testing

Tests can now be run using pytest:

```bash
# Run unit tests
python -m pytest tests/unit

# Run integration tests
python -m pytest tests/integration
```

## Developer Information

### Running Tests

Tests have been migrated from shell scripts to pytest format:

```bash
# Run all tests
pytest

# Run only unit tests
pytest tests/unit

# Run only integration tests
pytest tests/integration

# Run tests with specific markers
pytest -m "not integration"  # Skip integration tests
pytest -m "not cpp"          # Skip tests that require C++ components
```

### Type Checking

We're adding type hints to improve code reliability. You can run type checking with:

```bash
# Install mypy
pip install mypy

# Run type checking
mypy src/python
```

### Pre-commit Hooks

The project uses pre-commit hooks to enforce code quality:

```bash
# Install pre-commit
pip install pre-commit

# Install hooks
pre-commit install

# Run hooks manually
pre-commit run --all-files
```
