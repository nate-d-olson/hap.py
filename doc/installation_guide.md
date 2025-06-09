# Installation and Setup Guide

This guide provides comprehensive instructions for installing and setting up the modernized hap.py package.

## Quick Start

For users who just want to get started quickly:

```bash
# Clone repository
git clone https://github.com/nate-d-olson/hap.py.git
cd hap.py

# Install
pip install .

# Verify installation
hap.py --help
```

## Table of Contents

1. [Prerequisites](#prerequisites)
2. [Installation Methods](#installation-methods)
3. [Development Setup](#development-setup)
4. [Verification](#verification)
5. [Environment Configuration](#environment-configuration)
6. [Troubleshooting](#troubleshooting)
7. [Uninstallation](#uninstallation)

## Prerequisites

### System Requirements

**Operating Systems:**
* Linux (Ubuntu 18.04+, CentOS 7+, RHEL 7+)
* macOS 10.14+ (Mojave or later)
* Windows 10+ (experimental, via WSL2 recommended)

**Hardware:**
* **Minimum:** 8GB RAM, 2 CPU cores, 10GB disk space
* **Recommended:** 32GB RAM, 8+ CPU cores, 50GB disk space
* **Whole genome analysis:** 64GB RAM, 16+ CPU cores, 100GB disk space

### Software Prerequisites

**Required:**
* Python 3.8+ (recommended: Python 3.11)
* pip (Python package manager)
* Git

**Bioinformatics Tools (recommended):**
* bcftools (for VCF processing)
* samtools (for BAM/CRAM processing)  
* tabix (for indexed VCF files)

**Build Tools (for source installation):**
* CMake 3.10+
* C++ compiler (GCC 7+ or Clang 10+)
* Make

### Installing Prerequisites

**Ubuntu/Debian:**
```bash
sudo apt update
sudo apt install -y python3 python3-pip git cmake build-essential
sudo apt install -y bcftools samtools tabix  # Optional bioinformatics tools
```

**CentOS/RHEL:**
```bash
sudo yum update -y
sudo yum install -y python3 python3-pip git cmake gcc-c++ make
# For bioinformatics tools, enable EPEL repository
sudo yum install -y epel-release
sudo yum install -y bcftools samtools tabix
```

**macOS:**
```bash
# Install Homebrew if not available
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"

# Install prerequisites
brew install python@3.11 git cmake
brew install bcftools samtools htslib  # Optional bioinformatics tools
```

**Windows (WSL2):**
```bash
# Enable WSL2 and install Ubuntu
# Then follow Ubuntu instructions above
```

## Installation Methods

### Method 1: Standard pip Installation (Recommended)

This is the recommended method for most users:

```bash
# Clone the repository
git clone https://github.com/nate-d-olson/hap.py.git
cd hap.py

# Install the package
pip install .

# Verify installation
hap.py --version
hap.py --help
```

### Method 2: Development Installation

For developers or users who want to modify the code:

```bash
# Clone the repository
git clone https://github.com/nate-d-olson/hap.py.git
cd hap.py

# Install in editable mode with development dependencies
pip install -e .[dev]

# Set up pre-commit hooks (optional but recommended)
pre-commit install

# Verify installation
hap.py --help
pytest tests/unit/ -v  # Run unit tests
```

### Method 3: Environment-Isolated Installation

Using micromamba (recommended for development):

```bash
# Install micromamba if not available
curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj bin/micromamba

# Create environment
micromamba create -n happy-dev python=3.11
micromamba activate happy-dev

# Clone and install
git clone https://github.com/nate-d-olson/hap.py.git
cd hap.py
pip install -e .[dev]

# Verify installation
hap.py --help
```

Using conda/mamba:

```bash
# Create environment
conda create -n happy-dev python=3.11
conda activate happy-dev

# Install dependencies
conda install -c bioconda bcftools samtools tabix

# Clone and install hap.py
git clone https://github.com/nate-d-olson/hap.py.git
cd hap.py
pip install -e .[dev]
```

### Method 4: Build from Source

For advanced users who need custom compilation:

```bash
# Clone repository
git clone https://github.com/nate-d-olson/hap.py.git
cd hap.py

# Install build dependencies
pip install build wheel

# Build package
python -m build

# Install built package
pip install dist/hap_py-*.whl
```

## Development Setup

### Full Development Environment

For contributors and advanced developers:

```bash
# Create development environment
micromamba create -n happy-dev python=3.11
micromamba activate happy-dev

# Clone repository
git clone https://github.com/nate-d-olson/hap.py.git
cd hap.py

# Install all development dependencies
pip install -e .[dev,viz,docs]

# Set up development tools
pre-commit install

# Verify setup
hap.py --help
pytest tests/unit/ -v
black --check src/
ruff check src/
```

### IDE Configuration

**Visual Studio Code:**
```json
// .vscode/settings.json
{
    "python.linting.enabled": true,
    "python.linting.ruffEnabled": true,
    "python.formatting.provider": "black",
    "python.testing.pytestEnabled": true,
    "python.testing.pytestArgs": ["tests/"]
}
```

**PyCharm:**
* Set Python interpreter to your environment
* Configure pytest as test runner
* Install Black and Ruff plugins

## Verification

### Basic Verification

```bash
# Check installation
hap.py --version
hap.py --help

# Test basic functionality
python -c "import hap_py; print('✓ Package imports successfully')"

# Check entry points
which hap.py
which som.py
which pre.py
which qfy.py
```

### RTG Tools Verification

```bash
# Check RTG tools integration
hap.py --list-engines

# Should show 'vcfeval' in the list
# RTG tools are automatically downloaded if needed
```

### Comprehensive Testing

```bash
# Run unit tests (fast)
pytest tests/unit/ -v

# Run integration tests (slower, requires external tools)
pytest tests/integration/ -v

# Run all tests
pytest tests/ -v

# Run tests with coverage
pytest tests/unit/ --cov=hap_py --cov-report=html
```

### Test Installation with Example Data

```bash
# Navigate to example directory
cd example

# Run a simple test
hap.py chr21.refcalls.vcf.gz hc.vcf.gz \
    -f performance.confident.bed.gz \
    -o test_run \
    -r chr21.fa

# Check outputs
ls test_run.*
# Should see: test_run.summary.csv, test_run.metrics.json, etc.
```

## Environment Configuration

### Environment Variables

**Optional but recommended:**

```bash
# Add to ~/.bashrc or ~/.zshrc

# Reference genome (for integration tests)
export HGREF="/path/to/reference/genome.fa"

# Temporary directory (optional)
export TMPDIR="/path/to/large/tmp/directory"

# RTG tools path (usually auto-detected)
export RTG_PATH="/path/to/rtg/tools"
```

### Configuration Files

Create `~/.happy_config.json` for persistent settings:

```json
{
    "default_engine": "vcfeval",
    "default_output_format": "ga4gh",
    "rtg_path": "/path/to/rtg/tools",
    "reference_genome": "/path/to/reference.fa",
    "temp_directory": "/path/to/tmp"
}
```

## Troubleshooting

### Common Installation Issues

**1. Permission Errors:**
```bash
# Use user installation
pip install --user .

# Or use virtual environment
python -m venv happy_env
source happy_env/bin/activate
pip install .
```

**2. Missing Dependencies:**
```bash
# Update pip and setuptools
pip install --upgrade pip setuptools wheel

# Install missing system packages
# Ubuntu/Debian:
sudo apt install python3-dev cmake build-essential

# CentOS/RHEL:
sudo yum install python3-devel cmake gcc-c++
```

**3. RTG Tools Issues:**
```bash
# Check RTG path
python -c "from hap_py.external.rtg_manager import get_rtg_path; print(get_rtg_path())"

# Manual RTG download (if automatic fails)
cd external/
wget https://github.com/RealTimeGenomics/rtg-tools/releases/download/3.12.1/rtg-tools-3.12.1-linux-x64.zip
unzip rtg-tools-3.12.1-linux-x64.zip
```

**4. Import Errors:**
```bash
# Check Python path
python -c "import sys; print('\\n'.join(sys.path))"

# Reinstall package
pip uninstall hap_py
pip install -e .
```

### Runtime Issues

**1. Memory Errors:**
```bash
# Monitor memory usage
htop

# Use smaller chromosomes for testing
hap.py truth.vcf query.vcf -o output -r ref.fa -l chr22
```

**2. Temporary File Issues:**
```bash
# Set temporary directory with more space
export TMPDIR=/path/to/large/tmp
mkdir -p $TMPDIR
```

**3. VCF Processing Errors:**
```bash
# Validate VCF files
bcftools view -h input.vcf
tabix -p vcf input.vcf.gz  # Create index if missing
```

### Getting Help

1. **Check logs:** Most tools provide detailed logging with `-v` or `--verbose`
2. **Run diagnostics:** Use `hap.py --help` to check available options
3. **Test environment:** Run unit tests to verify setup
4. **Review documentation:** Check relevant docs in `doc/` directory
5. **Report issues:** Open an issue on the project repository

## Uninstallation

### Remove Package

```bash
# Uninstall hap.py
pip uninstall hap_py

# Remove development environment (if using micromamba)
micromamba env remove -n happy-dev

# Remove configuration files (optional)
rm ~/.happy_config.json
```

### Clean Up

```bash
# Remove cloned repository (if no longer needed)
rm -rf /path/to/hap.py

# Remove downloaded RTG tools (optional)
rm -rf ~/.cache/hap_py/rtg-tools*
```

## Next Steps

After successful installation:

1. **Read the user guide:** [doc/happy.md](happy.md)
2. **Try examples:** Use data in `example/` directory
3. **Review migration guide:** [doc/migration_guide.md](migration_guide.md) (if coming from original hap.py)
4. **Explore advanced features:** [doc/quantify.md](quantify.md) and [doc/ga4gh_compliance.md](ga4gh_compliance.md)
5. **Set up development environment:** If you plan to contribute

## Support

For installation support:

* **Documentation:** Review all files in `doc/` directory
* **Test your setup:** Run the verification steps above
* **Check examples:** Try the provided example data
* **Report issues:** If problems persist, report on the project repository with:
  * Operating system and version
  * Python version (`python --version`)
  * Installation method used
  * Complete error messages
  * Steps to reproduce the issue
