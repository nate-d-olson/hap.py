# Migration Guide from Original hap.py

This guide helps users migrate from the original Illumina hap.py to this modernized version.

## Overview

This modernized version of hap.py maintains the core functionality while providing:

* **Python 3.8+ compatibility** (original required Python 2.7)
* **Modern packaging** with `pip` installation
* **Enhanced quantify module** with GA4GH compliance and statistical confidence intervals
* **Improved error handling** and user experience
* **Active maintenance** and development

## Breaking Changes

### Removed Components

Some components from the original hap.py were not included in this modernized version:

* **som.py**: The somatic variant comparison tool was not migrated to Python 3
* **scmp engine**: The specialized comparison engines (scmp-distancebased and scmp-somatic) were not included

For detailed information about these removed components, alternatives, and how to request their restoration, please see the [removed components documentation](removed_components.md).

### Python Version Requirements

**Original hap.py:**
```bash
python2.7 /path/to/hap.py/bin/hap.py --help
```

**Modernized hap.py:**
```bash
# After pip installation
hap.py --help

# Or using python module syntax
python -m hap_py.hap --help
```

### Installation Method

**Original hap.py:**
```bash
git clone https://github.com/Illumina/hap.py.git
cd hap.py
make install-tools
```

**Modernized hap.py:**
```bash
git clone https://github.com/nate-d-olson/hap.py.git
cd hap.py
pip install .
```

### Command Line Interface

The basic CLI remains the same, but entry points are now standardized:

**Original hap.py:**
```bash
/path/to/hap.py/bin/hap.py truth.vcf query.vcf -o output -r ref.fa
/path/to/hap.py/bin/som.py truth.vcf query.vcf -o output -r ref.fa  # Note: not available in modernized version
/path/to/hap.py/bin/qfy.py truth.vcf query.vcf -o output -r ref.fa
```

**Modernized hap.py:**
```bash
hap.py truth.vcf query.vcf -o output -r ref.fa
# som.py is not available in this version (see "Removed Components" section below)
qfy.py truth.vcf query.vcf -o output -r ref.fa
```

### Package Structure

**Original structure:**
```
hap.py/
├── src/python/Haplo/
├── src/c++/
└── bin/
```

**Modernized structure:**
```
hap.py/
├── src/hap_py/
│   ├── haplo/          # Core functionality
│   ├── quantify/       # Enhanced quantify module
│   ├── external/       # RTG management
│   └── tools/          # Utility tools
└── tests/              # Comprehensive test suite
```

## Compatible Features

### Core Functionality

All core benchmarking functionality remains compatible:

```bash
# Basic comparison (compatible)
hap.py truth.vcf query.vcf -f confident.bed -o benchmark -r reference.fa

# Using vcfeval engine (compatible)
hap.py truth.vcf query.vcf -o benchmark -r reference.fa --engine=vcfeval

# Preprocessing (compatible)
pre.py input.vcf -o normalized.vcf -r reference.fa

# Somatic comparison (compatible)
som.py truth.vcf query.vcf -o somatic_benchmark -r reference.fa
```

### Output Formats

All standard output formats are maintained:

* `.summary.csv` - High-level metrics summary
* `.metrics.json` - Detailed metrics in JSON format
* `.vcf.gz` - Annotated VCF files with benchmarking decisions
* `.bed` - Confident regions and FP regions

## Enhanced Features

### GA4GH Compliance

New in the modernized version:

```bash
# Enable GA4GH compliance mode
hap.py truth.vcf query.vcf -o output --quantify-method=ga4gh

# With stratification
hap.py truth.vcf query.vcf -o output --quantify-method=ga4gh \
    --ga4gh-stratification=stratification.bed
```

### Enhanced ROC Analysis

The modernized quantify module provides:

* **Bootstrap confidence intervals** for statistical rigor
* **Quality score stratification** across multiple thresholds
* **Multi-threshold analysis** at standard quality cutoffs
* **Comprehensive ROC curves** with precision-recall analysis

```bash
# ROC analysis is enabled by default
hap.py truth.vcf query.vcf -o benchmark -r reference.fa

# Output includes:
# benchmark.roc.tsv
# benchmark.quality_stratification.tsv
# benchmark.multi_threshold.tsv
```

## Troubleshooting Migration

### Environment Setup

If you encounter issues, ensure proper environment setup:

```bash
# Recommended: Use micromamba for environment management
micromamba create -n happy-dev python=3.11
micromamba activate happy-dev
pip install -e .[dev]
```

### RTG Tools

RTG tools are automatically managed in the modernized version:

```bash
# Check RTG availability
hap.py --list-engines

# RTG is automatically downloaded and configured
# No manual RTG installation required
```

### Common Issues

1. **"Command not found" errors:**
   ```bash
   # Ensure pip installation completed successfully
   pip list | grep hap
   which hap.py
   ```

2. **Python import errors:**
   ```bash
   # Verify package installation
   python -c "import hap_py; print(hap_py.__version__)"
   ```

3. **RTG-related errors:**
   ```bash
   # Check RTG path detection
   python -c "from hap_py.external.rtg_manager import get_rtg_path; print(get_rtg_path())"
   ```

## Development Migration

### For Developers

**Original development setup:**
```bash
export PYTHONPATH=/path/to/hap.py/src/python:$PYTHONPATH
```

**Modernized development setup:**
```bash
micromamba activate happy-dev
pip install -e .[dev]
pre-commit install
```

### Testing

**Original testing:**
```bash
cd /path/to/hap.py
bash src/sh/test_haplopy.sh
```

**Modernized testing:**
```bash
# Unit tests
pytest tests/unit/ -v

# Integration tests
pytest tests/integration/ -v

# All tests
pytest tests/ -v
```

### Code Quality

The modernized version includes automated code quality tools:

```bash
# Format code
black src/ tests/

# Lint code
ruff check src/ tests/ --fix

# Type checking
mypy src/hap_py/

# Pre-commit hooks (automatic on commit)
pre-commit run --all-files
```

## Performance Considerations

### Memory Usage

The modernized version maintains similar memory characteristics:

* **Whole genome:** 20-64GB RAM typical
* **Whole exome:** 8GB RAM sufficient
* **Small regions:** Desktop system compatible

### Processing Speed

* **Core algorithms:** Performance maintained or improved
* **RTG integration:** Optimized for faster startup
* **Python 3:** Modern optimizations provide better performance

## Support and Documentation

### Documentation

* **Main README:** [README.md](../README.md) - Updated installation and usage
* **Quantify Module:** [doc/quantify.md](quantify.md) - Enhanced features documentation
* **GA4GH Compliance:** [doc/ga4gh_compliance.md](ga4gh_compliance.md) - Standards compliance
* **ROC Analysis:** [doc/roc_analysis_guide.md](roc_analysis_guide.md) - Statistical analysis guide

### Getting Help

1. **Check the documentation** in the `doc/` directory
2. **Run tests** to verify installation: `pytest tests/unit/ -v`
3. **Report issues** on the project repository
4. **Review logs** for detailed error information

## Conclusion

The modernized hap.py maintains full compatibility for core use cases while providing significant improvements in:

* **Ease of installation** with pip
* **Development experience** with modern Python tools
* **Statistical rigor** with enhanced quantify module
* **Standards compliance** with GA4GH support
* **Long-term maintenance** with Python 3 support

For most users, migration involves simply changing the installation method and using the new command-line entry points. Advanced features like GA4GH compliance and enhanced ROC analysis are opt-in additions that don't affect existing workflows.
