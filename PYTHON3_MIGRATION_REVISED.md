# Revised Python 3 Migration Strategy for hap.py

## Core Functionality Focus

This revised strategy focuses on migrating the essential parts of hap.py to Python 3 while modernizing the codebase to improve maintainability and installation experience.

### Prioritized Components

1. **Core hap.py Functionality:**
   - vcfeval comparison engine integration
   - Stratified performance metrics using TSV input
   - Essential preprocessing capabilities
   - Core output generation

2. **Components to Defer/Remove:**
   - Somatic variant calling support (som.py)
   - bamstats.py script
   - scmp and xcmp comparison engines
   - Non-essential utility scripts

## Technical Approach

### Phase 1: Core Migration and Simplification (2-3 weeks)

1. **Create Feature Branches:**
   - Create `feature/somatic-support` branch to preserve somatic functionality
   - Create `feature/alternative-engines` branch for scmp and xcmp engines
   - Move non-core components to these branches before beginning migration

2. **Minimize External Dependencies:**
   - Replace custom C++ VCF processing with modern pysam (≥0.15)
   - Evaluate pybedtools for interval operations
   - Maintain only necessary C++ components for vcfeval integration

3. **Update Core Files:**
   - Focus on main hap.py script and vcfeval integration first
   - Update Tools/ module components needed for core functionality
   - Update Haplo/ module focusing on vcfeval.py and supporting files

4. **Optimize Preprocessing:**
   - Redesign preprocessing to use fewer temporary files
   - Implement smarter chunking based on variant density
   - Use modern multiprocessing approaches

### Phase 2: Modernization and Testing (2-3 weeks)

1. **API Modernization:**
   - Implement proper package structure
   - Add comprehensive type hints
   - Update docstrings to Google style
   - Use pathlib for file handling

2. **VCF Compatibility:**
   - Ensure compatibility with VCF spec v4.5
   - Create robust error handling for malformed VCFs
   - Add validation for VCF inputs

3. **Testing Infrastructure:**
   - Implement pytest-based test suite
   - Create integration tests for core functionality
   - Add CI pipeline for automated testing

### Phase 3: Packaging and Documentation (1-2 weeks)

1. **Modern Installation:**
   - Create proper pyproject.toml/setup.py
   - Implement pip-installable package
   - Provide containerized distribution option

2. **Documentation:**
   - Update user documentation focused on core workflows
   - Create clear error messages and troubleshooting guides
   - Document API for programmatic usage

## Implementation Details

### Dependency Management

1. **Modern Python Libraries:**
   - pysam ≥0.15 for VCF/BAM handling
   - pybedtools for interval operations
   - attrs/dataclasses for data structures
   - click for command-line interface

2. **Retained C++ Components:**
   - Necessary vcfeval integration points
   - High-performance evaluation algorithms (where Python would be too slow)

### Preprocessing Optimization

1. **Improved Data Chunking:**
   - Use variant density estimation for optimal chunk sizes
   - Minimize file I/O with better streaming approaches
   - Implement lazy evaluation for performance metrics

2. **Memory Efficiency:**
   - Use generators instead of materializing full lists
   - Implement proper resource management
   - Add memory profiling to detect leaks

### VCF Compatibility

1. **Robust Parsing:**
   - Add more permissive VCF parsing options
   - Implement auto-detection and fixes for common issues
   - Support both strict and lenient parsing modes

2. **Format Standardization:**
   - Normalize variant representations
   - Handle multiallelic variants consistently
   - Standardize genotype representation

## Migration Checklist

1. [ ] Create feature branches for deferred functionality
2. [ ] Audit core dependencies and identify modern replacements
3. [ ] Update main hap.py script for Python 3
4. [ ] Modernize vcfeval integration
5. [ ] Implement improved preprocessing pipeline
6. [ ] Update stratification handling
7. [ ] Create comprehensive test suite
8. [ ] Implement modern packaging
9. [ ] Update documentation

## Resources

- [Truvari project](https://github.com/spiralgenetics/truvari) - Reference for modern genomics code
- [pysam documentation](https://pysam.readthedocs.io/) - For VCF/BAM handling
- [pybedtools documentation](https://daler.github.io/pybedtools/) - For interval operations
- [VCF specification v4.5](https://samtools.github.io/hts-specs/VCFv4.5.pdf)
