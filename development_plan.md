# Development Plan for hap.py Modernization

## Current Status

We are in Phase 1 of the modernization effort, focusing on build system improvements and dependency management. Key progress includes:

- ✅ Created modern CMake configuration with CMakePresets.json
- ✅ Added package configuration and Find modules
- ✅ Set up CPack configuration
- ✅ Created pyproject.toml with updated dependencies
- ✅ Modernized C++ build with target-based approach
- ✅ Created dedicated CMake modules for standardization
- ✅ Added macOS specific build enhancements
- ✅ Created comprehensive build system documentation
- ✅ Documented overall architecture and components

## Architecture Overview

The hap.py toolkit consists of several interconnected components structured in a layered architecture, designed for modular operation and extensibility. Understanding this architecture is critical for efficient debugging and future development.

### Core Components

1. **Python Frontend Layer**
   - `hap.py`: Main entry point for diploid precision/recall evaluation
     - Orchestrates workflow from input parsing to result generation
     - Manages parallel execution across chromosomes/regions
     - Handles session information and logging
   - `som.py`: Somatic variant comparison (allele-based)
     - Performs bcftools-based allele comparison, ignoring genotypes
     - Specialized for somatic mutation calling evaluation
   - `qfy.py`: Quantification engine for variant statistics
     - Processes comparison outputs into precision/recall metrics
     - Generates stratified performance tables and ROC curves
     - Creates final report summaries in CSV format
   - `pre.py`: Variant preprocessing and normalization
     - Performs variant decomposition, left-alignment, and normalization
     - Converts between chromosome naming conventions
     - Filters variants based on specified criteria
   - Core Python modules in `/src/python/`:
     - `Haplo/`: Diploid comparison utilities
     - `Somatic/`: Somatic comparison utilities
     - `Tools/`: Shared helper functions for I/O, parallelization

2. **C++ Processing Engine**
   - Core library `libhaplotypes` (static) in `/src/c++/lib/`
     - `haplotypes/`: Core variant representation and comparison logic
     - `variant/`: VCF record processing and manipulation
     - `helpers/`: Utility functions shared across components
   - Key executables:
     - `xcmp`: VCF comparison by genotypes (truth vs query)
       - Block-based variant processing with haplotype comparison
       - Genotype-aware variant matching
     - `scmp`: VCF comparison by alleles (ignores zygosity)
       - Allele presence/absence testing
       - Distance-based or exact matching modes
     - `preprocess`: Variant normalization
       - Left-shifting and decomposition of complex variants
       - Reference-based validation and correction
     - `quantify`: Variant counting for statistics
       - Categorizes variants by type (SNP/INDEL/etc.)
       - Computes TP/FP/FN metrics with stratification
     - `hapcmp`: Haplotype-based comparison
       - Graph-based haplotype enumeration and alignment
       - Complex variant resolution through sequence comparison
     - `multimerge`: VCF merging utilities
       - Combines multiple VCF files with specialized options
     - Additional utilities: `hapenum`, `dipenum`, `roc`, `vcfcheck`

3. **Integration Layer**
   - Cython bindings for C++/Python interoperability (though not explicitly visible in repo)
     - Memory management between languages
     - Type conversion between Python and C++ data structures
   - Command-line interfaces in Python call the appropriate C++ executables
   - Temporary file management and serialization/deserialization logic
   - Performance-critical operations delegated to C++ components

### Supporting Components

4. **External Dependencies** (`/external/`)
   - HTSlib: VCF/BAM file handling
     - Core dependency for all variant operations
     - Provides indexing, compression, and fast access
   - Boost: C++ utilities (version 1.58.0 subset)
     - Program options parsing
     - File system operations
     - String manipulation
   - RTG vcfeval: Additional comparison engine
     - Alternative to native hap.py comparison
     - Used for cross-validation of results
   - samtools: SAM/BAM processing
     - Reference sequence handling
     - Alignment file operations
   - bcftools: VCF manipulation
     - VCF merging, filtering, and conversion
     - Header manipulation and annotation

5. **Utility Components**
   - Test framework
     - Unit tests in `/src/c++/test/`
     - Integration tests in shell scripts (`/src/sh/run_integration_test.sh`)
     - Example datasets in `/example/`
   - Build system
     - CMake-based (CMakeLists.txt in multiple directories)
     - External dependency management (`/external/make_dependencies.sh`)
     - Python integration via `install.py`
   - Docker containerization
     - Multiple Dockerfiles for different environments
     - Ensures consistent runtime behavior across systems
   - Continuous Integration pipeline (Jenkins)

### Data Flow & Processing Pipeline

```
                                 ┌─────────────────────┐
                                 │  Input VCF Files    │
                                 │  Truth & Query      │
                                 └─────────┬───────────┘
                                           │
                                           ▼
┌───────────────────────────────────────────────────────────────────────────┐
│ pre.py / preprocess                                                       │
│ ┌─────────────────┐ ┌──────────────┐ ┌────────────────┐ ┌──────────────┐ │
│ │ Decomposition   │ │ Left-shifting│ │ Normalization  │ │ Chr renaming │ │
│ └────────┬────────┘ └───────┬──────┘ └───────┬────────┘ └──────┬───────┘ │
└──────────┼────────────────┬─┼──────────────┬─┼──────────────┬──┼─────────┘
           │                │ │              │ │              │  │
           ▼                ▼ ▼              ▼ ▼              ▼  ▼
┌───────────────────────────────────────────────────────────────────────────┐
│ C++ Comparison Engine                                                     │
│ ┌───────────────────────┐ ┌───────────────────────┐ ┌───────────────────┐ │
│ │ xcmp                  │ │ scmp                  │ │ hapcmp            │ │
│ │ ┌───────────────────┐ │ │ ┌───────────────────┐ │ │                   │ │
│ │ │ Genotype-based    │ │ │ │ Allele-based      │ │ │ ┌─────────────┐   │ │
│ │ │ Block processing  │◄┼─┼─┤ Distance-based    │◄┼─┼─┤ Haplotype   │   │ │
│ │ │ Reference models  │ │ │ │ Complex variants  │ │ │ │ enumeration │   │ │
│ │ └───────────────────┘ │ │ └───────────────────┘ │ │ └─────────────┘   │ │
│ └───────────────────────┘ └───────────────────────┘ └───────────────────┘ │
└────────────────────────────────────┬──────────────────────────────────────┘
                                     │
                                     ▼
┌───────────────────────────────────────────────────────────────────────────┐
│ qfy.py / quantify                                                         │
│ ┌──────────────────┐ ┌───────────────┐ ┌─────────────────┐ ┌────────────┐ │
│ │ TP/FP/FN counts  │ │ Stratification│ │ ROC curve data  │ │ Report     │ │
│ │ by variant type  │ │ by regions    │ │ generation      │ │ formatting │ │
│ └──────────────────┘ └───────────────┘ └─────────────────┘ └────────────┘ │
└────────────────────────────────────┬──────────────────────────────────────┘
                                     │
                                     ▼
                             ┌─────────────────┐
                             │ Output Files    │
                             │ - Summary CSV   │
                             │ - Extended CSV  │
                             │ - ROC data      │
                             │ - VCF with tags │
                             └─────────────────┘
```

### Key Design Principles

1. **Modular Architecture**
   - Clear separation of concerns between components
   - Well-defined interfaces between Python and C++ layers
   - Each executable focuses on a specific task in the pipeline

2. **Performance Optimization**
   - Memory-intensive operations implemented in C++
   - Parallel processing across chromosomes/regions
   - Efficient handling of large VCF files and complex regions

3. **Flexible Comparison Strategies**
   - Multiple comparison engines (xcmp, scmp, hapcmp, vcfeval)
   - Variant matching at different levels of granularity:
     - Genotype-aware (diploid precision/recall)
     - Allele-based (somatic comparison)
     - Sequence-level (haplotype comparison for complex variants)

4. **Robust Error Handling**
   - Input validation in preprocessing steps
   - Graceful handling of malformed variants
   - Detailed logging of processing steps and decisions

5. **Extensibility**
   - Well-defined plugin points for new comparison methods
   - Stratification framework for custom analysis regions
   - Configurable output formats and metrics

### Technical Debt & Modernization Targets

- C++/Python interface modernization (move from raw C API to pybind11)
- Dependency on deprecated Python 2 libraries and syntax
- Inconsistent error handling across components
- Limited type hinting and documentation
- Mix of different build systems and packaging approaches

Understanding these architectural elements is crucial for effective debugging, enhancement, and maintenance of the codebase during the modernization process.

## Phased Development Plan

### Phase 1: Build System & Dependencies (Current)

#### Build System Modernization

- [x] Create CMakePresets.json for better IDE integration
- [x] Add package configuration files (hap.py-config.cmake.in)
- [x] Create modern Find modules for dependencies
- [x] Setup CPack configuration
- [x] Update main CMakeLists.txt with modern practices
- [x] Create dedicated CMake modules for C++ standards and settings
- [x] Create Cython support module for Python-C++ integration
- [x] Add centralized version handling in HappyProperties.cmake
- [x] Add macOS-specific enhancements and code signing support
- [ ] Convert install.py functionality to CMake
- [ ] Add proper install targets
- [ ] Resolve HTSlib detection issues
- [ ] Fix zlib version conflicts

#### Python Package Updates

- [x] Create pyproject.toml
- [ ] Update Python package structure
- [ ] Convert setuptools configuration
- [ ] Add proper entry points
- [ ] Create setup.py wrapper for backward compatibility

### Phase 2: Python 3 Migration

#### Core Conversion

- [ ] Run 2to3 on Python modules
- [ ] Fix syntax changes (print, exceptions)
  - Convert all print statements to print() function
  - Update exception handling to use "as" syntax
  - Replace deprecated **future** imports
- [ ] Update string handling
  - Fix unicode/bytes conversions
  - Add explicit encoding for file operations
  - Update string formatting to f-strings where appropriate
- [ ] Fix integer division
  - Replace / with // for integer division
  - Add explicit casts where needed
- [ ] Update dict/iterator methods
  - Replace iteritems/iterkeys/itervalues with items/keys/values
  - Update obsolete dictionary methods
  - Fix list/iterator behavior differences

#### Dependency Updates

- [ ] Update Python package versions
  - Update pysam to latest compatible version
  - Modernize numpy/scipy/pandas usage
  - Replace deprecated modules with maintained alternatives
  - Add type hint support packages
- [ ] Fix API changes in dependencies
  - Address breaking changes in Python libraries
  - Update Cython bindings for compatibility
  - Fix file handling for Python 3 path objects
- [ ] Update C++ library versions
  - Update HTSlib to latest compatible version
  - Modernize Boost usage (1.74+)
  - Update compression libraries
- [ ] Test compatibility with genomic file formats
  - Validate VCF/BCF parsing in Python 3
  - Test FASTA/BAM file handling
  - Verify compression/decompression in modernized code

### Phase 3: Testing Infrastructure

#### Test Framework

- [ ] Convert shell tests to pytest
- [ ] Add CTest integration
- [ ] Create test data management
- [ ] Set up coverage reporting
- [ ] Add integration test framework

#### CI/CD Setup

- [ ] Create GitHub Actions workflow
- [ ] Configure test automation
- [ ] Set up automated releases
- [ ] Add documentation builds

### Phase 4: Modernization & Cleanup

#### Code Quality

- [ ] Add type hints
- [ ] Update docstrings
- [ ] Apply Black formatting
- [ ] Fix linting issues

#### Performance

- [ ] Profile memory usage
- [ ] Optimize large dataset handling
- [ ] Improve parallelization
- [ ] Benchmark key operations

## Dependencies

### Python

- pysam ≥0.21.0
- numpy ≥1.24.0
- pytest ≥7.0.0
- black ≥23.0.0
- mypy ≥1.0.0
- sphinx ≥7.0.0

### C++

- CMake ≥3.15 (3.20+ recommended)
- Boost ≥1.74.0
  - Required components: filesystem, program_options, system, regex, test
  - Optional components: thread, chrono, iostreams
- HTSlib ≥1.17
  - Requires zlib, libbz2, liblzma
  - Optional: libcurl for remote access
- zlib ≥1.2.13
- jsoncpp ≥1.9.5
- Compiler requirements:
  - GCC 9+ or Clang 10+
  - C++11 support required (C++14 recommended)

## Progress Tracking

| Component          | Status | Notes                                        |
|--------------------|--------|----------------------------------------------|
| CMake Config       | 🟢     | Modern structure implemented                  |
| C++ Build System   | 🟢     | Modernized with dedicated modules            |
| Python Packaging   | 🟡     | pyproject.toml created, more integration needed |
| Build System       | 🟡     | Need install.py replacement and HTSlib fixes |
| Python 3 Port      | ⚪️     | Not started                                  |
| Test Framework     | 🟡     | C++ tests integrated with CTest              |
| Documentation      | 🟡     | Architecture documented, need API docs       |

Legend:

- ⚪️ Not Started
- 🟡 In Progress
- 🟢 Complete  

5. Code Modernization  
   - Add type hints and docstrings  
   - Refactor modules for modularity  
   - Apply code formatting with Black  
6. Test Framework Modernization  
   - Convert bash tests to pytest  
   - Retain integration tests for end-to-end validation  
7. Packaging & Installation  
   - Switch to setuptools or Poetry  
   - Ensure pip-installable package  
8. Documentation & Cleanup  
   - Update user/developer docs  
   - Remove deprecated or unused code

## Dependencies

### Python Dependencies

- pysam ≥0.x, numpy ≥1.x, …

### C++ Dependencies

- htslib ≥1.x  
- Boost ≥1.x

### External Tools

- RTG tools  
- samtools

## Progress Tracking

| Component    | Py3 Conversion | Tests Passing | Notes |
|--------------|----------------|---------------|-------|
| Core         | ⬜             | ⬜            |       |
| IO           | ⬜             | ⬜            |       |
| C++ Bindings | ⬜             | ⬜            |       |
| CLI          | ⬜             | ⬜            |       |

## Challenges & Solutions

- Deep Python–C++ integration: update Cython interfaces  
- External tool API changes: pin compatible versions
