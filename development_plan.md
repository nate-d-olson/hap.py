# Development Plan for hap.py Modernization

## Overview
This document outlines the roadmap for migrating hap.py from Python 2 to Python 3, modernizing the codebase, and improving usability.

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

## Development Steps
1. Initial Assessment & Environment Setup  
   - Create a reproducible Docker environment  
   - Document current functionality with examples  
   - Establish baseline tests using the existing Python 2 code  
2. Python 2 → Python 3 Conversion  
   - Run 2to3 on core modules  
   - Fix syntax and dependency API changes  
3. Dependency Updates  
   - Update C++ libraries (htslib, Boost)  
   - Upgrade Python packages to Python 3–compatible versions  
   - Adjust CMake and build scripts  
4. Iterative Testing & Debugging  
   - Track test results and error logs  
   - Debug common failures in unit and integration tests  
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
