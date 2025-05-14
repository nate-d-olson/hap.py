# hap.py - Bioinformatics Variant Calling Benchmarking Tool

hap.py is a bioinformatics tool for benchmarking small variant calls. It analyzes and compares variant calls from different variant callers against a truth set to measure performance metrics like precision and recall. The tool is widely used in the genomics community despite its legacy codebase.

## Project Goal

We are modernizing this codebase by:

1. Converting from Python 2 to Python 3
2. Updating outdated dependencies
3. Implementing modern Python packaging standards
4. Improving installation process
5. Adding type hints and better documentation
6. Replacing ad-hoc testing with proper frameworks
7. Optimizing performance for genomic datasets

## Tech Stack

- Python 3.6+ (migrating from Python 2)
- C++ (core algorithms and performance-critical components)
- CMake (build system)
- Cython (Python/C++ interface)
- Bcftools/htslib (for VCF/BAM handling)
- Boost C++ libraries (subset)
- Pandas, NumPy (data analysis)

## Key Considerations

- Bioinformatics data tends to be large - consider memory usage
- Maintain backward compatibility with existing workflows
- Preserve core algorithms while modernizing the code around them
- Focus on reliability and correctness - this is scientific software
- Make installation process more user-friendly

## Development Environment

- Use `install.py` to build the tool in a temporary directory
- Always run the full test suite after changes
- Implement proper version control practices
- Follow the coding standards in the instructions files
