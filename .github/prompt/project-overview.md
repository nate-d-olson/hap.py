# hap.py - Bioinformatics Variant Calling Benchmarking Tool

hap.py is a bioinformatics tool for benchmarking small variant calls. It analyzes and compares variant calls from different variant callers against a truth set to measure performance metrics like precision and recall. The tool is widely used in the genomics community despite its legacy codebase.

## Project Goal

We are modernizing this codebase by:

1. Updating outdated dependencies
1. Implementing modern Python packaging standards
1. Improving installation process
1. Adding type hints and better documentation
1. Replacing ad-hoc testing with proper frameworks
1. Optimizing performance for genomic datasets

## Tech Stack

- Python 3.6+ (migrating from Python 2)
- C++ (core algorithms and performance-critical components)
- CMake (build system)
- Cython (Python/C++ interface)

## Key Considerations

- Maintain backward compatibility with existing workflows
- Preserve core algorithms while modernizing the code around them
- Focus on reliability and correctness - this is scientific software
- Make installation process more user-friendly

## Development Environment

- Use `install.py` to build the tool in a temporary directory
- Always run the full test suite after changes
- Implement proper version control practices
- Follow the coding standards in the instructions files
