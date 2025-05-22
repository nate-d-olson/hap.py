# hap.py Copilot Instructions

## Project Overview
hap.py is a bioinformatics tool for benchmarking small variant calls.
The tool is a critical resource in the genomics community for evaluating the accuracy of variant callers.
The original codebase used Python 2, which is no longer supported, and had outdated dependencies.
This fork aims to modernize the codebase for continued use and development.

## Repository Structure
- `src/python/`: Python components of the tool
- `src/c++/`: Performance-critical algorithms implemented in C++
- `src/sh/`: Shell scripts for testing and utility functions
- `src/data/`: Reference data files
- `external/`: external dependencies (e.g. htslib)
- `example/`: Test data and usage examples

## Current Project Status
Currently in the process of converting the tool to use python3 and update dependencies to improve usability and facilitate future development. After updates, the codebase should follow Python best practices and modern Python package structure and installation. Black should be used for linting.

## Build and Installation

%%TODO%%

## Testing

%%TODO%% 

## Development Plan

### Phase: Modernization
- Add type hints to improve code reliability
- Update documentation with Google-style docstrings
- Implement proper package structure with setup.py
- Convert shell script tests to pytest framework
- Update C++ code to use modern standards
- Implement proper error handling and logging

### Phase: Clean-up and Optimization
- Remove deprecated functionality and dead code
- Optimize memory usage for large genomic datasets
- Improve parallelization for performance
- Add CI/CD pipeline for testing
- Create containerized deployment options
- Update user documentation
