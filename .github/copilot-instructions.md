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
Currently the install process is handled by the `install.py` script. Which accepts a build directory as a command line argument. Use a temp directory when building and testing. Do not make edits to the build directory as they will be lost next time the install script is run. The tool uses Cython. CMake is used for installing C++ dependencies and libraries and building source code.

## Testing
After the install completes, the install script tests the install using a series of bash and Python scripts in the `src/sh` directory which are executed by the `src/sh/run_tests.sh` script. The test scripts should be run without failing on error so that all errors can be identified. As the scripts contain multiple tests with ambiguous names, the standard error and standard out should be recorded and used to determine the individual failing tests and specific parts of the codebase that are failing the test.

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
