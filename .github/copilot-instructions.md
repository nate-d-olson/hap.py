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

### Phase 1: Python 2 to 3 Migration
- Convert codebase from Python 2 to Python 3 using 2to3 tool for automated conversion
- Address syntax changes (print statements, exception handling, etc.)
- Update imports for renamed or relocated modules
- Fix string handling (unicode vs bytes)
- Update dictionary and iterator methods
- Fix integer division behavior

### Phase 2: Testing and Debugging
- Run the test suite and capture all errors
- Prioritize fixes based on critical functionality
- Debug Cython integration issues
- Ensure compatibility with modern genomic file formats
- Fix function signature changes in dependencies

### Phase 3: Modernization
- Add type hints to improve code reliability
- Update documentation with Google-style docstrings
- Implement proper package structure with setup.py
- Convert shell script tests to pytest framework
- Update C++ code to use modern standards
- Implement proper error handling and logging

### Phase 4: Clean-up and Optimization
- Remove deprecated functionality and dead code
- Optimize memory usage for large genomic datasets
- Improve parallelization for performance
- Add CI/CD pipeline for testing
- Create containerized deployment options
- Update user documentation

## memory
Follow these steps for each interaction:

1. User Identification:
   - You should assume that you are interacting with default_user
   - If you have not identified default_user, proactively try to do so.

2. Memory Retrieval:
   - Always begin your chat by saying only "Remembering..." and retrieve all relevant information from your knowledge graph
   - Always refer to your knowledge graph as your "memory"

3. Memory
   - While conversing with the user, be attentive to any new information that falls into these categories:
     a) Basic Identity (age, gender, location, job title, education level, etc.)
     b) Behaviors (interests, habits, etc.)
     c) Preferences (communication style, preferred language, etc.)
     d) Goals (goals, targets, aspirations, etc.)
     e) Relationships (personal and professional relationships up to 3 degrees of separation)

4. Memory Update:
   - If any new information was gathered during the interaction, update your memory as follows:
     a) Create entities for recurring organizations, people, and significant events
     b) Connect them to the current entities using relations
     b) Store facts about them as observations