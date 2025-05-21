# hap.py Code Cleanup Session Summary

## Project Context

- hap.py is a bioinformatics tool for benchmarking small variant calls
- The codebase is being modernized from Python 2 to Python 3
- We need to identify and remove code not utilized by the core functions

## Core Files

- `src/python/hap.py` - Main variant comparison tool
- `src/python/pre.py` - Preprocessing tool
- `src/python/qfy.py` - Quality filtering tool

## Task Description

The goal is to remove any code or functionality not utilized by the core functions to streamline future development and installation. This requires:

1. Building a mental model/dependency graph of the code structure
2. Identifying which modules are imported and used by the core files
3. Tracing indirect dependencies
4. Identifying unused modules and components
5. Safely removing unused code without breaking core functionality

## Initial Analysis Progress

- Started examining core files (hap.py, pre.py, qfy.py) to understand their imports and dependencies
- Began building a dependency graph to trace what's actually used
- Looked at additional related files like cnx.py and ovc.py to understand their relationship to core functions
- Started examining the Haplo and Tools modules to understand the package structure
- Began tracing imports to identify which components are actually used

## Next Steps

1. Complete the dependency graph/mental model of the codebase
2. Identify with confidence which modules are unused
3. Create a removal plan with specific files and components to be removed
4. Implement the removal in stages, testing after each stage
5. Update any import statements or references as needed
6. Document the changes made

## Development Guidelines

- Follow PEP 8 and use Black for code formatting with 88 character line limit
- Include type hints and docstrings in any modified code
- Test all changes thoroughly using the project's test suite
- Use pre-commit hooks for quality checks
- Document all removals with explanations of why the code was deemed unnecessary
