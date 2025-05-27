---
applyTo: "**"
---
# Project General Coding Standards

## Code Formatting and Quality Tools

- Follow language-specific style guides (Python: PEP 8, C++: C++11 style)
- Use Black for Python code formatting with 88 character line limit
- New Python code should include type hints and docstrings
- Use consistent indentation (4 spaces, no tabs)
- Use meaningful variable and function names

### Pre-commit Hooks
The project uses pre-commit hooks to enforce code quality standards:

- **Installation**:
  ```bash
  pip install pre-commit
  pre-commit install
  ```

- **Supported hooks**:
  - `pre-commit-hooks`: Basic checks (trailing whitespace, merge conflicts, etc.)
  - `black`: Automatic code formatting
  - `ruff`: Fast linter with Python 3 compatibility checks
  - `isort`: Import sorting with Black compatibility
  - `pyupgrade`: Automatic upgrades to Python 3.7+ syntax
  - `mypy`: Optional static type checking

- **Running hooks**:
  - Automatically on commit: `git commit -m "Your message"`
  - Manually on all files: `pre-commit run --all-files`
  - On specific files: `pre-commit run --files path/to/file.py`
  - Single hook: `pre-commit run black --files path/to/file.py`

- **CI Integration**:
  - Pre-commit hooks are part of the CI pipeline
  - Failed hooks will cause the CI build to fail

## Error Handling
- Use appropriate exception handling in Python (try/except)
- Use proper error handling in C++ (avoid silent failures)
- Log errors with contextual information about inputs and state
- Validate function inputs, especially for genomic coordinates and data
- Provide meaningful error messages that help diagnose issues

## Code Organization
- Keep functions small and focused on a single task
- Group related functions into logical modules
- Use appropriate access modifiers in C++ (public/private/protected)
- Keep dependencies explicit and minimize global state
- Follow a consistent naming convention for each language

## Testing Best Practices

### Unit Tests
- Each function should have corresponding unit tests
- Use pytest fixtures for common setup/teardown operations
- Add type annotations to test functions for better validation
- Use descriptive test names that explain what's being tested
- Include tests for edge cases and error handling
- When mocking objects:
  - Ensure decorator parameter order matches function signature
  - Use proper parameter names in mock objects to match actual objects
  - Add assertions to verify mock calls are made as expected

### Integration Tests
- Verify interactions between components work correctly
- Set up proper test fixtures and cleanup
- Use relative paths or environment variables for file paths
- Use the `pathlib` module for cross-platform compatibility
- Ensure all required tools and dependencies are available
- Properly handle temporary files and directories
- Use the `rtg_executable` fixture to access RTG tools
- Use the `reference_file` fixture for reference FASTA files
- When testing with `hap.py`, use `--engine-vcfeval-path` with the correct RTG path
- For temporary directories, use `tempfile.mkdtemp()` or pytest's `tmp_path` fixture
- Make sure tests create and clean up their own test data
- Use patching carefully; ensure order of decorators matches function arguments
- For diagnosing integration test failures, check `integration_test_output.txt`
- Handle both string and bytes data types appropriately

## Python 3 Compatibility Notes
- Use type annotations for function parameters and return values
- Convert byte strings to Unicode where appropriate: `ensure_str()` and `ensure_bytes()`
- Use `pathlib.Path` for file path operations instead of `os.path`
- Replace Python 2 dictionary methods with Python 3 equivalents
- Ensure all package folders have proper `__init__.py` files
- Use proper relative imports: `from . import module` instead of `import module`
- Handle iterator differences: `list(map(...))` instead of `map(...)`
