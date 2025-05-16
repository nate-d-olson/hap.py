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
