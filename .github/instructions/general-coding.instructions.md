---
applyTo: "**"
---
# Project General Coding Standards

## Code Formatting

- Follow language-specific style guides (Python: PEP 8, C++: C++11 style)
- Use Black for Python code formatting with 88 character line limit
- New Python code should include type hints and docstrings
- Use consistent indentation (4 spaces, no tabs)
- Use meaningful variable and function names

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
