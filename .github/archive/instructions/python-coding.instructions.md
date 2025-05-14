---
applyTo: "**/*.py"
---
# Python Coding Standards for hap.py

## Python Version
- All code should be compatible with Python 3.6+
- No Python 2-specific code should remain in the codebase

## Code Formatting
- Use Black for code formatting
- Maximum line length: 88 characters (Black default)
- Use 4 spaces for indentation (no tabs)
- Follow PEP 8 style guidelines

## Documentation
- All functions, methods, and classes should have docstrings
- Use Google-style docstrings
- Include parameter types, return types, and descriptions
- Document exceptions that may be raised

## Type Hinting
- Add type hints to all new code
- Gradually add type hints to existing code during refactoring
- Use typing module for complex types

## Imports
- Group imports in the following order:
  1. Standard library imports
  2. Related third-party imports
  3. Local application/library specific imports
- Sort imports alphabetically within each group
- Use absolute imports over relative imports

## Error Handling
- Use specific exception types rather than generic exceptions
- Provide context in exception messages
- Use try/except blocks appropriately
- Log exceptions with proper context

## Testing
- Write unit tests for new functionality
- Update tests when modifying existing code
- Aim for high test coverage
- Use pytest for testing framework

## Dependencies
- Keep dependencies minimal and up-to-date
- Document all dependencies in requirements files
- Use virtual environments for development

## Performance Considerations
- Consider memory usage for bioinformatics data
- Optimize computation-heavy functions
- Use appropriate data structures for genomic data
