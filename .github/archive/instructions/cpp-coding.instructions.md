---
applyTo: "**/*.{cpp,h,hh,hpp}"
---
# C++ Coding Standards for hap.py

## Code Structure
- Use modern C++ features (C++11 or newer) where appropriate
- Follow consistent naming conventions
- Keep functions small and focused
- Minimize global state

## Code Formatting
- Use 4 spaces for indentation (no tabs)
- Maximum line length: 100 characters
- Use consistent brace style (prefer opening brace on the same line)
- Use consistent spacing around operators

## Documentation
- Document all public APIs
- Use doxygen-style comments for functions, classes, and methods
- Include parameter descriptions, return values, and exceptions

## Memory Management
- Prefer smart pointers (std::unique_ptr, std::shared_ptr) over raw pointers
- Clean up resources appropriately
- Be mindful of ownership semantics

## Error Handling
- Use exceptions for exceptional conditions
- Return error codes or optional values for expected failure modes
- Document error conditions

## Performance Considerations
- Be mindful of performance in critical bioinformatics algorithms
- Consider memory usage with genomic data
- Optimize hot paths where necessary

## Compatibility
- Ensure code compiles on target platforms (Linux, macOS)
- Be mindful of compiler differences
- Use CMake for cross-platform building

## Dependencies
- Minimize external dependencies
- Document third-party library usage
- Use consistent versions of dependencies

## Testing
- Write unit tests for C++ code
- Test edge cases and error conditions
- Update tests when changing existing code
