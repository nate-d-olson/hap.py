---
applyTo: "**"
---
# Development Workflow Guidelines

## %%TODO%% Install and Build Process

## %%TODO%% Debugging and Fixing Build Issues

## Development Workflow

### Making Code Changes

1. Branch approach:

   * Create a fix-build-system branch for CMake and build fixes
   * Create a python3-migration branch for Python code updates
   * Work on these branches in parallel when possible

2. Testing workflow:

   * Test build system fixes with existing Python 2 code first
   * Test Python 3 updates with mock C++ interfaces when needed
   * Integrate both changes when each part is stable

3. Code quality workflow:

  * Use pre-commit for automated checks: `pre-commit run --files <changed_files>`
  * Fix any issues reported by pre-commit
  * Run specific hooks as needed:

```bash
# Format specific files with Black
pre-commit run black --files src/python/path/to/file.py

# Run specific Python 3 compatibility checks
pre-commit run pyupgrade --files src/python/path/to/file.py

# Run linting with auto-fixes
pre-commit run ruff --files src/python/path/to/file.py
```

### Testing Process

1. Run build verification tests after any CMake changes
2. Run unit tests for components you modified
3. Run integration tests to ensure components work together
4. Capture and analyze test failures to identify root causes
5. Document any new test cases you add

### Code Review

1. Self-review your changes before submission
2. Address all comments from reviewers
3. Verify changes with appropriate tests
4. Update documentation to reflect changes

### Debugging

1. Use logging to trace execution flow
2. For build issues, use CMake's verbose output and message() commands
3. Test fixes thoroughly before committing
4. Document root causes of significant bugs for future reference

## Release Process

1. Update version numbers in relevant files
2. Update the RELEASES.md file with changes
3. Ensure all tests pass
4. Build and verify the release package
5. Create a tagged release in the repository

## Documentation

1. Update API documentation for any changed functions
2. Keep the main README up-to-date
3. Document any new features or significant changes
4. Ensure installation instructions remain accurate
