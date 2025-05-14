---
applyTo: "**"
---
# Development Workflow Guidelines

## Getting Started

1. Clone the repository
2. Set up a Python virtual environment
3. Run the installation script with a temporary build directory:
   ```bash
   python install.py /tmp/happy-build
   ```
4. Verify the installation with the test suite:
   ```bash
   cd /tmp/happy-build
   src/sh/run_tests.sh
   ```

## Development Workflow

### Making Code Changes

1. **Branch appropriately** - Create feature branches for specific changes
2. **Test locally** - Always run tests before committing
3. **Document changes** - Update relevant documentation
4. **Follow standards** - Adhere to the coding standards in these instructions

### Testing Process

1. Run unit tests for components you modified
2. Run integration tests to ensure components work together
3. Capture and analyze test failures to identify root causes
4. Document any new test cases you add

### Code Review

1. Self-review your changes before submission
2. Address all comments from reviewers
3. Verify changes with appropriate tests
4. Update documentation to reflect changes

### Debugging

1. Use logging to trace execution flow
2. Create minimal reproducible examples for bugs
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
