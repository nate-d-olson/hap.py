---
applyTo: "**"
---
# Development Workflow Guidelines

## Environment Setup

**ALWAYS activate the micromamba environment before starting any development work:**
```bash
micromamba activate happy-dev
```

## %%TODO%% Install and Build Process

## %%TODO%% Debugging and Fixing Build Issues

## Development Workflow

### Making Code Changes

1. **Environment activation:**
   * Always start with: `micromamba activate happy-dev`
   * Verify correct environment: `which python` should show `/Users/nolson/micromamba/envs/happy-dev/bin/python`

2. Branch approach:

   * Create a fix-build-system branch for CMake and build fixes
   * Create a python3-migration branch for Python code updates
   * Work on these branches in parallel when possible

3. Testing workflow:

   * Test build system fixes with existing Python 2 code first
   * Test Python 3 updates with mock C++ interfaces when needed
   * Integrate both changes when each part is stable

4. Code quality workflow:

  * Use pre-commit for automated checks: `pre-commit run --files <changed_files>`
  * Fix any issues reported by pre-commit
  * Run specific hooks as needed:

```bash
# Always ensure environment is active first
micromamba activate happy-dev

# Format specific files with Black
pre-commit run black --files src/python/path/to/file.py

# Run specific Python 3 compatibility checks
pre-commit run pyupgrade --files src/python/path/to/file.py

# Run linting with auto-fixes
pre-commit run ruff --files src/python/path/to/file.py
```

### Testing Process

1. **Always activate environment first:** `micromamba activate happy-dev`
2. Run build verification tests after any CMake changes
3. Run unit tests for components you modified
4. Run integration tests to ensure components work together
5. Capture and analyze test failures to identify root causes
6. Document any new test cases you add

### Code Review

1. Self-review your changes before submission
2. Address all comments from reviewers
3. Verify changes with appropriate tests
4. Update documentation to reflect changes

### Debugging

1. **Ensure correct environment:** `micromamba activate happy-dev`
2. Use logging to trace execution flow
3. For build issues, use CMake's verbose output and message() commands
4. Test fixes thoroughly before committing
5. Document root causes of significant bugs for future reference

### Debugging Test Failures

#### Unit Test Failures

1. **Activate environment:** `micromamba activate happy-dev`
2. Identify the specific failing test and error message
3. Examine the test's expected vs. actual behavior
4. Check for:
   - Decorator parameter order in patched tests
   - Type annotation issues
   - Implementation vs. test expectation mismatches
   - Inconsistent behavior between strict and non-strict modes

#### Integration Test Failures

1. **Activate environment:** `micromamba activate happy-dev`
2. Run integration tests with output capture:
   ```bash
   micromamba activate happy-dev
   pytest tests/integration/ -v | tee integration_test_output.txt
   ```

2. For RTG-related failures:
   - Ensure RTG tools is properly located and accessible
   - Verify the `get_rtg_path()` function in `tests/conftest.py` correctly locates RTG tools
   - Check that `--engine-vcfeval-path` points to the correct RTG executable (usually at `/Users/nolson/hap.py-modern-claude4/hap.py/external/rtg-tools-3.12.1/rtg`)
   - Ensure tests are passing the RTG path via the `rtg_executable` fixture
   - Look for "SDF directory already exists" errors from RTG format commands
   - Validate that `@patch` decorators for RTG-related functions are applied in the correct order

3. For temporary file/directory issues:
   - Confirm proper cleanup in test fixtures
   - Use `tempfile.mkdtemp()` instead of `tempfile.NamedTemporaryFile()` for directories
   - Look for issues in `vcfeval.py` related to SDF template directory creation
   - Check for permission or cross-filesystem issues with temporary directories
   - Ensure test fixtures are using `tmp_path` to create isolated test directories

4. For reference file errors:
   - Verify reference FASTA files have proper `.fai` indexes
   - Check that VCF files have appropriate `.tbi` indexes
   - Ensure tests are using the `reference_file` fixture from `conftest.py`
   - Use `pathlib.Path` for cross-platform path handling
   - Check for missing environment variables like `HGREF` that might be needed

5. For VCF parsing and validation issues:
   - Look for header validation errors in `_check_header` method
   - Check for proper detection of required fields like FILTER
   - Ensure VCF preprocessing handles AC field values correctly
   - Verify that the `normalize_variant` method behaves consistently

6. For missing modernized tools:
   - Check if placeholder scripts like `multimerge` need implementation
   - Ensure binary wrapper scripts are available in the `build/bin` directory
   - Look for "has been replaced with Python modules" error messages

7. Systematic debugging approach:
   - Start by fixing the SDF template directory issues in `vcfeval.py`
   - Then address reference file availability using the `reference_file` fixture
   - Update tests to use proper RTG path detection
   - Fix VCF header validation issues
   - Finally, address any remaining test-specific failures

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
