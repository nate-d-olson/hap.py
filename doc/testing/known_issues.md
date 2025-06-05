# Known Testing Issues

This document outlines known issues with the hap.py test suite that require attention.

## Integration Test Issues

### 1. Integration Tests Hanging/Unresponsive During Execution

**Priority**: HIGH  
**Status**: Unresolved

**Problem**: The integration test suite appears to hang or become unresponsive when running with `pytest tests/integration/`. This makes it difficult to complete the full test validation process.

**Observed Behavior**:
- Running `pytest tests/integration/ -v` causes the terminal to become unresponsive
- Individual tests may work but the full suite execution hangs
- Commands like `python -m pytest tests/integration/test_fastasize.py -v` don't return results in reasonable time

**Potential Causes**:
- Resource contention issues with temporary file creation
- Blocking I/O operations that don't timeout properly
- Process deadlocks in subprocess calls to external tools
- Insufficient cleanup of background processes

**Suggested Solutions**:
1. Add timeout mechanisms to all subprocess calls
2. Implement proper cleanup in test teardown methods
3. Use pytest-timeout plugin to prevent hanging tests
4. Review temporary file handling patterns

### 2. GA4GH Integration Tests Failing with Import/Setup Errors

**Priority**: MEDIUM  
**Status**: Partially resolved (unit tests fixed)

**Problem**: GA4GH-related integration tests fail with import errors and module setup issues.

**Observed Behavior**:
- Import errors for GA4GH compliance modules
- Setup/teardown issues in GA4GH integration tests
- Some floating point precision issues in F1 score calculations (fixed in unit tests)

**Resolution**: Unit tests have been fixed with proper floating point precision handling. Integration tests may still need work.

### 3. VCF Processing Integration Tests Failing

**Priority**: MEDIUM  
**Status**: Unresolved

**Problem**: Various VCF processing integration tests are failing, potentially due to:
- File format compatibility issues
- Path resolution problems with external tools
- Data validation failures

**Recommended Action**: These issues should be addressed systematically by repository maintainers with appropriate development environment setup.

## Test Suite Status

- **Unit Tests**: ✅ 71 passed, 5 skipped, 0 failed
- **Integration Tests**: ⚠️ Issues prevent full validation
- **Performance Tests**: 🔄 Require integration test fixes first

## For Developers

When working on these issues:
1. Always use the `happy-dev` micromamba environment
2. Ensure RTG tools are properly configured
3. Check `.github/instructions/debugging_guide.md` for common solutions
4. Run unit tests first to ensure base functionality works
5. Use individual integration test files for debugging rather than the full suite
