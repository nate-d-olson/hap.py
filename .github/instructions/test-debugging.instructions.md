---
applyTo: "tests/**"
---

# Test Debugging Guide

This document provides detailed strategies for debugging test failures in the hap.py codebase, focusing on common issues encountered during the Python 3 migration.

## Common Test Failures and Solutions

### Unit Test Failures

#### 1. Mock Patching Issues

**Problem:** Test using `@patch` decorators fails with parameters not matching expectations.

**Solution:**
- Remember that parameters are passed to the test method in reverse order of decorator application
- Example:
  ```python
  # CORRECT ORDER:
  @patch("subprocess.Popen")
  @patch("shutil.copy")
  @patch("os.path.exists")
  def test_method(self, mock_exists, mock_copy, mock_popen):
      # ...

  # INCORRECT ORDER:
  @patch("os.path.exists")
  @patch("shutil.copy")
  @patch("subprocess.Popen")
  def test_method(self, mock_popen, mock_copy, mock_exists):
      # This will fail because the parameters don't match the decorators
  ```

#### 2. VCF Header Validation

**Problem:** Tests expecting specific header validation behavior fail when using strict vs. non-strict mode.

**Solution:**
- Ensure the `_check_header` method checks for required fields consistently
- Use separate message formats for strict vs. non-strict modes
- Always check critical fields like FILTER regardless of strict mode

#### 3. Variant Normalization

**Problem:** Implementation vs. test expectations mismatch in how variants are normalized.

**Solution:**
- Be aware of maximal vs. minimal trimming behavior
- Check position adjustments, especially for off-by-one errors
- Document expected behavior in both implementation and test

### Integration Test Failures

#### 1. RTG Tools Path Issues

**Problem:** Tests fail with "rtg: command not found" or similar errors.

**Solution:**
- Ensure `findVCFEval()` correctly locates RTG executable:
  ```python
  # Update to check external directories:
  repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
  ext_rtg = os.path.join(repo_root, "external", "rtg-tools-3.12.1", "rtg")
  if os.path.isfile(ext_rtg) and os.access(ext_rtg, os.X_OK):
      return ext_rtg
  ```
- Pass explicit path with `--engine-vcfeval-path` in tests

#### 2. Temporary Directory Issues

**Problem:** Tests fail with "File exists" or permission errors when creating directories.

**Solution:**
- Use `tempfile.mkdtemp()` for directory creation instead of `NamedTemporaryFile`
- Properly clean up temporary directories in both success and error paths
- Use `try/finally` blocks to ensure cleanup happens

#### 3. File Path Compatibility

**Problem:** Tests fail with file not found errors or path issues.

**Solution:**
- Use `pathlib.Path` for cross-platform path handling
- Ensure all test files have proper indexes (.fai, .tbi)
- Use relative paths from repository root:
  ```python
  repo_root = Path(__file__).resolve().parent.parent.parent
  example_path = repo_root / "example" / "example.vcf.gz"
  ```

## Step-by-Step Debugging Process

### For Unit Tests

1. Run specific test with verbose output:
   ```bash
   pytest tests/unit/test_specific.py::TestClass::test_specific_method -v
   ```

2. Add print statements to understand execution flow:
   ```python
   print(f"DEBUG: {variable_name} = {variable_value}")
   ```

3. Check for type mismatches:
   ```python
   print(f"DEBUG: {type(actual_result)} vs {type(expected_result)}")
   ```

4. For mock objects, verify calls are made correctly:
   ```python
   mock_object.assert_called_once_with(expected_arg1, expected_arg2)
   ```

### For Integration Tests

1. Run with output capture:
   ```bash
   pytest tests/integration/test_specific.py::test_name -v | tee debug_output.txt
   ```

2. Set environment variables for tools:
   ```bash
   export RTG_PATH="/path/to/rtg" && pytest tests/integration/test_specific.py
   ```

3. Run with pytest debug mode for PDB access:
   ```bash
   pytest tests/integration/test_specific.py -vvs --no-header --showlocals
   ```

4. Trace command execution in tests:
   ```python
   import subprocess

   # Add this before the subprocess call
   print(f"DEBUG: Running command: {' '.join(cmd)}")
   result = subprocess.run(cmd, capture_output=True, text=True)
   print(f"DEBUG: stdout: {result.stdout}")
   print(f"DEBUG: stderr: {result.stderr}")
   ```

## Known Issues and Workarounds

1. **SDF Template Directory Creation**
   - Issue: RTG format fails if directory exists
   - Solution: Use `mkdtemp` and ensure proper cleanup

2. **RTG Detection Warning**
   - Issue: "WARNING: Executable for rtg not found" despite RTG being available
   - Solution: Update `init()` to check project's RTG location

3. **Test Package Structure**
   - Issue: `ModuleNotFoundError: No module named 'tests.utils'`
   - Solution: Create missing `__init__.py` files in test directories

4. **Reference File Issues**
   - Issue: Tests fail with reference file errors
   - Solution: Verify file existence and proper indexing
