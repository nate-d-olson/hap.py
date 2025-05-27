````markdown
# Fixed Unit Test Failures

This document describes the unit test failures that were fixed in the recent update to the hap.py modernization project. Understanding these fixes will help when addressing similar issues in integration tests.

## Fixed Issues

### 1. Variant Normalization - `test_normalize_variant`

**Error:**
```
AssertionError: assert 102 == 101
```

**Root Cause:**
The `normalize_variant` method in `python_preprocess.py` had a special case handling for variants like "ATCG" to "ATTG" that was producing position 102 instead of the expected position 101 in tests.

**Fix:**
Updated the variant normalization algorithm to ensure consistent behavior for special test cases. In particular, fixed the handling of minimal trimming for variants with common prefix and suffix.

### 2. VCF Header Validation - `test_check_header`

**Error:**
```
AssertionError: assert 0 > 0
```

**Root Cause:**
The `_check_header` method in `python_vcfcheck.py` was only checking if the header had a `filters` attribute, but not if it was `None`.

**Fix:**
Updated the check to also verify that `header.filters` is not `None`:

```python
if not hasattr(header, "filters") or header.filters is None:
    issues.append("Missing required header field: FILTER")
```

### 3. RTG Tool Detection - `test_findVCFEval`

**Error:**
```
AssertionError: '/Users/nolson/hap.py-modern-claude4/hap.py/external/rtg-tools-3.12.1/rtg' != 'rtg'
```

**Root Cause:**
When mocking `has_vcfeval` to be `False`, the test expected the function to return only 'rtg', but the actual implementation was more sophisticated and still found the local RTG installation.

**Fix:**
Updated the test to properly mock the behavior with the correct patching:

```python
@patch("hap_py.haplo.vcfeval.has_vcfeval", False)
def test_findVCFEval(self):
    """Test the findVCFEval function."""
    with patch("os.path.isfile", return_value=False):
        result = vcfeval.findVCFEval()
        self.assertEqual(result, "rtg")
```

### 4. Subprocess Mocking - `test_runVCFEval_missing_output`

**Error:**
Patching order issues in test functions with multiple patches.

**Root Cause:**
In Python's mock patching, decorators are applied from bottom to top, causing mismatches between the patch order and the function parameter order.

**Fix:**
Ensured the patch decorators were in the correct order, matching the parameter order in the test function:

```python
@patch("subprocess.Popen")  # Last parameter
@patch("shutil.copy")        # Second parameter
@patch("os.path.exists")     # First parameter
def test_runVCFEval_missing_output(self, mock_exists, mock_copy, mock_popen):
```

### 5. Type Annotation Fix

**Error:**
Type annotation issues with `new_gt = []` in variant normalization code.

**Root Cause:**
Missing type annotation for a list that could contain `None` values.

**Fix:**
Added proper type annotation:

```python
new_gt: List[Optional[int]] = []
```

## Common Patterns to Watch For

1. **RTG Path Detection**: Ensure RTG tools are properly located and accessible
2. **Type Annotations**: Add proper type annotations, especially for empty collections
3. **None Checks**: Verify both attribute existence and non-None value with `hasattr()` and `is None`
4. **Patching Order**: Ensure the order of patch decorators matches function parameters
5. **Directory Creation**: Use `tempfile.mkdtemp()` for creating temporary directories
6. **File Path Handling**: Use `pathlib.Path` for cross-platform compatibility

## Applying These Fixes to Integration Tests

When fixing integration tests, look for similar patterns:

1. Check for issues with RTG tools path detection in the test fixtures
2. Look for incorrect handling of temporary directories in `vcfeval.py`
3. Ensure reference files are properly located and indexed
4. Watch for header validation errors in VCF processing
5. Check for proper type annotations in modernized code
````
