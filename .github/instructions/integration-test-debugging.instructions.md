---
applyTo: "tests/integration/**"
---
# Integration Test Debugging Guide

## Integration Test Overview

The hap.py modernization project includes extensive integration tests that verify the functionality of the modernized Python 3 code. These tests cover core features like variant comparison, preprocessing, and filtering. The tests are located in `tests/integration/` and use example data from `example/integration/`.

## Common Integration Test Failures

### 1. RTG Tools SDF Directory Conflicts

**Error Pattern**:
```
Error: The directory "/var/folders/.../vcfeval.sdf.xyz" already exists. Please remove it first or choose a different directory.
```

**Root Cause**: The SDF template creation logic in `vcfeval.py` may create temporary directories that conflict with RTG's expectations.

**Fix Strategy**:
1. Fix the temporary directory creation logic in `vcfeval.py` to use `tempfile.mkdtemp()` properly
2. Ensure the directory is properly created and cleaned up between test runs
3. Use a unique naming pattern to avoid conflicts

**Example Fix**:
```python
# Create a unique temporary directory for the SDF template
template_dir = tempfile.mkdtemp(
    dir=args.scratch_prefix,
    prefix="vcfeval.sdf.",
)
args.engine_vcfeval_template = template_dir
```

### 2. Reference File Availability

**Error Pattern**:
```
FileNotFoundError: Please specify a valid reference path using -r.
```

**Root Cause**: Tests not specifying reference files when required.

**Fix Strategy**:
1. Update tests to use the `reference_file` fixture from `conftest.py`
2. Ensure proper reference file paths are provided to tests
3. Validate reference files have proper indexes (.fai files)

**Example Fix**:
```python
@pytest.fixture(scope="session")
def reference_file(example_dir):
    """Provide reference file for tests."""
    # Check environment variable first
    if "HGREF" in os.environ:
        return os.environ["HGREF"]

    # Use chr21.fa from example directory
    ref_path = example_dir / "chr21.fa"
    if ref_path.exists():
        return str(ref_path)
```

### 3. RTG Path Detection

**Error Pattern**:
```
WARNING: Executable for rtg not found
```

**Root Cause**: Tests unable to find the RTG executable.

**Fix Strategy**:
1. Update the `get_rtg_path()` function in `conftest.py` to correctly detect RTG tools
2. Add support for environment variable RTG_PATH
3. Add diagnostic output to help identify RTG path issues
4. Use `shutil.which("rtg")` to check for RTG in the PATH

**Example Fix**:
```python
def get_rtg_path():
    """Get the RTG tools path for integration tests."""
    project_root = Path(__file__).parent.parent

    # Check environment variable first
    if "RTG_PATH" in os.environ and Path(os.environ["RTG_PATH"]).exists():
        return os.environ["RTG_PATH"]

    # Try external directory (modernized location)
    rtg_path = project_root / "external" / "rtg-tools-3.12.1" / "rtg"
    if rtg_path.exists():
        print(f"Found RTG tools at {rtg_path}")
        return str(rtg_path)
```

### 4. VCF Header Validation

**Error Pattern**:
```
ERROR: Error checking file: Invalid header
Missing required header field: FILTER
```

**Root Cause**: VCF header validation may be too strict or inconsistent.

**Fix Strategy**:
1. Update the `_check_header` method in `VCFChecker` to check both `hasattr` and `is None`
2. Ensure tests that expect FILTER validation use strict mode
3. Update mock headers in tests to have appropriate fields

**Example Fix**:
```python
if not hasattr(header, "filters") or header.filters is None:
    issues.append("Missing required header field: FILTER")
```

### 5. Variant Normalization Issues

**Error Pattern**:
```
AssertionError: assert 102 == 101
```

**Root Cause**: Inconsistent variant normalization behavior, particularly with edge cases.

**Fix Strategy**:
1. Examine the `normalize_variant` implementation for edge cases
2. Check for off-by-one errors in position calculations
3. Update test expectations to match actual behavior if appropriate

### 6. Missing Modernized Tools

**Error Pattern**:
```
ERROR: multimerge has been replaced with Python modules
```

**Root Cause**: Placeholder scripts for tools that haven't been fully modernized.

**Fix Strategy**:
1. Implement Python equivalents for missing tools
2. Update tests to use the Python implementations directly
3. Skip tests for functionality that isn't yet available

## Strategic Test Debugging Process

### Step 1: Identify Failure Patterns

Run integration tests to identify common errors:
```bash
micromamba activate happy
pytest tests/integration/ -v | tee integration_test_output.txt
```

Review `integration_test_output.txt` to identify common error patterns.

### Step 2: Prioritize Fixes

1. **First Fix**: SDF Template Directory Conflicts (in `vcfeval.py`)
2. **Second Fix**: RTG Path Detection (in `conftest.py`)
3. **Third Fix**: Reference File Availability (in test fixtures)
4. **Fourth Fix**: VCF Validation Issues (in `python_vcfcheck.py`)
5. **Fifth Fix**: Tool-specific Implementations (as needed)

### Step 3: Update Test Configuration

Ensure proper test fixtures and configuration:
```python
@pytest.fixture
def run_happy(rtg_executable, reference_file, tmp_path):
    """Run hap.py with proper configuration."""
    def _run(truth_vcf, query_vcf, out_prefix, extra_args=None):
        cmd = [
            "python", "bin/hap.py", truth_vcf, query_vcf,
            "-r", reference_file,
            "--engine", "vcfeval",
            "--engine-vcfeval-path", rtg_executable,
            "-o", os.path.join(tmp_path, out_prefix),
        ]
        if extra_args:
            cmd.extend(extra_args)
        subprocess.check_call(cmd)
    return _run
```

### Step 4: Fix and Test Iteratively

1. Make one fix at a time and run specific tests to verify the fix works
2. Prioritize fixes that address multiple failures
3. Use the pytest `-k` option to run specific tests:
```bash
pytest tests/integration/test_integration.py -k "test_name" -v
```

### Step 5: Clean Up and Document

1. Remove any temporary files created during debugging
2. Document any edge cases or special handling required
3. Update tests to reflect correct behavior

## Test Verification Checklist

Before committing fixes:

1. ✅ Verify RTG tools are correctly located
2. ✅ Confirm reference files are available and indexed
3. ✅ Check that temporary files are cleaned up
4. ✅ Ensure all tests are passing or explicitly skipped
5. ✅ Run linting and formatting tools: `pre-commit run --all-files`
