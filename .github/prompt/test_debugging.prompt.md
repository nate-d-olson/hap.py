# hap.py Test Debugging Strategy

## Overview
This document provides a systematic approach for debugging and fixing failing tests in the modernized hap.py package. The focus is on sustainable fixes that align with Python 3 best practices rather than quick workarounds.

## Prerequisites

### Environment Setup
```bash
# Create the environment if it doesn't exist
micromamba create -n happy python=3.9 -y

# Activate the environment
micromamba activate happy

# Install development dependencies
pip install -e '.[dev,cpp]'

# Verify installation
python -c "import hap_py; print(f'hap.py version: {hap_py.__version__}')"
python -c "import sys; print(f'Python version: {sys.version}')"
```

### Dependency Verification
```bash
# Check for required tools
python -c "
import shutil
import sys

tools = ['bcftools', 'samtools', 'tabix']
missing = [tool for tool in tools if not shutil.which(tool)]
if missing:
    print(f'Missing tools: {missing}', file=sys.stderr)
    print('Please install these tools before proceeding', file=sys.stderr)
    sys.exit(1)
"

# Verify RTG tools installation
python -c "
from pathlib import Path
import sys

rtg_candidates = [
    Path('external/rtg-tools-3.12.1/rtg')  # Updated to match actual location
]

rtg_found = False
for rtg in rtg_candidates:
    if rtg.exists():
        print(f'Found RTG tools: {rtg.absolute()}')
        rtg_found = True
        break

if not rtg_found:
    print('RTG tools not found in expected locations!', file=sys.stderr)
    print('Tests using vcfeval will fail without RTG tools', file=sys.stderr)
"
```

## Running Tests with Proper Tracking

### Unit Tests
```bash
# Run all unit tests with output capture
pytest tests/unit/ -v 2>&1 | tee unit_test_output.txt

# Count failures
echo "Failed tests: $(grep -c 'FAILED' unit_test_output.txt)"
```

### Integration Tests
```bash
# Run all integration tests with output capture and reasonable timeout
python -c "
import subprocess
import time

start_time = time.time()
print(f'Starting integration tests at {time.strftime(\"%Y-%m-%d %H:%M:%S\")}')

try:
    # Set 60-minute timeout for all tests
    proc = subprocess.run(
        ['pytest', 'tests/integration/', '-v'],
        timeout=3600,
        capture_output=True,
        text=True
    )
    with open('integration_test_output.txt', 'w') as f:
        f.write(proc.stdout)
        f.write(proc.stderr)
    print(f'Tests completed in {(time.time() - start_time)/60:.1f} minutes')
except subprocess.TimeoutExpired:
    print('Tests timed out after 60 minutes!')

# Count failures
import subprocess
result = subprocess.run('grep -c \"FAILED\" integration_test_output.txt',
                       shell=True, capture_output=True, text=True)
print(f'Failed tests: {result.stdout.strip()}')
"
```

## Systematic Analysis of Test Failures

### 1. Categorize Failures
```bash
# Group failures by error type
python -c "
import re
from collections import defaultdict
from pathlib import Path

error_patterns = {
    'rtg_not_found': r'(rtg: command not found|rtg.*executable.*not found)',
    'reference_missing': r'(reference.*not found|valid reference path)',
    'temp_dir_conflict': r'directory.*already exists',
    'vcf_validation': r'Invalid header|failed to parse VCF',
    'unimplemented': r'(has been replaced|not yet available|not yet implemented)',
    'python_error': r'AttributeError|ImportError|TypeError'
}

# Read the test output
output = Path('integration_test_output.txt').read_text()

# Find all test failures
failures = re.findall(r'(test_\w+\.py::\w+).*FAILED', output)

# Categorize failures
categories = defaultdict(list)
for test in failures:
    test_failure = output.split(test)[1].split('FAILED')[0]

    matched = False
    for category, pattern in error_patterns.items():
        if re.search(pattern, test_failure, re.IGNORECASE):
            categories[category].append(test)
            matched = True
            break

    if not matched:
        categories['other'].append(test)

# Print summary
print('Failure Categories:')
for category, tests in categories.items():
    print(f'{category}: {len(tests)} tests')
    for test in tests[:5]:  # Show first 5
        print(f'  - {test}')
    if len(tests) > 5:
        print(f'  - ...and {len(tests) - 5} more')
"
```

### 2. Prioritize Fixes
```bash
# Identify most impactful fixes
python -c "
from collections import Counter
from pathlib import Path
import re

# Read the test output
output = Path('integration_test_output.txt').read_text()

# Extract error messages
error_msgs = re.findall(r'E\s+(Error:|Exception:|Fatal:|WARNING:.*failed)', output)
error_counts = Counter(error_msgs)

print('Top 10 Most Common Errors:')
for error, count in error_counts.most_common(10):
    print(f'{count:3d} × {error[:100]}')
"
```

## Implementation Guide for Key Fixes

### 1. RTG Tools Integration

**Sustainable Fix:** Rather than hardcoding paths, implement a robust RTG tools finder that works across environments:

```python
# Add to conftest.py
@pytest.fixture
def rtg_path():
    """Return path to RTG tools with robust fallbacks."""
    from pathlib import Path
    import shutil
    import os

    # Check environment variable first
    if "RTG_PATH" in os.environ and Path(os.environ["RTG_PATH"]).exists():
        return os.environ["RTG_PATH"]

    # Check common locations - updated to match actual locations in project
    rtg_candidates = [
        Path(__file__).parent.parent / "external" / "rtg-tools-3.12.1" / "rtg",
        Path("external/rtg-tools-3.12.1/rtg"),
    ]

    for rtg in rtg_candidates:
        if rtg.exists():
            return str(rtg.absolute())

    # If not found in known locations, check PATH
    rtg = shutil.which("rtg")
    if rtg:
        return rtg

    pytest.skip("RTG tools not found - required for this test")
```

**Usage in Tests:**
```python
def test_vcfeval_example(rtg_path, tmp_path, reference_file):
    """Test vcfeval with proper path handling."""
    result = subprocess.run([
        "python", "-m", "hap_py.hap_py",
        "--engine", "vcfeval",
        "--engine-vcfeval-path", rtg_path,
        "-r", reference_file,
        "-o", str(tmp_path / "output"),
        "truth.vcf", "query.vcf"
    ], capture_output=True, text=True)
    # Test assertions...
```

### 2. Temporary File Management

**Sustainable Fix:** Use pytest's `tmp_path` fixture consistently to avoid cleanup issues:

```python
def test_with_temp_files(tmp_path):
    """Test using temporary directories properly."""
    # Create output directory
    output_dir = tmp_path / "output"
    output_dir.mkdir()

    # Files will be created inside the temporary directory
    output_vcf = output_dir / "output.vcf"

    # Run your test command
    result = subprocess.run([
        "python", "-m", "hap_py.hap_py",
        # other arguments...
        "-o", str(output_vcf),
        # more arguments...
    ], capture_output=True, text=True)

    # Verify results
    assert output_vcf.exists()
    assert result.returncode == 0

    # No explicit cleanup needed - pytest handles tmp_path cleanup automatically
```

### 3. Reference File Management

**Sustainable Fix:** Create a fixture that handles reference file detection:

```python
# Add to conftest.py
@pytest.fixture
def reference_file():
    """Return path to a reference genome file."""
    from pathlib import Path
    import os

    # Check environment variable first
    if "HGREF" in os.environ and Path(os.environ["HGREF"]).exists():
        return os.environ["HGREF"]

    # Check common locations
    ref_candidates = [
        Path("example/happy/resources/chr21.fa"),
        Path("example/happy/resources/test.fa"),
        Path("example/data/reference.fa"),
    ]

    for ref in ref_candidates:
        if ref.exists():
            return str(ref.absolute())

    pytest.skip("Reference genome not found - required for this test")
```

### 4. Python Implementation of C++ Components

**Decision Framework:**

1. For critical functionality that blocks many tests:
   - Implement Python equivalent with same interface
   - Log accurate information about performance differences

2. For less critical functionality:
   - Skip tests with clear markers
   - Document in a central "implementation roadmap" file

**Example Implementation Pattern:**
```python
def multimerge_py(*args, **kwargs):
    """Python implementation of multimerge C++ tool.

    This is a reimplementation of the C++ multimerge tool in Python.
    It may have different performance characteristics but maintains
    functional equivalence for testing purposes.
    """
    # Implementation here
    pass

# In the wrapper script:
def main():
    """Entry point for multimerge command."""
    import sys
    from hap_py.utils.multimerge import multimerge_py

    try:
        sys.exit(multimerge_py(*sys.argv[1:]))
    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(1)
```

### 5. Tracking Fixed Tests

Create a tracking file `.github/TEST_STATUS.md` to track fixed tests:

```markdown
# Test Status Tracking

## Fixed Tests

| Date | Test | Issue | Fix Description |
|------|------|-------|----------------|
| 2023-06-01 | test_integration.py::test_vcfeval | RTG path not found | Added rtg_path fixture |
| 2023-06-01 | test_happy_pg.py::test_happy_pg_basic | Missing reference | Added reference_file fixture |

## Known Issues

| Test | Issue | Status | Priority |
|------|-------|--------|----------|
| test_gvcf_homref.py::test_homref | Requires multimerge | Not implemented | Medium |
```

## Verification Process

After implementing fixes:

1. **Verify fixes incrementally:**
   ```bash
   # Run previously failing tests one by one
   pytest tests/integration/test_specific.py::test_specific -v
   ```

2. **Regression testing:**
   ```bash
   # Run all tests to ensure no regressions
   pytest tests/unit/ tests/integration/ -v
   ```

3. **Code quality validation:**
   ```bash
   # Ensure code meets quality standards
   pre-commit run --all-files
   ```

4. **Update documentation:**
   - Update README with new fixtures
   - Document any new Python implementations
   - Update TEST_STATUS.md

## Best Practices Checklist

- [ ] Fix uses proper Python 3 idioms
- [ ] Fix is maintainable (not hard-coded)
- [ ] Fix maintains compatibility with original behavior
- [ ] Fix has proper error handling
- [ ] Fix has been verified with tests
- [ ] Fix is documented
- [ ] Fix addresses root cause, not just symptoms
