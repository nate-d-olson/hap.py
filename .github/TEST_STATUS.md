# hap.py Test Status Tracking

## Current Status Summary

- **Unit Tests**: Most passing, remaining failures relate to Python/C++ porting
- **Integration Tests**: Several failures due to path issues, reference handling, and unimplemented C++ components
- **RTG Integration**: Mostly fixed but requires standardized path handling

## Fixed Test Categories

| Category | Issue | Fix Implemented | Status |
|----------|-------|-----------------|--------|
| RTG Tools Path Detection | "rtg: command not found" errors | Updated `findVCFEval()` in vcfeval.py | ✅ Complete |
| SDF Template Directory | "directory already exists" errors | Replaced with `tempfile.mkdtemp()` | ✅ Complete |
| Test Package Structure | Module import errors | Added missing `__init__.py` files | ✅ Complete |
| RTG Detection Warning | False warnings about missing RTG | Updated `init()` in `__init__.py` | ✅ Complete |

## Known Issues

| Test Pattern | Issue | Status | Priority |
|--------------|-------|--------|----------|
| test_gvcf_homref.py and similar | Requires `multimerge` functionality | Needs Python implementation | High |
| Various reference file tests | Missing reference file configuration | Needs standardized fixture | Medium |
| VCF header validation | Overly strict validation | Needs updated validator | Low |

## Implementation Roadmap

### 1. Multimerge Python Implementation
The `multimerge` C++ tool needs a Python equivalent that maintains the same interface:

```python
def multimerge_py(*args, **kwargs):
    """Python implementation of multimerge C++ tool."""
    # TODO: Implement functionality
    pass
```

### 2. Standard RTG Path Handling
All tests should use the `rtg_path` fixture instead of hardcoded paths:

```python
def test_example(rtg_path, tmp_path):
    result = subprocess.run([
        "python", "-m", "hap_py.hap_py",
        "--engine", "vcfeval",
        "--engine-vcfeval-path", rtg_path,
        # other arguments...
    ])
```

### 3. Reference File Handling
Use the `reference_file` fixture consistently:

```python
def test_with_reference(reference_file, tmp_path):
    result = subprocess.run([
        "python", "-m", "hap_py.hap_py",
        "-r", reference_file,
        # other arguments...
    ])
```
