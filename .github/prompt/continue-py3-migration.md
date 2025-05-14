# Continue Python 3 Migration: hap.py

## Migration Status Summary

Based on current progress (Updated May 14, 2025):

- **Overall Progress**: 43.8% of files fully migrated (21/48 files)
- **Remaining Issues**: 90 issues across 27 files
- **Priority Areas**: Haplo module (32 issues), Tools module (25 issues), Somatic module (12 issues)

### Issues by Type
1. String/Unicode Issues (78 occurrences)
   - Unicode to str conversion
   - File I/O encoding specification

2. Exception Syntax (12 occurrences)
   - Update from `except X, e` to `except X as e`

### Completed Work

- **Cython Integration**:
  - Created `src/cmake/CythonSupport.cmake` with Python 3 support
  - Implemented proper module structure in `src/python/Haplo/cython/`
  - Added fallback mechanism with mock implementations
  - Updated build system integration in CMakeLists files
  - Implemented string handling with proper encoding/decoding

- **Build System**:
  - Created Python 3 compatible `CMakeLists.txt.py3`
  - Added Python version detection
  - Updated paths for Python 3 module installation

- **Code Fixes**:
  - Applied automated fixes for print statements, exception syntax
  - Updated string handling in core utility modules
  - Fixed file I/O operations with encoding specifications

## Next Steps

### 1. Continue Automated Migration

```bash
# Run the enhanced migration script on remaining Python files
python3 automate_py3_updates.py --paths "src/python" --apply

# Generate updated report of remaining issues
python3 check_py3_issues.py --generate-report
```

### 2. Test Cython Integration

```bash
# Test with mock implementation
HAPLO_USE_MOCK=1 python3 test_cython_mock.py

# Test with actual Cython modules (if build successful)
python3 test_cython_module_py3.py

# Run comprehensive C++ integration tests
bash test_cpp_integration.sh
```

### 3. Build and Verify

```bash
# Build with Python 3 installer
python3 install_py3.py /tmp/happy-py3-build

# Run the test suite
cd /tmp/happy-py3-build
src/sh/run_tests.sh

# Capture and analyze failures
src/sh/run_tests.sh 2>&1 | tee test_output.log
grep -A 5 "FAILED" test_output.log
```

### 4. Targeted Fixes

Focus on priority modules in this order:
1. Fix remaining Haplo module issues (string/bytes handling)
2. Update Tools module with proper encoding
3. Fix Somatic module integrations
4. Update main scripts (hap.py, som.py, qfy.py)

### 5. Update Documentation

- Update `PYTHON3_MIGRATION.md` with progress
- Document any module-specific considerations
- Update test documentation for Python 3 compatibility

## Common Issues to Watch For

- String vs bytes handling, especially in file I/O
- Unicode encoding/decoding in VCF parsing
- Memory management in Cython modules
- Path handling differences between Python 2 and 3
- Iterator behavior changes

---

_For development workflow details, see `.github/instructions/workflow.instructions.md`_
_For Cython-specific guidance, see `.github/instructions/cython-modernization.instructions.md`_
