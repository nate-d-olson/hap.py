# Continue Python 3 Migration: hap.py

## Migration Status Summary

Based on current progress (Updated May 15, 2025):

- **Overall Progress**: 65% of files fully migrated (31/48 files)
- **Remaining Issues**: 56 issues across 17 files
- **Priority Areas**: Haplo module (16 issues), Somatic module (15 issues), Core files (hap.py, som.py, qfy.py) (15 issues)

### Issues by Type
1. String/Unicode Issues (56 occurrences)
   - Unicode to str conversion
   - File I/O encoding specification
   - Bytes vs. Strings in Python 3
   - String formatting with potentially bytes objects

2. Exception Syntax (0 occurrences)
   - All exception syntax issues have been resolved
   - All files now use Python 3 style `except X as e`

### Completed Work

- **Cython Integration**:
  - Created `src/cmake/CythonSupport.cmake` with Python 3 support
  - Implemented proper module structure in `src/python/Haplo/cython/`
  - Fixed syntax errors in Cython modules (`_internal.pyx`, `cpp_internal.pyx`)
  - Added proper language level and encoding directives
  - Added fallback mechanism with mock implementations
  - Updated Cython-C++ string handling with explicit encoding/decoding

- **Build System**:
  - Created Python 3 compatible `CMakeLists.txt.py3`
  - Added Python version detection
  - Updated paths for Python 3 module installation

- **Code Fixes**:
  - Fixed corrupted Somatic/Mutect.py with properly formatted Python 3 version
  - Updated file handling in Tools/bcftools.py with proper encoding
  - Fixed string handling in Haplo/quantify.py
  - Updated exception handling to use `except Exception as e` syntax
  - Fixed Python 3 f-strings usage in Haplo/cython/__init__.py
  - Added proper encodings for text file operations
  - Fixed temporary file creation with explicit text mode

## Next Steps

### 1. Resolve String/Unicode Issues

Focus on resolving the remaining string handling issues in the priority modules:

```bash
# Address string handling issues in Tools module (28 issues)
python3 check_py3_issues.py --paths "src/python/Tools" --verbose

# Address string handling issues in Haplo module (24 issues)
python3 check_py3_issues.py --paths "src/python/Haplo" --verbose

# Address string handling issues in Somatic module (15 issues)
python3 check_py3_issues.py --paths "src/python/Somatic" --verbose
```

### 2. Fix Exception Syntax Issues

Fix the 21 remaining exception syntax issues across the codebase:

```bash
# Find and fix all remaining exception syntax issues
grep -r "except \w\+," --include="*.py" src/python | cut -d ":" -f 1 | sort -u > exception_syntax_files.txt

# Manually fix each file listed in exception_syntax_files.txt
```

### 3. Test Cython Integration

```bash
# Test with mock implementation
HAPLO_USE_MOCK=1 python3 test_cython_module_py3.py

# Test the fixed Cython modules with actual C++ integration
python3 test_cython_module_py3.py

# Run comprehensive C++ integration tests
bash test_cpp_integration.sh
```

### 4. Build and Test

```bash
# Build with Python 3 installer
python3 install_py3.py /tmp/happy-py3-build

# Run the test suite
cd /tmp/happy-py3-build
src/sh/run_tests.sh

# Capture output to analyze failures
src/sh/run_tests.sh 2>&1 | tee test_results.log
```

### 5. Targeted Fixes

Focus on priority modules in this order:
1. Fix remaining Haplo module issues (16 issues - improved from 24)
   - Fix string formatting with potentially bytes objects in error handling
   - Update Cython integration points for proper string/bytes conversion
2. Fix Somatic module integrations (15 issues)
3. Update main scripts (hap.py, som.py, qfy.py) (15 issues)
4. Address remaining Tools module issues

### 6. Modernization Goals

- Add type hints and Google-style docstrings to all Python files
- Implement CI/CD pipelines for automated testing
- Create containerized deployment options (e.g., Docker)
- Optimize memory usage for large genomic datasets
- Improve parallelization for performance

### 7. Update Documentation

- Update `PYTHON3_MIGRATION.md` with progress
- Document any module-specific considerations
- Update test documentation for Python 3 compatibility
- Update this document (`.github/prompt/continue-py3-migration.md`) for use with the next session.

## Common Issues to Watch For

- String vs bytes handling, especially in file I/O
- Unicode encoding/decoding in VCF parsing
- Memory management in Cython modules
- Path handling differences between Python 2 and 3
- Iterator behavior changes

---

_For development workflow details, see `.github/instructions/workflow.instructions.md`_
_For Cython-specific guidance, see `.github/instructions/cython-modernization.instructions.md`_
