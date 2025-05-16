# Continue Python 3 Migration: hap.py

## Migration Status Summary

Based on current progress (Updated May 16, 2025):

- **Overall Progress**: 71% of files fully migrated (35/49 files)
- **Remaining Issues**: 42 issues across 14 files
- **Priority Areas**: Haplo module (16 issues), Tools module (12 issues), Somatic module (8 issues), Core files (hap.py, som.py, qfy.py) (6 issues)

### Issues by Type
1. String/Unicode Issues (42 occurrences)
   - Unicode to str conversion
   - File I/O encoding specification
   - Bytes vs. Strings in Python 3
   - String formatting with potentially bytes objects

2. Exception Syntax (0 occurrences)
   - All exception syntax issues have been resolved
   - All files now use Python 3 style `except X as e`

3. File Truncation Issues (0 occurrences)
   - All truncated files have been identified and fixed
   - Full content has been restored with Python 3 compatibility

### Completed Work

- **Truncated Files Recovery**:
  - Identified and fixed 6 files truncated during Python 3 migration
  - Restored full content in Haplo/happyroc.py, Somatic/Pisces.py, Somatic/Strelka.py,
    Somatic/Varscan2.py, Tools/bcftools.py, and Tools/roc.py
  - Applied Python 3 compatibility fixes during restoration
  - Created tools to detect and fix file truncation (`check_file_truncation.py` and `fix_truncated_files.py`)
  - Documented full recovery process in TRUNCATED_FILES_REPORT.md

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

### 1. Resolve Remaining String/Unicode Issues

Focus on resolving the remaining string handling issues in the priority modules:

```bash
# Address string handling issues in Haplo module (16 issues)
python3 check_py3_issues.py --paths "src/python/Haplo" --verbose

# Address string handling issues in Tools module (12 issues)
python3 check_py3_issues.py --paths "src/python/Tools" --verbose

# Address string handling issues in Somatic module (8 issues)
python3 check_py3_issues.py --paths "src/python/Somatic" --verbose
```

### 2. Verify Integration of Fixed Truncated Files

Run integration tests to verify the fixed truncated files work correctly:

```bash
# Test truncated files integration
python3 test_fixed_truncated_files.py
```

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

1. Fix remaining Haplo module issues (16 issues)
   - Fix string formatting with potentially bytes objects in error handling
   - Update Cython integration points for proper string/bytes conversion
   
2. Fix remaining Tools module issues (12 issues)
   - Update file handling with proper encoding
   - Fix string vs bytes handling

3. Address Somatic module issues (8 issues)
   - Complete integration testing of fixed truncated files 
   - Verify proper string handling

4. Update main scripts (hap.py, som.py, qfy.py) (6 issues)
   - Ensure command line argument handling is Python 3 compatible
   - Verify proper encoding for all file I/O

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
- File truncation during automated conversion (now fixed, but watch for similar issues)
- Incomplete module dependencies after restoration of truncated files

## Tools Added During Migration

- `check_file_truncation.py` - Detect files truncated during the migration process
- `fix_truncated_files.py` - Restore truncated files with Python 3 compatibility fixes
- `update_cython_for_py3.py` - Update Cython modules with language_level directive and string handling

## Troubleshooting Recommendations

### Handling File Truncation Issues

If more truncated files are discovered:

1. Use the `check_file_truncation.py` script to systematically compare line counts:
   ```bash
   python3 check_file_truncation.py --dir src/python --threshold 0.9
   ```

2. Use `fix_truncated_files.py` to restore truncated files:
   ```bash
   python3 fix_truncated_files.py --files file1.py file2.py
   ```

3. Manually verify the fixed files compile with Python 3:
   ```bash
   python3 -m py_compile file1.py file2.py
   ```

4. Test functionality of restored files with appropriate test cases

### Testing the Fixed Files

The recently fixed truncated files are particularly important to test thoroughly:

1. **Varscan2.py** - Test somatic variant calling functionality
2. **Strelka.py** - Verify variant parsing and integration
3. **bcftools.py** - Ensure VCF handling and manipulation works
4. **roc.py** - Validate ROC curve calculation and plotting
5. **happyroc.py** - Test ROC functionality specific to diploid calls
6. **Pisces.py** - Verify somatic variant calling integration

See [TRUNCATED_FILES_REPORT.md](TRUNCATED_FILES_REPORT.md) for detailed information about the recovery process.

---

_For development workflow details, see `.github/instructions/workflow.instructions.md`_
_For Cython-specific guidance, see `.github/instructions/cython-modernization.instructions.md`_
