# Continue Python 3 Migration: hap.py

## Migration Status Summary

Based on current progress (Updated May 16, 2025):

- **Overall Progress**: 71.4% of files fully migrated (35/49 files)
- **Remaining Issues**: 42 issues across 14 files
- **Priority Areas**: Haplo module (16 issues), Tools module (12 issues), Somatic module (8 issues), Core files (hap.py, som.py, qfy.py) (6 issues)

### Recently Migrated Modules

The following modules have been successfully migrated to Python 3:

- `Tools/ci.py` → `Tools/ci_py3.py`
- `Tools/roc.py` → `Tools/roc_py3.py`
- `Tools/sessioninfo.py` → `Tools/sessioninfo_py3.py`
- `Haplo/partialcredit.py` → `Haplo/partialcredit_py3.py`
- `Haplo/vcfeval.py` → `Haplo/vcfeval_py3.py`
- `Somatic/Varscan2.py` → `Somatic/Varscan2_py3.py`
- `Somatic/Pisces.py` → `Somatic/Pisces_py3.py`

### Issues by Type

1. String/Unicode Issues (42 occurrences)
   - Unicode to str conversion
   - File I/O encoding specification
   - Bytes vs. Strings in Python 3
   - String formatting with potentially bytes objects
   - Encoding/decoding operations

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
  - Created `Haplo/cython/CMakeLists.txt` for Cython-specific build configuration

- **Build System**:
  - Created Python 3 compatible `CMakeLists.txt.py3`
  - Added Python version detection
  - Updated paths for Python 3 module installation

- **Code Fixes**:
  - Fixed corrupted Somatic/Mutect.py with properly formatted Python 3 version
  - Updated file handling in Tools/bcftools.py with proper encoding
  - Fixed string handling in Haplo/quantify.py
  - Updated exception handling to use `except Exception as e` syntax
  - Fixed Python 3 f-strings usage in Haplo/cython/**init**.py
  - Added proper encodings for text file operations
  - Fixed temporary file creation with explicit text mode

## Structured Action Plan

### Phase 1: Complete String/Unicode Issues (1-2 Weeks)

1. **Focus on Haplo Module (16 issues)**
   - Fix string formatting with potentially bytes objects in error handling
   - Update Cython integration points for proper string/bytes conversion
   - Address file I/O encoding issues

   ```bash
   # Address string handling issues in Haplo module
   python3 check_py3_issues.py --paths "src/python/Haplo" --verbose
   ```

2. **Complete Tools Module (12 issues)**
   - Update file handling with proper encoding
   - Fix string vs bytes handling in command execution
   - Update temporary file handling

   ```bash
   # Address string handling issues in Tools module
   python3 check_py3_issues.py --paths "src/python/Tools" --verbose
   ```

3. **Finish Somatic Module (8 issues)**
   - Complete integration testing of fixed truncated files
   - Verify proper string handling in variant parsing

   ```bash
   # Address string handling issues in Somatic module
   python3 check_py3_issues.py --paths "src/python/Somatic" --verbose
   ```

4. **Update Core Files (6 issues)**
   - Ensure command line argument handling is Python 3 compatible
   - Verify proper encoding for all file I/O

### Phase 2: Testing and Integration (1 Week)

1. **Test Cython Integration**

   ```bash
   # Test with mock implementation
   HAPLO_USE_MOCK=1 python3 test_cython_module_py3.py

   # Test the fixed Cython modules with actual C++ integration
   python3 test_cython_module_py3.py

   # Run comprehensive C++ integration tests
   bash test_cpp_integration.sh
   ```

2. **Full Build and Test**

   ```bash
   # Build with Python 3 installer
   python3 install_py3.py /tmp/happy-py3-build

   # Run the test suite
   cd /tmp/happy-py3-build
   src/sh/run_tests.sh

   # Capture output to analyze failures
   src/sh/run_tests.sh 2>&1 | tee test_results.log
   ```

3. **Validate Output Consistency**
   - Compare outputs between Python 2 and Python 3 versions
   - Verify numerical results match within tolerance
   - Test with real genomic datasets

### Phase 3: Code Quality and Modernization (2-3 Weeks)

1. **Add Type Hints and Documentation**
   - Complete type annotations for all Python files
   - Add Google-style docstrings to all functions
   - Update module-level documentation

2. **Implement Modern Python Practices**
   - Use f-strings for string formatting
   - Implement context managers for resource handling
   - Use pathlib for file path operations
   - Leverage modern Python features (dataclasses, etc.)

3. **Code Structure Improvements**
   - Implement proper package structure
   - Create clean public APIs
   - Separate concerns (I/O, processing, analysis)

### Phase 4: Build System and Deployment (1-2 Weeks)

1. **Modernize Build System**
   - Update CMake configuration for modern practices
   - Implement proper package installation
   - Create reproducible builds

2. **CI/CD Integration**
   - Set up automated testing
   - Implement test coverage reporting
   - Create deployment pipelines

3. **Documentation Updates**
   - Update user documentation
   - Create migration guides for users
   - Document new features and APIs

## Tools and Resources

### Migration Tools

- **Automated Migration Tools**:
  - `automate_py3_updates.py`: Comprehensive tool for automating Python 2 to 3 updates
  - `check_py3_issues.py`: Tool for identifying Python 2/3 compatibility issues
  - `update_cython_modules_py3.py`: Tool for updating Cython modules for Python 3 compatibility
  - `fix_exception_syntax.py`: Fix Python 2 style exception syntax
  - `fix_string_unicode.py`: Fix string/unicode handling issues for Python 3
  - `update_cython_for_py3.py`: Update Cython modules with language_level directive

- **File Truncation Tools**:
  - `check_file_truncation.py`: Detect files truncated during the migration process
  - `fix_truncated_files.py`: Restore truncated files with Python 3 compatibility fixes

- **Testing Tools**:
  - `test_python3_compatibility_enhanced.py`: Comprehensive tool for testing Python 3 module compatibility
  - `test_cython_module_py3.py`: Tool for testing Cython modules with Python 3
  - `validate_cpp_integration_enhanced.py`: Tool for validating C++ integration with Python 3

### Using Pre-commit Hooks for Migration

Pre-commit hooks have been added to the project to automate and standardize the Python 3 migration:

1. **Installation**:

   ```bash
   pip install pre-commit
   pre-commit install
   ```

2. **Recommended Migration Workflow**:

   ```bash
   # After initial 2to3 conversion or manual edits
   pre-commit run --files src/python/module/file.py

   # For batch processing
   pre-commit run --files src/python/module/*.py

   # After fixes, verify code quality
   pre-commit run --all-files
   ```

3. **Most Useful Hooks for Migration**:
   - `pyupgrade`: Automatically upgrades Python 2 syntax to Python 3
   - `ruff`: Detects and fixes compatibility issues
   - `black`: Ensures consistent code style
   - `isort`: Properly orders imports for Python 3
   - `mypy`: Optional type checking

4. **Specific Migration Tasks**:
   - Fix string/bytes issues: `ruff --select=B --fix src/python/path/to/file.py`
   - Format code: `black src/python/path/to/file.py`
   - Sort imports: `isort src/python/path/to/file.py`

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

### Code Pattern Updates

1. **String Handling**

   ```python
   # Add to the top of each file:
   from __future__ import unicode_literals

   # Update file operations:
   with open(filename, 'rt', encoding='utf-8') as f:
       content = f.read()

   # Binary file handling:
   with open(filename, 'rb') as f:
       binary_data = f.read()  # Returns bytes in Python 3
   ```

2. **Cython Changes**

   ```cython
   # Add to all .pyx files:
   # cython: language_level=3
   # distutils: language=c++

   # String handling:
   cdef bytes py_bytes = python_str.encode('utf8')
   cdef char* c_str = py_bytes

   # C string to Python string:
   py_str = c_str.decode('utf8') if c_str != NULL else ""
   ```

### Verification Steps

1. Run automated checks:

   ```bash
   python3 check_py3_issues.py --paths src/python --generate-report
   ```

2. Build verification:

   ```bash
   python3 install_py3.py /tmp/happy-py3-build
   ```

3. Test execution:

   ```bash
   cd /tmp/happy-py3-build
   src/sh/run_tests.sh 2>&1 | tee test_results.log
   ```

4. Compare with Python 2 version:

   ```bash
   python compare_outputs.py --py2-out /tmp/py2-output --py3-out /tmp/py3-output
   ```

---

_For development workflow details, see `.github/instructions/workflow.instructions.md`_
_For Cython-specific guidance, see `.github/instructions/cython-modernization.instructions.md`_
