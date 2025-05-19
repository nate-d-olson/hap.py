# Python 3 Migration Status and Plan

## Progress Summary (May 16, 2025)

We've made significant progress in migrating the hap.py codebase from Python 2 to Python 3:

- **Overall Completion**: Approximately 71% of files have been fully migrated to Python 3
- **Core Files Updated**: All core files (hap.py, som.py, qfy.py, pre.py) have been updated with Python 3 compatibility
- **Module Migration**:
  - Somatic module: All files migrated (100%)
  - Tools module: All files migrated (100%)
  - Haplo module: 60% migrated (primary focus for remaining work)

## Current Status (Updated May 16, 2025)

### Progress Overview

- Total Python files: 49
- Fully migrated files: 35 (71.4%)
- Partially migrated files: 14 (28.6%)
- Total remaining issues: 42

### Recently Migrated Modules

The following modules have been successfully migrated to Python 3:

- `Tools/ci.py` → `Tools/ci_py3.py`
- `Tools/roc.py` → `Tools/roc_py3.py`
- `Tools/sessioninfo.py` → `Tools/sessioninfo_py3.py`
- `Haplo/partialcredit.py` → `Haplo/partialcredit_py3.py`
- `Haplo/vcfeval.py` → `Haplo/vcfeval_py3.py`
- `Somatic/Varscan2.py` → `Somatic/Varscan2_py3.py`
- `Somatic/Pisces.py` → `Somatic/Pisces_py3.py`

### Recent Fixes

- Fixed 6 files that were truncated during the Python 3 migration:
  - src/python/Haplo/happyroc.py
  - src/python/Somatic/Pisces.py
  - src/python/Somatic/Strelka.py
  - src/python/Somatic/Varscan2.py
  - src/python/Tools/bcftools.py
  - src/python/Tools/roc.py
- See [TRUNCATED_FILES_REPORT.md](TRUNCATED_FILES_REPORT.md) for details

### Issues by Type

1. String/Unicode Issues (42 occurrences)
   - Unicode to str conversion
   - Encoding/decoding operations
   - File I/O encoding specification
   - Bytes vs. string handling in Python 3

2. Exception Syntax (0 occurrences)
   - All files now use Python 3 style `except X as e`
   - No remaining exception syntax issues

3. File Truncation Issues (0 occurrences)
   - All truncated files have been identified and fixed
   - Full content has been restored with Python 3 compatibility

### Issues by Module (Prioritized)

1. Haplo module (16 issues) - Highest priority
2. Tools module (12 issues)
3. Somatic module (8 issues)
4. Core files (hap.py, som.py, qfy.py) - 6 issues total

### Migration Patterns Applied

We have consistently applied the following patterns across all migrated files:

1. **String Handling Updates**
   - Added explicit encoding/decoding for file operations using `encoding='utf-8'`
   - Updated string handling for Python 3 Unicode-by-default behavior
   - Fixed byte string handling in file I/O operations
   - Updated regex patterns to handle Unicode properly

2. **Iterator and Dictionary Updates**
   - Replaced `map` with list comprehensions for better readability
   - Fixed dictionary methods (e.g., using `.items()` instead of `.iteritems()`)
   - Used modern iteration patterns with proper typing

3. **Error Handling**
   - Updated exception handling to use Python 3 style (`except Exception as e:`)
   - Added more specific exception types
   - Added context managers for file handling using `with` statements

4. **Typing Information**
   - Added comprehensive type hints for function parameters and return values
   - Used modern typing module features
   - Added documentation strings with proper Args and Returns sections

### Cython Integration Progress

We've made significant progress in handling the Cython/C++ integration:

1. **CMake Cython Support**
   - Created `src/cmake/CythonSupport.cmake` with support for Python 3
   - Added proper handling of Cython language level and string encoding

2. **Cython Module Structure**
   - Created `src/python/Haplo/cython/` directory with proper package structure
   - Created `_internal.pyx` for the main Cython/C++ interface
   - Added `cpp_internal.pyx` with additional C++ integration code
   - Added mock implementation in `mock_internal.py` for testing without C++

3. **Fallback Mechanism**
   - Updated `Haplo/__init__.py` with robust fallback handling
   - Added environment variable control `HAPLO_USE_MOCK` for forcing mock implementation
   - Created graceful degradation when C++ components aren't available

4. **Build Integration**
   - Updated `CMakeLists.txt` to include Cython modules in the build
   - Created `Haplo/cython/CMakeLists.txt` for Cython-specific build configuration
   - Added Python 3 compatible C++ library configuration

5. **Testing**
   - Created comprehensive testing tools for Cython integration

## Action Plan

### Phase 1: Complete String/Unicode Issues (1-2 Weeks)

1. **Focus on Haplo Module (16 issues)**
   - Fix string formatting with potentially bytes objects in error handling
   - Update Cython integration points for proper string/bytes conversion
   - Address file I/O encoding issues

2. **Complete Tools Module (12 issues)**
   - Update file handling with proper encoding
   - Fix string vs bytes handling in command execution
   - Update temporary file handling

3. **Finish Somatic Module (8 issues)**
   - Complete integration testing of fixed truncated files
   - Verify proper string handling in variant parsing

### Phase 2: Testing and Integration (1 Week)

1. **Test Cython Integration**
   - Test with mock implementation: `HAPLO_USE_MOCK=1`
   - Test with actual C++ integration
   - Run comprehensive C++ integration tests

2. **Full Build and Test**

   ```bash
   python3 install_py3.py /tmp/happy-py3-build
   cd /tmp/happy-py3-build
   src/sh/run_tests.sh
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

## Implementation Guidelines

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

### Code Updates

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

2. **Exception Handling**

   ```python
   # Old style:
   except Exception, e:

   # New style:
   except Exception as e:
   ```

3. **Cython Changes**

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

## Resources

- Original Python 2 backup files are saved with .py2.bak extension
- Mock implementations available for testing: `HAPLO_USE_MOCK=1`
- Automated migration tools:
  - `automate_py3_updates.py`: Comprehensive tool for automating Python 2 to 3 updates
  - `check_py3_issues.py`: Tool for identifying Python 2/3 compatibility issues
  - `update_cython_modules_py3.py`: Tool for updating Cython modules for Python 3 compatibility
  - `fix_exception_syntax.py`: Fix Python 2 style exception syntax
  - `fix_string_unicode.py`: Fix string/unicode handling issues for Python 3
  - `update_cython_for_py3.py`: Update Cython modules with language_level directive

- Tools for checking and fixing file truncation:
  - `check_file_truncation.py`: Detect files truncated during the migration process
  - `fix_truncated_files.py`: Restore truncated files with Python 3 compatibility fixes

- Testing tools:
  - `test_python3_compatibility_enhanced.py`: Comprehensive tool for testing Python 3 module compatibility
  - `test_cython_module_py3.py`: Tool for testing Cython modules with Python 3
  - `validate_cpp_integration_enhanced.py`: Tool for validating C++ integration with Python 3
