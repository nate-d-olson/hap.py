# Python 3 Migration Status and Plan

## Current Status (Updated May 15, 2025)

### Progress Overview
- Total Python files: 49
- Fully migrated files: 35 (71.4%)
- Partially migrated files: 14 (28.6%)
- Total remaining issues: 42

### Issues by Type
1. String/Unicode Issues (56 occurrences)
   - Unicode to str conversion
   - Encoding/decoding operations
   - File I/O encoding specification
   - Bytes vs. string handling in Python 3

2. Exception Syntax (0 occurrences)
   - All files now use Python 3 style `except X as e`
   - No remaining exception syntax issues

### Issues by Module (Prioritized)
1. Haplo (16 issues) - Highest priority
2. Tools (25 issues)
3. Somatic (12 issues)
4. Core files (hap.py, som.py, qfy.py) - 13 issues total
5. Other utilities (ftx.py, cnx.py, etc.) - 8 issues total

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
   - Updated `Haplo/__init__.py.py3` with robust fallback handling
   - Added environment variable control `HAPLO_USE_MOCK` for forcing mock implementation
   - Created graceful degradation when C++ components aren't available

4. **Build Integration**
   - Updated `CMakeLists.txt.py3` to include Cython modules in the build
   - Created `Haplo/cython/CMakeLists.txt` for Cython-specific build configuration
   - Added Python 3 compatible C++ library configuration in `src/c++/lib/CMakeLists.txt.py3`

5. **Testing**
   - Created `test_cython_integration.py` for testing both mock and real implementations
   - Updated `py2to3_migrate.sh` with Cython module setup

## Action Plan

### Phase 1: Critical Fixes (Immediate)

1. **Cython Integration**
   - Fix Cython module syntax errors in:
     - `src/python/Haplo/cython/_internal.pyx`
     - `src/python/Haplo/cython/cpp_internal.pyx`
   - Add Python 3 language level directive to all .pyx files:
     ```cython
     # cython: language_level=3
     ```

2. **Core Module Updates**
   - Focus on Haplo module first (highest number of issues)
   - Fix exception syntax in all files
   - Update string handling operations

### Phase 2: Build System Updates

1. **CMake Configuration**
   - Update Python path handling in CMakeLists.txt
   - Remove hardcoded Python 2.7 references
   - Add Python 3 version detection

2. **Dependencies**
   - Test with updated dependencies from happy.requirements.py3.txt
   - Verify Cython module compatibility

### Phase 3: Testing Infrastructure

1. **Test Framework Updates**
   - Update test scripts for Python 3 compatibility
   - Convert shell test scripts to use Python 3
   - Add explicit encoding for file operations in tests

2. **Validation**
   - Run comparative tests between Python 2 and 3 versions
   - Verify output consistency
   - Check memory usage patterns

## Implementation Guidelines

### Code Updates
1. **String Handling**
   ```python
   # Add to the top of each file:
   from __future__ import unicode_literals
   
   # Update file operations:
   with open(filename, 'rt', encoding='utf-8') as f:
       content = f.read()
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
   ```

### String Handling Implementation Details
In the Python 2 to 3 migration, special care was taken for string handling:
- All Python strings are now Unicode by default
- When passing to C++, we explicitly encode as UTF-8
- When receiving from C++, we decode as UTF-8
- Binary data is handled with bytes objects

### Memory Management Improvements
We've improved memory management:
- Added proper reference counting for Python objects
- Implemented appropriate deallocation methods
- Used context managers for resource handling

### Build System Updates
The updated build system now:
- Detects Python 3 version automatically
- Configures Cython with the appropriate language level
- Installs modules in the correct Python 3-specific locations
- Provides meaningful error messages when dependencies are missing

### Testing
1. Use pytest for new tests
2. Add encoding checks to file I/O tests
3. Verify outputs match between versions

## Verification Steps

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
   bash test_happy_py3.sh
   ```

## Resources

- Original Python 2 backup files are saved with .py2.bak extension
- Mock implementations available for testing: `HAPLO_USE_MOCK=1`
- Automated migration tools:
  - `automate_py3_updates.py`
  - `check_py3_issues.py`
  - `update_cython_modules_py3.py`
  - `fix_exception_syntax.py` - Fix Python 2 style exception syntax (except X, y:)
  - `fix_string_unicode.py` - Fix string/unicode handling issues for Python 3
  - `update_cython_for_py3.py` - Update Cython modules with language_level directive and string handling
