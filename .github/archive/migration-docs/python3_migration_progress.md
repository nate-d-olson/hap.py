# hap.py Python 2 to 3 Migration Progress

## Current Status

We have made significant progress towards migrating the hap.py code base from Python 2 to Python 3. The following components have been addressed:

### Completed

1. **Core Python Modules** (Python 3 versions created with `.py3` extension)
   - `src/python/Tools/__init__.py.py3`
   - `src/python/Tools/bcftools.py.py3`
   - `src/python/Tools/fastasize.py.py3`
   - `src/python/Haplo/__init__.py.py3` 
   - `src/python/Haplo/quantify.py.py3`
   - `src/python/pre.py.py3`
   - `src/python/hap.py.py3`
   - `src/python/qfy.py.py3`

2. **Build System Updates**
   - `install_py3.py` (Python 3 compatible version of install.py)
   - `happy.requirements.py3.txt` (Updated dependencies for Python 3)
   - `external/make_dependencies_py3.sh` (Python 3 compatible dependency build script)
   - `CMakeLists.txt.py3` (Updated with modern CMake practices)
   - `external/patches/boost_modern_compiler.patch` (Fix for Boost to work with modern compilers)

3. **Utility Scripts**
   - `convert_py2to3.py` (Automation script for Python 2 to 3 conversion)
   - `diagnose_build.sh` (Diagnostic script to capture build logs and identify issues)
   - `py2to3_migrate.sh` (Orchestration script for the migration process)
   - `test_python3_compatibility.py` (Test script to validate Python 3 compatibility)

### Pending

1. **Remaining Python Modules**
   - ✅ Added Python 3 version of `Tools/bedintervaltree.py`
   - ✅ Added Python 3 version of `Somatic/Mutect.py`
   - Additional modules in `Tools` directory
   - Other modules in `Haplo` package
   - Remaining modules in `Somatic` package

2. **Cython Integration**
   - ✅ Created `CythonSupport.cmake` for Python 3 Cython integration
   - ✅ Implemented `_internal.pyx` Cython module for C++ integration
   - ✅ Added `mock_internal.py` for testing without C++ dependencies
   - ✅ Created proper Python package structure for Cython modules
   - ✅ Updated Haplo.__init__ with robust fallback mechanisms
   - ✅ Added string encoding/decoding handling for Python 3 compatibility
   - ✅ Created test_cython_integration.py for validation
   - ⬜ Build and test with actual C++ code

3. **Build System**
   - Test Python 3 build process end-to-end
   - Ensure all C++ components compile correctly
   - Test with different compilers and platforms

4. **Tests**
   - Update test scripts for Python 3 compatibility
   - Create comprehensive integration tests
   - Compare outputs between Python 2 and Python 3 versions

## Migration Strategy

Our migration approach has been methodical, focusing first on the core modules needed for basic functionality:

1. **Phase 1: Core Module Migration**
   - Convert critical modules to Python 3, preserving original files
   - Add proper type hints and update documentation
   - Create unit tests to validate behavior

2. **Phase 2: Build System Updates**
   - Update CMake configuration for modern practices
   - Modernize dependency handling
   - Fix compiler compatibility issues

3. **Phase 3: Comprehensive Testing**
   - Create a comprehensive test suite for all modules
   - Compare outputs between Python 2 and 3 versions
   - Ensure identical results for all test cases

4. **Phase 4: Full Deployment**
   - Replace Python 2 versions with Python 3 versions
   - Update documentation and examples
   - Create a continuous integration pipeline

## How to Use

### Testing Python 3 Compatibility

To test Python 3 compatibility of the converted modules without modifying the original files:

```bash
python3 test_python3_compatibility.py
```

### Migration

To migrate core modules from Python 2 to Python 3:

```bash
bash py2to3_migrate.sh --core
```

To migrate all modules:

```bash
bash py2to3_migrate.sh --all
```

### Running with Python 3

To run the installation with Python 3:

```bash
python3 install_py3.py [build_dir]
```

## Next Steps

1. Continue converting remaining Python modules to Python 3
2. Address Cython integration issues
3. Run the diagnostic build script to identify build failures
4. Complete a comprehensive test suite
5. Test with real-world datasets to ensure identical results

## Known Issues

1. String handling in Python 3 requires careful consideration of encoding/decoding
2. Some dependencies might not be available in compatible versions
3. Boost and other C++ libraries might need additional patches
4. Cython integration requires specific attention for string handling
