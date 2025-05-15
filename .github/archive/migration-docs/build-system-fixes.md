# Build System Fixes for hap.py Python 3 Migration

## Identified Issues

Based on our analysis of the build system, we've identified several potential issues that need to be addressed:

1. **External C++ Dependencies**
   - Build process uses outdated dependency versions (zlib 1.2.8, boost 1.58.0)
   - Hardcoded download URLs and paths in `external/make_dependencies.sh`
   - No fallback to system libraries when available

2. **Python 2 Dependencies**
   - Virtualenv setup is fixed to version 12.0.7 (not Python 3 compatible)
   - `urllib2` usage in `install.py` (Python 2 specific)
   - Python 2 specific print statements

3. **CMake Configuration**
   - CMake minimum version 2.8 (outdated)
   - Potential deprecated CMake commands
   - Boost configuration may not be compatible with modern Boost

4. **Cython Integration**
   - No language_level specification for Cython
   - Potential string encoding/decoding issues

## Fix Implementation Plan

### 1. Update Build Environment

1. **Update `install.py` for Python 3 compatibility**
   - Replace Python 2 print statements
   - Modernize Python code
   - Update virtualenv creation process

2. **Add system library detection**
   - Add options to use system libraries (boost, zlib, etc.)
   - Implement version checking for dependencies

3. **Update dependency versions**
   - Update to modern versions of dependencies
   - Add SHA256 verification for downloads

### 2. Update CMake Configuration

1. **Modernize CMake files**
   - Update minimum required version
   - Replace deprecated commands
   - Use modern CMake targets

2. **Fix Boost configuration**
   - Update Boost finding and configuration
   - Add version compatibility checks

3. **Improve error reporting**
   - Add more informative error messages
   - Add dependency checks

### 3. Update C++ Code

1. **Fix compilation warnings**
   - Address deprecated API usage
   - Fix standard compliance issues

2. **Ensure C++11 compatibility**
   - Update code for C++11 standard
   - Use modern C++ features where appropriate

### 4. Update Python/Cython Integration

1. **Update Cython code**
   - Add language_level=3 directive
   - Fix string encoding/decoding
   - Update memory management

2. **Update Python bindings**
   - Ensure compatibility with Python 3
   - Fix string handling

## Implementation Order

1. **Phase 1: Build System Fixes**
   - Update `install.py` for Python 3
   - Fix CMake configuration
   - Update external dependency handling

2. **Phase 2: Core Python Module Updates**
   - Update Tools modules
   - Update pre.py
   - Update qfy.py

3. **Phase 3: Cython Module Updates**
   - Update Cython directives
   - Fix string handling
   - Update memory management

4. **Phase 4: Main Script Updates**
   - Update hap.py main script
   - Update command-line argument handling

5. **Phase 5: Testing and Integration**
   - Run test suite
   - Fix test failures
   - Document changes
