---
applyTo: "**/*.{cmake,CMakeLists.txt}"
---
# Build System Fixes for hap.py

## Diagnosis Approach

1. **Capture Detailed Build Logs**
   - Run builds with verbose output: `python install.py /tmp/happy-build VERBOSE=1 2>&1 | tee build_log.txt`
   - Extract specific errors: `grep -A 10 "error:" build_log.txt > cpp_errors.txt`
   - Identify dependency failures: `grep -A 5 "Could not find" build_log.txt > missing_deps.txt`

2. **Common Build Issues**
   - Missing system dependencies (htslib, etc.)
   - Compiler flag incompatibilities
   - CMake version differences
   - Hardcoded paths in build scripts

## CMake Modernization

1. **Update External Library Integration**
   - Replace manual downloads with CMake FetchContent
   - Move from ExternalProject to modern approaches
   - Add version checking for dependencies

2. **Compiler Configuration**
   - Update C++ standard flags (ensure C++11 compatibility)
   - Use target-specific compile options
   - Add platform-specific compiler detection

3. **Dependency Handling**
   - Centralize version requirements
   - Add fallback strategies for missing system libraries
   - Use CMake's find_package with CONFIG mode when possible

## External Library Fixes

1. **htslib Integration**
   - Update to compatible htslib version
   - Fix include paths and library linking
   - Add option to use system htslib when available

2. **Boost Dependencies**
   - Update Boost usage for modern CMake
   - Address deprecated Boost components
   - Fix header-only vs compiled library usage

3. **zlib and Other Common Libraries**
   - Ensure consistent usage across components
   - Fix static vs dynamic linking issues
   - Address platform-specific differences

## Testing Build Fixes

1. **Incremental Testing**
   - Test external library builds separately: `python install.py --build-externals-only /tmp/happy-ext`
   - Test C++ components separately: `python install.py --build-cpp-only /tmp/happy-cpp`
   - Test full system builds after fixes: `python install.py /tmp/happy-full`

2. **Platform Verification**
   - Test on all target platforms (Linux, macOS)
   - Document platform-specific requirements
   - Create platform-specific build instructions