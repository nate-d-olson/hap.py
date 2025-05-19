# Continue Python 3 Migration: hap.py

## Migration Status Summary (May 17, 2025)

- **Overall Progress**: 85.7% of files fully migrated (42/49 files)
- **Remaining Issues**: 28 issues across 7 files
- **Priority Areas**:
  - Core functionality optimization (primary focus)
  - Testing and validation of Python 3 code
  - Documentation and deployment

### Recently Completed Tasks

- **Core Functionality Streamlining**:
  - Focused on vcfeval comparison engine and stratified metrics
  - Removed non-core components (som.py, bamstats.py, scmp and xcmp engines)
  - Created feature branches to preserve removed functionality:
    - `feature/alternative-engines` for scmp and xcmp
    - `feature/somatic-support` for somatic variant calling

- **Python 3 Compatibility**:
  - Updated shebang lines from `python` to `python3`
  - Fixed string handling (bytes vs Unicode)
  - Fixed subprocess execution and error handling
  - Improved file I/O with proper encoding
  - Added improved type hints for better code clarity

- **Build & Cython Integration**:
  - Improved NumPy detection in `CythonSupport.cmake`
  - Simplified CMakeLists to remove NumPy from `find_package(Python3 ...)`
  - Excluded `pybedtools` from build requirements to avoid wheel errors

### Next Steps

1. **Testing and Validation (Priority)**:
   - Run comprehensive tests on the updated codebase
   - Create unit tests for core functionality
   - Verify output matches Python 2 version within acceptable tolerance
   - Add test cases for edge cases in string handling

   ```bash
   ./test_py3_core.sh
   python3 -m pytest src/python/Haplo/tests/
   ```

2. **Performance Optimization**:
   - Profile and optimize critical paths
   - Improve memory usage for large VCF files
   - Enhance parallel processing capabilities
   - Optimize Python-C++ interface

   ```bash
   python3 -m cProfile -o happy_profile.prof src/python/hap.py [args]
   python3 -m pstats happy_profile.prof
   ```

3. **Further Code Modernization**:
   - Complete type annotations throughout codebase
   - Use `pathlib` for file operations
   - Implement context managers for resource handling
   - Use f-strings consistently for string formatting

4. **CI/CD Setup**:
   - Create GitHub Actions workflow for automated testing
   - Implement code quality checks in the CI pipeline
   - Add test coverage reporting
   - Set up automatic documentation updates

5. **Documentation Updates**:
   - Expand user documentation with Python 3 specific information
   - Create migration guide for users
   - Update API documentation
   - Add examples for common use cases
