# Python 3 Migration Action Plan

## Phase 0: Setup and Assessment (1 Week)

1. **Environment Setup**
   ```bash
   # Create Python 3 virtual environment
   python3 -m venv venv_py3
   source venv_py3/bin/activate
   pip install -r requirements-dev.txt
   ```

2. **Codebase Analysis**
   ```bash
   # Run 2to3 to identify incompatibilities
   2to3 -p src/python > py3_conversion_report.txt
   
   # Identify critical modules with most issues
   grep -E "RefactoringTool|@@" py3_conversion_report.txt | sort -k2 -r | head -20
   ```

3. **Dependency Analysis**
   ```bash
   # Verify Python 3 compatibility of dependencies
   pip install caniusepython3
   caniusepython3 -r requirements.txt
   ```

4. **Test Infrastructure Baseline**
   ```bash
   # Run existing tests with Python 2 to establish a baseline
   cd /tmp/happy-py2-baseline
   src/sh/run_tests.sh 2>&1 | tee py2_test_baseline.log
   
   # Analyze test coverage
   grep -A 3 "PASS\|FAIL\|ERROR" py2_test_baseline.log > test_status_summary.txt
   
   # Catalog test data being used
   find src/data -type f | sort > test_data_inventory.txt
   ```

## Pre-Implementation Validation

1. **Dependency Management**
   ```bash
   # Create separate requirements files for Python 3
   pip-compile --upgrade --output-file=requirements-py3.txt requirements.in
   pip-compile --upgrade --output-file=requirements-dev-py3.txt requirements-dev.in
   ```

2. **Genomic Output Validation Plan**
   - Leverage existing test datasets in `src/data` directory for validation
   - Define acceptable differences in numerical outputs (floating point tolerance)
   - Create automated comparison scripts for running existing tests with both Python 2 and 3
   ```bash
   # Example validation command using existing test data
   python validate_genomic_results.py --test-dir src/sh --data-dir src/data --py2-bin /tmp/happy-py2-baseline/bin/hap.py --py3-bin /tmp/happy-py3-build/bin/hap.py --tolerance 1e-10
   ```

3. **Migration Staging Plan**
   - Stage 1: Create Python 3 compatible branches without breaking Python 2
   - Stage 2: Fork main repository for Python 3-only development
   - Stage 3: Final transition with version number bump (e.g., 0.x → 1.0)
   
4. **Parallel Processing Review**
   - Audit multiprocessing code for Python 3 compatibility
   - Test process pool implementations with various worker counts
   - Verify GIL handling in Cython components

## Phase 1: Cython Module Updates (2 Weeks)

1. **Analyze and Update Core Cython Modules**
   ```bash
   # Run the update script on all Cython modules
   python update_cython_modules_py3.py --src-dir src/python
   
   # Test updated modules
   python test_cython_module_py3.py --build-dir /tmp/happy-py3-test
   ```

2. **Verify String Handling**
   - Test Unicode character handling in variant IDs
   - Test file path handling with non-ASCII characters
   - Update string conversion in C++ interfaces
   - Create tests for bytes vs string conversions
   - Update file I/O operations for proper text/binary mode

3. **Fix Memory Management Issues**
   - Review all PyObject* handling
   - Ensure proper reference counting
   - Update memory allocation patterns
   - Test with memory profiler to detect leaks

## Phase 2: Core Python Updates (2 Weeks)

1. **Update Standard Library Imports**
   ```bash
   # Update common imports
   find src/python -name "*.py" | xargs sed -i 's/import Queue/import queue as Queue/g'
   find src/python -name "*.py" | xargs sed -i 's/import ConfigParser/import configparser as ConfigParser/g'
   ```

2. **Fix Common Python 3 Issues**
   - Update dict methods (items() vs iteritems())
   - Fix integer division (/ vs //)
   - Update exception syntax (except Exception, e → except Exception as e)
   - Convert print statements to functions
   - Update relative imports
   - Fix comparisons and sorting functions

3. **Add Type Hints to Core Modules**
   ```bash
   # Install type checking tools
   pip install mypy
   
   # Add type hints to critical modules
   mypy --install-types src/python/Haplo/
   ```

## Phase 3: Build System Integration (2 Weeks)

1. **Update CMake Configuration**
   ```bash
   # Update CMake files for Python 3
   cp src/cmake/CythonSupport.py3.cmake src/cmake/CythonSupport.cmake
   # Update main CMakeLists.txt
   ```

2. **Create Modern Package Structure**
   - Add pyproject.toml
   - Update setup.py for Python 3
   - Create proper package namespace structure

3. **Test Platform-Specific Builds**
   ```bash
   # Test on macOS
   ./test_py3_build.sh
   
   # Test on Linux (in Docker)
   docker run -it --rm -v $(pwd):/src ubuntu:20.04 /src/test_py3_build.sh
   ```

4. **Update External Dependencies**
   ```bash
   # Test updated external dependencies script
   ./external/make_dependencies_py3.sh rebuild
   ```

## Phase 4: Testing and Integration (2 Weeks)

1. **Connect Migrated Modules**
   - Update `__init__.py` files
   - Test imports across modules
   - Fix remaining inter-module dependencies

2. **Leverage Existing Test Infrastructure**
   ```bash
   # Run the existing test suite with Python 3
   cd /tmp/happy-py3-build
   src/sh/run_tests.sh 2>&1 | tee py3_test_run.log
   
   # Compare test results between Python 2 and 3
   python compare_test_results.py --py2-log /tmp/happy-py2-baseline/py2_test_baseline.log --py3-log py3_test_run.log
   
   # Debug failing tests individually
   for test in $(grep -l "FAIL\|ERROR" /tmp/happy-py3-build/src/sh/*.log); do
     echo "Debugging test: $test"
     # Run individual test with more verbosity
     bash $test --debug 2>&1 | tee debug_$(basename $test).log
   done
   ```

3. **Test Data Validation**
   ```bash
   # Verify outputs with test data from src/data directory
   python validate_test_outputs.py --data-dir src/data --output-dir /tmp/happy-py3-build/test_outputs
   
   # Test with specific benchmark datasets
   src/sh/run_benchmark_tests.sh src/data/chr21.1M.vcf.gz 2>&1 | tee benchmark_results.log
   
   # Compare output formats and values
   python check_output_compatibility.py --py2-output /tmp/happy-py2-baseline/test_outputs --py3-output /tmp/happy-py3-build/test_outputs
   ```

4. **Performance Testing**
   ```bash
   # Run performance comparisons with test data
   python compare_performance_py2_py3.py --vcf src/data/PG_performance.vcf.gz
   ```

5. **Validate Real Data Processing**
   ```bash
   # Test with example data using both runtimes
   time /tmp/happy-py2-baseline/bin/hap.py src/data/test/truth.vcf.gz src/data/test/query.vcf.gz -o /tmp/py2_output
   time /tmp/happy-py3-build/bin/hap.py src/data/test/truth.vcf.gz src/data/test/query.vcf.gz -o /tmp/py3_output
   
   # Compare outputs between Python 2 and 3 versions
   python compare_outputs.py --py2-out /tmp/py2_output --py3-out /tmp/py3_output --tolerance 1e-10
   ```

## Phase 5: Finalization and Documentation (1 Week)

1. **Create Unified Installer**
   ```bash
   # Update install.py with Python version detection and support
   python update_installer.py
   
   # Test installation with Python 3
   python3 install.py /tmp/happy-py3-build
   ```

2. **Update Documentation**
   - Update README.md with Python 3 compatibility information
   - Create migration guide for users
   - Update API documentation with new types
   - Create examples for common usage patterns
   - Document any breaking changes

3. **Final Validation**
   ```bash
   # Run full test suite with the established test infrastructure
   cd /tmp/happy-py3-build
   src/sh/run_tests.sh
   
   # Create test report
   src/sh/generate_test_report.sh > test_report.txt
   
   # Test with real-world datasets from src/data
   find src/data -name "*.vcf.gz" -type f | sort | xargs -I{} \
     bash -c 'echo "Testing {}" && /tmp/happy-py3-build/bin/hap.py {} src/data/test/query.vcf.gz -o test_output_$(basename {})'
   ```

## Immediate Next Steps

1. Set up Python 3 virtual environment for development
2. Run 2to3 analysis on the codebase to identify hotspots
3. Verify all testing scripts in src/sh directory are compatible with Python 3
   ```bash
   # Check syntax compatibility of test scripts
   find src/sh -name "*.py" | xargs 2to3 -p > test_scripts_py3_compatibility.txt
   
   # Ensure test data access is compatible
   python verify_test_data_access.py --data-dir src/data
   ```
4. Create unit tests for critical components before modification
5. Apply Cython updates to the core modules:
   ```
   src/python/Haplo/cython/_internal.pyx
   src/python/Haplo/cython/cpp_internal.pyx
   src/python/Haplo/hapcompare.pyx
   ```

6. Create a prototype build with Python 3:
   ```bash
   python install_py3.py /tmp/happy-py3-build
   ```

## Key Success Criteria

1. All tests pass with Python 3
2. Performance is comparable to Python 2 version (within 10%)
3. Memory usage is stable or improved
4. All features work as expected
5. Documentation is updated and accurate
6. Installation is smooth on all supported platforms
7. Type checking passes with minimal warnings

## Risk Mitigation

1. Keep Python 2 version working throughout the migration
2. Create mock implementations for all C++ components
3. Add fallback mechanisms for critical functionality
4. Use extensive logging for debugging
5. Test on multiple platforms
6. Set up CI/CD for Python 3 builds
7. Create automated conversion validation tests
8. Document all changes and their impacts
9. Create a test data compatibility layer to handle any format changes between Python 2 and 3
