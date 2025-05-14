

## 4. Update 

Add a section specifically on build system modernization:

```markdown
---
applyTo: "**"
---
# Modernization Strategy for hap.py

## Build System Modernization (Priority)

- Fix external library build process
  - Update CMake files for better dependency management
  - Replace outdated build commands with modern equivalents
  - Fix platform-specific issues (Mac vs Linux differences)
  - Document system requirements clearly

- Address C++ compilation issues
  - Update compiler flags for modern compilers
  - Fix deprecated code patterns
  - Ensure C++11 compatibility across codebase
  - Address warnings that cause build failures

- Create reproducible builds
  - Pin dependency versions where appropriate
  - Document exact dependency requirements
  - Enable incremental builds for faster development

- Improve error reporting
  - Add better error messages for common build failures
  - Create troubleshooting guide for build issues
  - Add validation checks before building

## Python 2 to 3 Migration

## Python 2 to 3 Migration
- Replace print statements with print() function
  - Use the `print("text", var)` syntax instead of `print "text", var`
  - Use parentheses for all print statements
  - Consider using f-strings for complex formatting: `print(f"Value: {var}")`
- Update string handling (unicode vs bytes)
  - Convert all `unicode()` to `str()`
  - Replace `str()` with `bytes()` when dealing with binary data
  - Use explicit encoding/decoding: `.encode('utf-8')` and `.decode('utf-8')`
  - Update string literals with b'' prefix for bytes
- Replace xrange with range
  - Remove all instances of `xrange()` and replace with `range()`
  - Consider efficiency implications for very large ranges
- Update dictionary methods
  - Replace `.iteritems()` with `.items()`
  - Replace `.iterkeys()` with `.keys()`
  - Replace `.itervalues()` with `.values()`
  - Update `.has_key()` with `in` operator
  - Note that dict views are now iterator objects, not lists
- Fix division behavior (/ vs //)
  - Update division operations to use `//` for integer division
  - Use `/` for floating point division
  - Add explicit casts where needed: `int(a/b)` or `float(a)`
- Update exception syntax
  - Replace `except Exception, e:` with `except Exception as e:`
  - Update multiple exception handling: `except (TypeError, ValueError) as e:`
- Replace obsolete modules
  - `ConfigParser` → `configparser`
  - `Queue` → `queue`
  - `StringIO`/`cStringIO` → `io.StringIO` or `io.BytesIO`
  - `urllib2` → `urllib.request`
  - `urlparse` → `urllib.parse`
  - `thread` → `_thread`
  - `dbm.*` modules → unified `dbm` interface

## Dependency Updates
- Replace deprecated library functions with modern equivalents
- Update minimum versions of dependencies
- Remove unnecessary dependencies
- Use pip/setuptools modern practices

## Code Architecture
- Convert ad-hoc scripts to well-organized modules
  - Refactor utility functions into dedicated modules
  - Create explicit imports instead of relative imports
  - Use `__init__.py` files to define public APIs
- Implement proper package structure
  - Follow standard Python package layout (src-based)
  - Use namespace packages where appropriate
  - Separate code from configuration and data
  - Implement entry points for command-line scripts
- Separate concerns (data processing, analysis, reporting)
  - Create dedicated modules for I/O operations
  - Separate algorithm implementation from presentation
  - Factor out reusable components from application logic
  - Implement configuration handling separate from business logic
- Create clear APIs between components
  - Define explicit interfaces for component interaction
  - Use type hints for function signatures
  - Document public APIs with docstrings
  - Create integration tests for component boundaries
  - Minimize global state and side effects

## Testing Infrastructure
- Convert bash test scripts to Python test framework
  - Replace shell scripts with pytest fixtures and tests
  - Organize tests by component and functionality
  - Implement proper setup/teardown procedures
  - Create helper utilities for test data generation
  - Add parameterized tests for edge cases
- Implement automated testing with pytest
  - Utilize pytest fixtures for test resources
  - Add pytest.ini configuration file
  - Create conftest.py files for shared fixtures
  - Implement test markers for categorization
  - Add parallel test execution support
- Add CI/CD pipeline integration
  - Set up GitHub Actions or other CI system
  - Create test matrices for different platforms
  - Configure automatic test execution on PR
  - Implement staged testing (unit → integration → system)
  - Add build artifact generation
- Track test coverage
  - Configure pytest-cov for coverage reporting
  - Set minimum coverage thresholds
  - Create coverage reports by component
  - Exclude appropriate files from coverage metrics
  - Integrate coverage reporting with CI

## Build System
- Modernize CMake usage
  - Use target-based approach instead of global variables
  - Leverage modern CMake modules like FetchContent
  - Use CMakePresets.json for standard build configurations
  - Configure proper package exports with CMake package config files
- Create proper Python package installation
  - Use pyproject.toml for modern Python packaging
  - Ensure proper integration between C++ and Python components
  - Implement consistent Cython module building
- Simplify dependency management
  - Centralize version requirements in a single location
  - Provide fallback mechanisms for missing system dependencies
  - Use standard CMake modules to find dependencies
- Make installation more user-friendly
  - Support both development and deployment installations
  - Provide clear error messages for missing dependencies
  - Ensure cross-platform compatibility

## Documentation
- Create comprehensive API documentation
  - Add Google-style docstrings to all public functions
  - Generate API documentation with Sphinx
  - Document parameter types and return values
  - Add examples in docstrings
  - Create module-level documentation
- Add user guides and examples
  - Write step-by-step tutorials for common use cases
  - Create Jupyter notebook examples
  - Document command-line options and configuration
  - Add troubleshooting guides for common issues
  - Create visual diagrams for complex workflows
- Document installation process clearly
  - Update README with modern installation instructions
  - Document environment requirements
  - Add platform-specific installation notes
  - Document dependency requirements
  - Create quick-start guide for new users
- Add contributor guidelines
  - Document code style requirements
  - Create pull request templates
  - Document test requirements for new code
  - Add development environment setup instructions
  - Create architectural overview documentation

## Performance Improvements
- Identify and optimize bottlenecks
  - Profile code execution with standard tools
  - Measure memory and CPU usage
  - Create benchmarking suite for key operations
  - Focus optimization on critical paths
  - Use profiling-guided optimization
- Improve memory usage for large genomic datasets
  - Implement streaming operations where possible
  - Use memory-mapped files for large data
  - Optimize data structures for genomic data
  - Implement chunking for large operations
  - Consider compression for in-memory representations
- Optimize file I/O operations
  - Use buffered I/O operations
  - Implement efficient file parsing
  - Use indexed file formats where appropriate
  - Add caching for frequently accessed data
  - Optimize compression/decompression operations
- Consider parallelization where appropriate
  - Identify parallelizable operations
  - Implement multi-threading for computation
  - Use process pools for independent tasks
  - Balance parallelism with overhead
  - Add asynchronous processing for I/O bound tasks
