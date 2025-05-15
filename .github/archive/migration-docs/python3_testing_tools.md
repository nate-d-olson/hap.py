# Python 3 Testing Tools

This document explains the various tools created to test and validate the Python 3 migration of hap.py.

## Overview

The migration from Python 2 to Python 3 for hap.py involves updating not just Python code but also the C++ integration through Cython. To ensure a smooth migration, we've created several specialized testing tools.

## Testing Tools

### 1. Cython Module Testing

**`test_cython_module_py3.py`**

This script tests if the Cython modules can be loaded in Python 3 and verifies basic functionality, particularly focusing on string handling and memory management.

```bash
python test_cython_module_py3.py --build-dir /path/to/build/dir [--module module_name]
```

Options:
- `--build-dir`: Path to the build directory
- `--module`: Optional specific module to test
- `--verbose`: Enable verbose output

### 2. Cython Module Updating

**`update_cython_modules_py3.py`**

This script updates Cython module files to ensure they're compatible with Python 3. It adds language level directives and fixes string handling.

```bash
python update_cython_modules_py3.py --src-dir /path/to/src/dir [--dry-run]
```

Options:
- `--src-dir`: Path to the source directory
- `--dry-run`: Don't make any changes, just report what would be done

### 3. Comprehensive Cython Integration Testing

**`test_cython_integration_py3.sh`**

This script builds and tests Cython modules to verify their compatibility with Python 3.

```bash
bash test_cython_integration_py3.sh [/path/to/build/dir]
```

This script:
1. Creates a minimal Cython test module
2. Builds it using CMake and the updated CythonSupport.py3.cmake
3. Tests string handling, memory management, and NumPy integration

### 4. Main hap.py Testing

**`test_happy_py3.sh`**

This script tests the main hap.py functionality with Python 3.

```bash
bash test_happy_py3.sh [/path/to/build/dir]
```

This script:
1. Builds a minimal test installation
2. Tests basic functionality like help and version display
3. Tests running with the mock implementation

### 5. C++ Integration Validation

**`validate_cpp_integration_enhanced.py`**

This script validates the integration between Python 3 modules and C++ components.

```bash
python validate_cpp_integration_enhanced.py --build-dir /path/to/build --test-type [variant|reference|all] [--use-mock]
```

Options:
- `--build-dir`: Path to the build directory
- `--test-type`: Type of test to run
- `--use-mock`: Use mock implementation instead of real C++ components

## Using These Tools During Migration

During the Python 3 migration process, we recommend using these tools in the following order:

1. First, run `update_cython_modules_py3.py` to update Cython modules
2. Then run `test_cython_module_py3.py` to verify basic loading works
3. Use `test_cython_integration_py3.sh` to test comprehensive Cython functionality
4. Finally, run `test_happy_py3.sh` to test the full application

For continuous integration and build verification, all these tests can be chained together using `test_py3_build.sh`.

## Mock Implementation Strategy

For testing without requiring the full C++ build, we've implemented a mock implementation strategy:

1. Each Cython module has a corresponding Python mock implementation
2. The `__init__.py` files attempt to import the Cython modules but fall back to the mock implementation if they're not available
3. Set `HAPLO_USE_MOCK=1` to force using mock implementations

This allows developing and testing the Python 3 migration even on systems where C++ components cannot be built.

## CMake Integration

For proper Cython and Python 3 integration, we've created an improved CMake module:

**`CythonSupport.py3.cmake`**

This file provides:
1. Better detection of Python 3 and NumPy
2. Proper language level setting for Cython modules
3. Platform-specific extension handling

## Next Steps

After these testing tools confirm Python 3 compatibility, the next steps are:

1. Complete the migration of all remaining Python modules
2. Finalize the build system updates
3. Create a unified installer for both Python 2 and 3
4. Update all documentation
