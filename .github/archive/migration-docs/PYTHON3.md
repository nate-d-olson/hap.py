# Python 3 Migration Guide for hap.py

This document provides instructions for using the Python 3 version of hap.py.

## Prerequisites

- Python 3.6 or newer
- CMake 3.10 or newer
- C++11 compatible compiler
- Necessary system libraries (zlib, bzip2, etc.)

## Installation

For a standard installation using Python 3:

```bash
# Clone the repository
git clone https://github.com/your-org/hap.py.git
cd hap.py

# Create a Python 3 virtual environment (recommended)
python3 -m venv venv
source venv/bin/activate

# Install dependencies
pip install -r happy.requirements.py3.txt

# Update Cython modules for Python 3 compatibility
python update_cython_modules_py3.py --src-dir src/python

# Build and install hap.py
python3 install_py3.py /path/to/build/dir
```

## Testing

After installation, you can verify the Python 3 compatibility:

```bash
# Test Cython module loading
python test_cython_module_py3.py --build-dir /path/to/build/dir

# Test hap.py functionality
bash test_happy_py3.sh /path/to/build/dir

# Run the comprehensive integration test
bash test_cython_integration_py3.sh /path/to/build/dir
```

## Key Changes from Python 2 Version

- Updated string handling for Unicode compatibility
- Modern Python packaging structure with type hints
- Improved error handling and reporting
- Better integration with C++ components through Cython
- Updated dependency requirements
- More robust fallback mechanisms when C++ components are unavailable

## Using with Mock C++ Implementation

During development or in environments where C++ compilation is not possible,
you can use the mock implementation:

```bash
# Set environment variable before importing
export HAPLO_USE_MOCK=1

# Run hap.py scripts
python3 /path/to/build/dir/bin/hap.py [arguments]
```

## Testing the Integration between Python 3 and C++

To verify that the Python 3 code correctly integrates with the C++ components,
you can use the validation script:

```bash
# Test with actual C++ implementation
python3 validate_cpp_integration.py --build-dir /path/to/build/dir

# Test with mock implementation
python3 validate_cpp_integration.py --build-dir /path/to/build/dir --use-mock
```

## Migrated Python Modules

The following modules have been migrated to Python 3:

- Core components: `hap.py`, `qfy.py`, `pre.py`
- Utility modules in `Tools` package including `bcftools.py`, `fastasize.py`, `bedintervaltree.py`
- Key analysis modules in `Haplo` package
- Statistical components in `Somatic` package

Each Python module has been updated to use:

1. Consistent string handling (bytes vs. unicode) for Python 3
2. Modern iterator patterns (replacing `xrange` with `range`, etc.)
3. Type annotations for better code documentation and IDE support
4. Context managers for file handling
5. Updated exception handling

## Build System Improvements

The build system has been modernized to support Python 3:

1. CMake integration with Cython via `CythonSupport.cmake`
2. Better dependency management with version checks
3. Improved error reporting during build failures
4. Fallback mechanisms when components cannot be built

## Testing the Installation

To verify your installation:

```bash
# Test Cython integration
python3 test_cython_integration.py --test-both

# Run the test suite
cd /path/to/build/dir
src/sh/run_tests.sh
```

## Known Issues

- Some advanced features may not be available in the mock implementation
- Performance might differ between Python 2 and Python 3 versions
- Certain dependencies may require specific versions for compatibility

## Reporting Problems

If you encounter issues, please run the diagnostic script:

```bash
./diagnose_build.sh /path/to/build/dir
```

This will generate detailed logs in the `build_logs` directory that can help diagnose build problems.
