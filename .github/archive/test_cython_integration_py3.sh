#!/bin/bash
# test_cython_integration_py3.sh - Test Cython modules with Python 3
#
# This script builds and tests the Cython modules for the hap.py Python 3 migration.
# It verifies string handling, memory management, and basic functionality.

set -e

# Set directory variables
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
BUILD_DIR="${1:-/tmp/happy-py3-cython-test}"
LOG_DIR="${SCRIPT_DIR}/build_logs"
VERBOSE="${VERBOSE:-0}"

# Print colored status messages
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${BLUE}=== hap.py Python 3 Cython Integration Test ===${NC}"
echo "Build directory: ${BUILD_DIR}"

# Create log directory if it doesn't exist
mkdir -p "${LOG_DIR}"

# Function to check if Python 3 and required packages are available
check_python() {
    if ! command -v python3 &> /dev/null; then
        echo -e "${RED}Error: Python 3 is required but not found in PATH${NC}"
        exit 1
    fi

    PYTHON_VERSION=$(python3 -c "import sys; print(f'{sys.version_info.major}.{sys.version_info.minor}')")
    echo -e "${GREEN}Using Python ${PYTHON_VERSION}${NC}"

    # Check for required packages
    python3 -c "import numpy, cython" &> /dev/null || {
        echo -e "${RED}Error: Required Python packages not found. Please install:${NC}"
        echo "  pip install numpy cython"
        exit 1
    }

    # Get Python site-packages directory for later use
    PYTHON_SITE_PACKAGES=$(python3 -c "import site; print(site.getsitepackages()[0])")
    echo "Python site-packages: ${PYTHON_SITE_PACKAGES}"
}

# Function to build Cython modules
build_cython_modules() {
    echo -e "${BLUE}Building Cython modules...${NC}"

    # Create build directory if it doesn't exist
    mkdir -p "${BUILD_DIR}"

    # Copy CythonSupport.py3.cmake to the build directory
    cp "${SCRIPT_DIR}/src/cmake/CythonSupport.py3.cmake" "${BUILD_DIR}/CythonSupport.cmake"

    # Create a minimal CMake file for testing
    TEST_CMAKE_FILE="${BUILD_DIR}/CythonTest.cmake"

    cat > "${TEST_CMAKE_FILE}" << EOF
cmake_minimum_required(VERSION 3.10)
project(HAPPY_CYTHON_TEST LANGUAGES CXX C)

# Include our custom Cython support
include(CythonSupport.cmake)

# Set up Python and paths
find_package(Python3 COMPONENTS Interpreter Development NumPy REQUIRED)

# Add Cython modules for testing
add_cython_module(
    _cython_test
    ${SCRIPT_DIR}/src/python/cython_test.pyx
)
EOF

    # Create a test Cython module
    mkdir -p "${SCRIPT_DIR}/src/python"
    TEST_CYTHON_FILE="${SCRIPT_DIR}/src/python/cython_test.pyx"

    cat > "${TEST_CYTHON_FILE}" << EOF
# cython: language_level=3
# distutils: language=c++
"""
Test Cython module for Python 3 compatibility testing
"""
import numpy as np
cimport numpy as np

from cpython.ref cimport PyObject
from libc.stdlib cimport malloc, free
from libc.string cimport memcpy, strlen
from libcpp.string cimport string
from libcpp.vector cimport vector

# Test string handling functions
def test_string_handling():
    """Test string handling with Python 3"""
    # Python string to C string
    cdef bytes py_bytes
    cdef char* c_str

    s = "Hello, world!"
    py_bytes = s.encode('utf-8')
    c_str = py_bytes

    # C string to Python string
    result = c_str.decode('utf-8')

    # Test with non-ASCII characters
    unicode_str = "Unicode: äöü 你好 🧬"
    py_bytes = unicode_str.encode('utf-8')
    c_str = py_bytes
    result2 = c_str.decode('utf-8')

    return {
        "orig": s,
        "roundtrip": result,
        "unicode_orig": unicode_str,
        "unicode_roundtrip": result2
    }

# Test memory management
def test_memory_management():
    """Test memory management with Python 3"""
    cdef char* buffer
    cdef Py_ssize_t size = 100

    # Allocate memory
    buffer = <char*>malloc(size * sizeof(char))
    if buffer == NULL:
        raise MemoryError("Failed to allocate memory")

    try:
        # Use the buffer
        for i in range(size-1):
            buffer[i] = 65 + (i % 26)  # ASCII A-Z
        buffer[size-1] = 0  # Null terminate

        result = buffer.decode('utf-8')
    finally:
        # Clean up
        free(buffer)

    return result

# Test NumPy integration
def test_numpy_integration():
    """Test NumPy integration with Python 3"""
    # Create a NumPy array
    cdef np.ndarray[np.float64_t, ndim=1] arr = np.linspace(0, 10, 10)

    # Modify the array
    cdef Py_ssize_t i
    for i in range(arr.shape[0]):
        arr[i] = arr[i] * 2

    return arr.tolist()

# Test basic functionality
def test_basic_functionality():
    """Test basic functionality with Python 3"""
    results = {
        "string_test": test_string_handling(),
        "memory_test": test_memory_management(),
        "numpy_test": test_numpy_integration()
    }
    return results
EOF

    # Configure and build
    cd "${BUILD_DIR}"
    echo -e "${BLUE}Configuring CMake...${NC}"
    cmake -S . -B build -G "Unix Makefiles" -DCMAKE_BUILD_TYPE=Debug

    echo -e "${BLUE}Building Cython modules...${NC}"
    cmake --build build -j$(sysctl -n hw.ncpu 2>/dev/null || echo 2)

    # Copy the built module to a location where Python can find it
    SITE_PACKAGES_DIR="${BUILD_DIR}/site-packages"
    mkdir -p "${SITE_PACKAGES_DIR}"
    cp "${BUILD_DIR}/build/lib/_cython_test."* "${SITE_PACKAGES_DIR}/"

    echo -e "${GREEN}Cython modules built successfully${NC}"
}

# Function to test Cython modules
test_cython_modules() {
    echo -e "${BLUE}Testing Cython modules...${NC}"

    # Create test Python script
    TEST_PY_FILE="${BUILD_DIR}/test_cython.py"

    cat > "${TEST_PY_FILE}" << EOF
#!/usr/bin/env python3
"""
Test the built Cython module
"""
import os
import sys
import json

# Add site-packages to Python path
site_packages_dir = os.path.join("${BUILD_DIR}", "site-packages")
sys.path.insert(0, site_packages_dir)

try:
    # Import the Cython module
    import _cython_test
    print("Successfully imported _cython_test module")

    # Test string handling
    print("\nTesting string handling...")
    string_results = _cython_test.test_string_handling()
    print(json.dumps(string_results, indent=2))
    assert string_results["orig"] == string_results["roundtrip"]
    assert string_results["unicode_orig"] == string_results["unicode_roundtrip"]
    print("String handling test passed")

    # Test memory management
    print("\nTesting memory management...")
    memory_result = _cython_test.test_memory_management()
    print(f"Memory test result: {memory_result}")
    assert len(memory_result) == 99  # 100-1 (null terminator)
    print("Memory management test passed")

    # Test NumPy integration
    print("\nTesting NumPy integration...")
    numpy_result = _cython_test.test_numpy_integration()
    print(f"NumPy test result: {numpy_result}")
    assert len(numpy_result) == 10
    print("NumPy integration test passed")

    # Test basic functionality
    print("\nTesting basic functionality...")
    basic_result = _cython_test.test_basic_functionality()
    print("Basic functionality tested successfully")

    print("\nAll tests passed!")
    sys.exit(0)
except ImportError as e:
    print(f"Error importing Cython module: {e}")
    sys.exit(1)
except Exception as e:
    print(f"Test failed: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)
EOF

    # Make the test script executable
    chmod +x "${TEST_PY_FILE}"

    # Run the test
    echo -e "${BLUE}Running Cython module tests...${NC}"
    python3 "${TEST_PY_FILE}"

    if [ $? -eq 0 ]; then
        echo -e "${GREEN}Cython module tests passed!${NC}"
    else
        echo -e "${RED}Cython module tests failed!${NC}"
        exit 1
    fi
}

# Main script

# Check Python environment
check_python

# Build Cython modules
build_cython_modules

# Test Cython modules
test_cython_modules

echo -e "${GREEN}All tests completed successfully!${NC}"

# Create a minimal test report
TEST_REPORT="${LOG_DIR}/cython_integration_test_report.txt"
cat > "${TEST_REPORT}" << EOF
Cython Integration Test Report
$(date)

Python version: ${PYTHON_VERSION}
Build directory: ${BUILD_DIR}

All tests passed successfully.

Next steps:
1. Apply the Cython module updates to the real hap.py codebase
2. Update the CMake files to use CythonSupport.py3.cmake
3. Test with real genomic data files
EOF

echo -e "${BLUE}Test report written to ${TEST_REPORT}${NC}"
