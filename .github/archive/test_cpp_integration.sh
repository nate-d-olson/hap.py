#!/bin/bash
# test_cpp_integration.sh - Test C++ integration with Python 3
#
# This script helps test the integration between Python 3 modules and C++ components
# in the hap.py project. It builds the C++ components if needed, sets up a test
# environment, and runs the validation scripts.

set -e

# Set directory variables
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
BUILD_DIR="${1:-/tmp/happy-py3-test}"
LOG_DIR="${SCRIPT_DIR}/build_logs"

# Print colored status messages
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${BLUE}=== hap.py Python 3 C++ Integration Test ===${NC}"

# Create log directory if it doesn't exist
mkdir -p "${LOG_DIR}"

# Function to check if Python 3 is available
check_python() {
    if ! command -v python3 &> /dev/null; then
        echo -e "${RED}Error: Python 3 is required but not found in PATH${NC}"
        exit 1
    fi

    PYTHON_VERSION=$(python3 -c "import sys; print(f'{sys.version_info.major}.{sys.version_info.minor}')")
    echo -e "${GREEN}Using Python ${PYTHON_VERSION}${NC}"

    # Check for required packages
    python3 -c "import numpy, pandas, cython" &> /dev/null || {
        echo -e "${RED}Error: Required Python packages not found. Please install:${NC}"
        echo "  pip install numpy pandas cython"
        exit 1
    }
}

# Function to build C++ components
build_cpp() {
    echo -e "${BLUE}Building C++ components...${NC}"

    # Create build directory if it doesn't exist
    mkdir -p "${BUILD_DIR}"

    # Use Python 3 install script
    python3 "${SCRIPT_DIR}/install_py3.py" "${BUILD_DIR}" 2>&1 | tee "${LOG_DIR}/cpp_build.log"

    if [ ${PIPESTATUS[0]} -ne 0 ]; then
        echo -e "${RED}Failed to build C++ components. See ${LOG_DIR}/cpp_build.log for details${NC}"

        # Try to extract specific errors
        grep -A 10 "error:" "${LOG_DIR}/cpp_build.log" > "${LOG_DIR}/cpp_errors.log" || true
        echo -e "${YELLOW}Extracted error messages to ${LOG_DIR}/cpp_errors.log${NC}"

        exit 1
    fi

    echo -e "${GREEN}Successfully built C++ components${NC}"
}

# Function to run tests with mock implementation
test_mock() {
    echo -e "${BLUE}Testing with mock C++ implementation...${NC}"

    # Set environment variable to use mock
    export HAPLO_USE_MOCK=1

    # Run enhanced validation script
    python3 "${SCRIPT_DIR}/validate_cpp_integration_enhanced.py" --build-dir "${BUILD_DIR}" --use-mock --verbose

    if [ $? -ne 0 ]; then
        echo -e "${RED}Mock implementation tests failed${NC}"
        return 1
    else
        echo -e "${GREEN}Mock implementation tests passed${NC}"
        return 0
    fi
}

# Function to run tests with actual C++ implementation
test_real() {
    echo -e "${BLUE}Testing with real C++ implementation...${NC}"

    # Clear environment variable
    unset HAPLO_USE_MOCK

    # Run enhanced validation script
    python3 "${SCRIPT_DIR}/validate_cpp_integration_enhanced.py" --build-dir "${BUILD_DIR}" --verbose

    if [ $? -ne 0 ]; then
        echo -e "${RED}Real implementation tests failed${NC}"
        return 1
    else
        echo -e "${GREEN}Real implementation tests passed${NC}"
        return 0
    fi
}

# Function to test Python 3 modules
test_py3_modules() {
    echo -e "${BLUE}Testing Python 3 modules...${NC}"

    python3 "${SCRIPT_DIR}/test_python3_compatibility_enhanced.py" --build-dir "${BUILD_DIR}" --verbose

    if [ $? -ne 0 ]; then
        echo -e "${RED}Python 3 module tests failed${NC}"
        return 1
    else
        echo -e "${GREEN}Python 3 module tests passed${NC}"
        return 0
    fi
}

# Main execution
echo -e "${YELLOW}Test build directory: ${BUILD_DIR}${NC}"

# Check Python
check_python

# Build C++ components if not already built or if requested
if [ ! -d "${BUILD_DIR}/lib" ] || [ "$2" == "--rebuild" ]; then
    build_cpp
else
    echo -e "${YELLOW}Using existing build in ${BUILD_DIR}${NC}"
    echo -e "${YELLOW}(use --rebuild to force rebuild)${NC}"
fi

# Run tests
SUCCESS=0

echo -e "${BLUE}Running integration tests...${NC}"

# Test Python 3 modules first
test_py3_modules || SUCCESS=1

# Test with mock implementation
test_mock || SUCCESS=1

# Test with real implementation
test_real || SUCCESS=1

# Print summary
if [ $SUCCESS -eq 0 ]; then
    echo -e "${GREEN}All tests passed!${NC}"
    exit 0
else
    echo -e "${RED}Some tests failed. See logs for details.${NC}"
    exit 1
fi
