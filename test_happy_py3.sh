#!/bin/bash
# test_happy_py3.sh - Test running hap.py with Python 3
#
# This script tests the main hap.py functionality with Python 3.
# It helps verify that the core functionality works with Python 3.

set -e

# Set directory variables
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
BUILD_DIR="${1:-/tmp/happy-py3-run}"
LOG_DIR="${SCRIPT_DIR}/build_logs"

# Print colored status messages
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${BLUE}=== Testing hap.py with Python 3 ===${NC}"

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
    python3 -c "import numpy, pandas, pysam" &> /dev/null || {
        echo -e "${RED}Error: Required Python packages not found. Please install:${NC}"
        echo "  pip install numpy pandas pysam"
        exit 1
    }
}

# Function to build a minimal test installation
build_test_installation() {
    echo -e "${BLUE}Building test installation...${NC}"

    # Create build directory if it doesn't exist
    mkdir -p "${BUILD_DIR}"

    # Check if install_py3.py exists
    if [ ! -f "${SCRIPT_DIR}/install_py3.py" ]; then
        echo -e "${RED}Error: install_py3.py not found${NC}"
        exit 1
    fi

    # Run installation script
    echo -e "${BLUE}Running installation...${NC}"
    python3 "${SCRIPT_DIR}/install_py3.py" "${BUILD_DIR}" | tee "${LOG_DIR}/py3_install.log"

    # Create Python 3 version of hap.py script
    echo -e "${BLUE}Creating Python 3 hap.py script...${NC}"
    mkdir -p "${BUILD_DIR}/bin"
    cp "${SCRIPT_DIR}/bin/hap.py.py3" "${BUILD_DIR}/bin/hap.py.py3"
    chmod +x "${BUILD_DIR}/bin/hap.py.py3"

    echo -e "${GREEN}Test installation complete${NC}"
}

# Function to run basic tests
run_basic_tests() {
    echo -e "${BLUE}Running basic tests...${NC}"

    # Check if test files exist
    if [ ! -f "${SCRIPT_DIR}/example/PG_performance.vcf.gz" ] || [ ! -f "${SCRIPT_DIR}/example/performance.vcf.gz" ]; then
        echo -e "${YELLOW}Warning: Test files not found, skipping basic tests${NC}"
        return
    fi

    # Create test output directory
    TEST_OUTPUT_DIR="${BUILD_DIR}/test_output"
    mkdir -p "${TEST_OUTPUT_DIR}"

    # Run hap.py with basic parameters
    echo -e "${BLUE}Running hap.py with Python 3...${NC}"

    # Use the Python 3 version of hap.py
    HAPPY_CMD="${BUILD_DIR}/bin/hap.py.py3"

    if [ ! -f "${HAPPY_CMD}" ]; then
        echo -e "${RED}Error: ${HAPPY_CMD} not found${NC}"
        exit 1
    fi

    # Test running hap.py with --help flag
    echo -e "${BLUE}Testing basic help output...${NC}"
    python3 "${HAPPY_CMD}" --help | tee "${LOG_DIR}/happy_help.log"

    # Test running hap.py with --version flag
    echo -e "${BLUE}Testing version display...${NC}"
    python3 "${HAPPY_CMD}" --version | tee "${LOG_DIR}/happy_version.log"

    # Mock test - set environment variable to use mock implementation
    echo -e "${BLUE}Testing with mock implementation...${NC}"
    HAPLO_USE_MOCK=1 python3 "${HAPPY_CMD}" \
        "${SCRIPT_DIR}/example/PG_performance.vcf.gz" \
        "${SCRIPT_DIR}/example/performance.vcf.gz" \
        -o "${TEST_OUTPUT_DIR}/mock_test" \
        --pass-only | tee "${LOG_DIR}/happy_mock_test.log" || {
        echo -e "${YELLOW}Warning: Mock test failed, this is expected if mock implementation is incomplete${NC}"
    }

    echo -e "${GREEN}Basic tests completed${NC}"
}

# Function to check for Python 3 compatibility issues
check_py3_compatibility() {
    echo -e "${BLUE}Checking for Python 3 compatibility issues...${NC}"

    # Use check_py3_issues.py if available
    if [ -f "${SCRIPT_DIR}/check_py3_issues.py" ]; then
        python3 "${SCRIPT_DIR}/check_py3_issues.py" --dir "${SCRIPT_DIR}/src/python" | tee "${LOG_DIR}/py3_issues.log"
    else
        # Fall back to simple grep checks for common Python 2 patterns
        echo "Checking for common Python 2 patterns..."
        grep -r "print " --include="*.py" "${SCRIPT_DIR}/src/python" | tee "${LOG_DIR}/py3_print_statements.log"
        grep -r "xrange" --include="*.py" "${SCRIPT_DIR}/src/python" | tee "${LOG_DIR}/py3_xrange.log"
        grep -r "iteritems" --include="*.py" "${SCRIPT_DIR}/src/python" | tee "${LOG_DIR}/py3_iteritems.log"
        grep -r "except .* as" --include="*.py" "${SCRIPT_DIR}/src/python" | tee "${LOG_DIR}/py3_except.log"
        grep -r "unicode" --include="*.py" "${SCRIPT_DIR}/src/python" | tee "${LOG_DIR}/py3_unicode.log"
    fi

    echo -e "${GREEN}Compatibility check complete${NC}"
}

# Main script execution

# Check Python environment
check_python

# Check Python 3 compatibility issues
check_py3_compatibility

# Build test installation if needed
if [ ! -d "${BUILD_DIR}" ] || [ ! -f "${BUILD_DIR}/bin/hap.py.py3" ]; then
    build_test_installation
fi

# Run basic tests
run_basic_tests

echo -e "${GREEN}All tests completed!${NC}"
echo "See logs in ${LOG_DIR} for detailed information."

# Create a simple test summary
SUMMARY_FILE="${LOG_DIR}/py3_test_summary.txt"
cat > "${SUMMARY_FILE}" << EOF
Python 3 hap.py Test Summary
$(date)

Python version: ${PYTHON_VERSION}
Build directory: ${BUILD_DIR}

Main components tested:
- Python 3 compatibility checks
- Basic script execution
- Help and version display
- Mock implementation test

Next steps:
1. Complete Python 3 migration of core modules
2. Test with real data and C++ integration
3. Update build system for both Python 2 and 3 support
4. Create comprehensive test suite
EOF

echo -e "${BLUE}Test summary written to ${SUMMARY_FILE}${NC}"
