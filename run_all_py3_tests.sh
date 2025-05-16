#!/bin/bash
# run_all_py3_tests.sh - Run all Python 3 migration tests for hap.py
#
# This script runs the full suite of tests for the Python 3 migration of hap.py.
# It verifies all aspects of the migration, from Cython modules to full application.

set -e

# Set directory variables
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
BUILD_DIR="${1:-/tmp/happy-py3-full-test}"
LOG_DIR="${SCRIPT_DIR}/build_logs"
TEST_REPORT="${LOG_DIR}/py3_migration_test_report.md"

# Print colored status messages
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${BLUE}=== hap.py Python 3 Migration Full Test Suite ===${NC}"

# Create log directory if it doesn't exist
mkdir -p "${LOG_DIR}"

# Check if Python 3 is available
if ! command -v python3 &> /dev/null; then
    echo -e "${RED}Error: Python 3 is required but not found in PATH${NC}"
    exit 1
fi

PYTHON_VERSION=$(python3 -c "import sys; print(f'{sys.version_info.major}.{sys.version_info.minor}')")
echo -e "${GREEN}Using Python ${PYTHON_VERSION}${NC}"

# Create test report header
cat > "${TEST_REPORT}" << EOF
# Python 3 Migration Test Report

Date: $(date)
Python Version: ${PYTHON_VERSION}
Build Directory: ${BUILD_DIR}

## Test Results

EOF

# Function to run a test and add result to report
run_test() {
    local test_name=$1
    local test_cmd=$2
    local test_desc=$3
    local log_file="${LOG_DIR}/${test_name}.log"

    echo -e "${BLUE}Running test: ${test_name}${NC}"
    echo "${test_desc}"

    # Add to report
    cat >> "${TEST_REPORT}" << EOF
### ${test_name}

${test_desc}

\`\`\`bash
${test_cmd}
\`\`\`

EOF

    # Run the test
    set +e
    eval "${test_cmd}" > "${log_file}" 2>&1
    local result=$?
    set -e

    # Add result to report
    if [ ${result} -eq 0 ]; then
        echo -e "${GREEN}✓ Test passed${NC}"
        cat >> "${TEST_REPORT}" << EOF
**Result: PASS** ✓

EOF
    else
        echo -e "${RED}✗ Test failed${NC}"
        cat >> "${TEST_REPORT}" << EOF
**Result: FAIL** ✗

Error summary:
\`\`\`
$(tail -n 20 "${log_file}")
\`\`\`

EOF
    fi

    return ${result}
}

# Check if all required test scripts are available
missing_scripts=0
for script in update_cython_modules_py3.py test_cython_module_py3.py test_cython_integration_py3.sh test_happy_py3.sh; do
    if [ ! -f "${SCRIPT_DIR}/${script}" ]; then
        echo -e "${RED}Error: Required test script ${script} not found${NC}"
        missing_scripts=$((missing_scripts + 1))
    fi
done

if [ ${missing_scripts} -gt 0 ]; then
    echo -e "${RED}Error: ${missing_scripts} required test scripts not found. Cannot proceed.${NC}"
    exit 1
fi

# Run tests in sequence
tests_failed=0

# Test 1: Python 3 compatibility check
run_test "py3_compatibility_check" \
    "python3 ${SCRIPT_DIR}/check_py3_issues.py --dir ${SCRIPT_DIR}/src/python" \
    "Check for common Python 2 to 3 compatibility issues" || tests_failed=$((tests_failed + 1))

# Test 2: Update Cython modules
run_test "update_cython_modules" \
    "python3 ${SCRIPT_DIR}/update_cython_modules_py3.py --src-dir ${SCRIPT_DIR}/src/python" \
    "Update Cython modules for Python 3 compatibility" || tests_failed=$((tests_failed + 1))

# Test 3: Test Cython module loading
run_test "test_cython_modules" \
    "python3 ${SCRIPT_DIR}/test_cython_module_py3.py --build-dir ${BUILD_DIR}" \
    "Test loading Cython modules in Python 3" || tests_failed=$((tests_failed + 1))

# Test 4: Comprehensive Cython integration test
run_test "cython_integration_test" \
    "bash ${SCRIPT_DIR}/test_cython_integration_py3.sh ${BUILD_DIR}/cython_test" \
    "Test comprehensive Cython integration with Python 3" || tests_failed=$((tests_failed + 1))

# Test 5: Build test
run_test "py3_build_test" \
    "bash ${SCRIPT_DIR}/test_py3_build.sh ${BUILD_DIR}/build_test" \
    "Test full build with Python 3" || tests_failed=$((tests_failed + 1))

# Test 6: hap.py test
run_test "happy_py3_test" \
    "bash ${SCRIPT_DIR}/test_happy_py3.sh ${BUILD_DIR}/happy_test" \
    "Test running hap.py with Python 3" || tests_failed=$((tests_failed + 1))

# Test 7: Mock implementation test
run_test "mock_implementation_test" \
    "HAPLO_USE_MOCK=1 python3 ${SCRIPT_DIR}/test_cython_module_py3.py --build-dir ${BUILD_DIR}" \
    "Test using mock implementations instead of C++ components" || tests_failed=$((tests_failed + 1))

# Test 8: C++ integration validation
if [ -f "${SCRIPT_DIR}/validate_cpp_integration_enhanced.py" ]; then
    run_test "cpp_integration_validation" \
        "python3 ${SCRIPT_DIR}/validate_cpp_integration_enhanced.py --build-dir ${BUILD_DIR}/build_test" \
        "Validate C++ integration with Python 3" || tests_failed=$((tests_failed + 1))
fi

# Add summary to report
cat >> "${TEST_REPORT}" << EOF
## Summary

Total tests: 8
Passed: $((8 - tests_failed))
Failed: ${tests_failed}

EOF

if [ ${tests_failed} -eq 0 ]; then
    cat >> "${TEST_REPORT}" << EOF
**All tests passed successfully!** ✓

## Next Steps

1. Complete the migration of any remaining Python modules
2. Finalize the build system updates
3. Create a unified installer for both Python 2 and 3
4. Update all documentation
EOF
    echo -e "${GREEN}All tests passed successfully!${NC}"
else
    cat >> "${TEST_REPORT}" << EOF
**Some tests failed.** ✗

## Recommended Actions

1. Review the logs in ${LOG_DIR} for detailed error information
2. Fix the identified issues
3. Re-run the tests
EOF
    echo -e "${RED}${tests_failed} tests failed. See ${TEST_REPORT} for details.${NC}"
fi

echo -e "${BLUE}Test report written to ${TEST_REPORT}${NC}"

exit ${tests_failed}
