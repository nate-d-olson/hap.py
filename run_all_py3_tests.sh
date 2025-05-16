#\!/bin/bash
#
# Run all Python 3 tests for hap.py
#
# This script will run all the Python 3 tests for the hap.py codebase.
# It assumes that the Python 3 environment is already set up and that
# the necessary dependencies are installed.
#
# Usage: ./run_all_py3_tests.sh [--build-dir /path/to/build]

set -e

# Parse command line arguments
BUILD_DIR=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        --build-dir)
            BUILD_DIR="$2"
            shift 2
            ;;
        *)
            echo "Unknown argument: $1"
            echo "Usage: $0 [--build-dir /path/to/build]"
            exit 1
            ;;
    esac
done

# Set default build directory if not specified
if [ -z "$BUILD_DIR" ]; then
    BUILD_DIR="build"
    echo "Using default build directory: $BUILD_DIR"
    mkdir -p "$BUILD_DIR"
fi

# Check if Python 3 is available
if \! command -v python3 &> /dev/null; then
    echo "Error: Python 3 is not available. Please install Python 3 before running this script."
    exit 1
fi

echo "Running all Python 3 tests for hap.py..."
echo "----------------------------------------"

# Function to run a test and report its status
run_test() {
    local test_name="$1"
    local test_command="$2"
    
    echo -n "Running $test_name... "
    if eval "$test_command"; then
        echo "PASSED"
        return 0
    else
        echo "FAILED"
        return 1
    fi
}

# Keep track of failures
FAILURES=0

# Test 1: Core Python 3 functionality
if [ -f ./test_py3_core.sh ]; then
    run_test "Core Python 3 functionality" "./test_py3_core.sh" || ((FAILURES++))
else
    echo "Skipping Core Python 3 functionality test (test_py3_core.sh not found)"
fi

# Test 2: Cython module loading
if [ -f ./test_cython_module_py3.py ]; then
    run_test "Cython module loading" "python3 test_cython_module_py3.py --build-dir $BUILD_DIR" || ((FAILURES++))
else
    echo "Skipping Cython module loading test (test_cython_module_py3.py not found)"
fi

# Test 3: Cython integration
if [ -f ./test_cython_integration_py3.sh ]; then
    run_test "Cython integration" "./test_cython_integration_py3.sh" || ((FAILURES++))
else
    echo "Skipping Cython integration test (test_cython_integration_py3.sh not found)"
fi

# Test 4: String handling
if [ -f ./test_string_handling.py ]; then
    run_test "String handling" "python3 test_string_handling.py" || ((FAILURES++))
else
    echo "Skipping String handling test (test_string_handling.py not found)"
fi

# Test 5: Python 3 compatibility
if [ -f ./test_python3_compatibility_enhanced.py ]; then
    run_test "Python 3 compatibility" "python3 test_python3_compatibility_enhanced.py" || ((FAILURES++))
else
    echo "Skipping Python 3 compatibility test (test_python3_compatibility_enhanced.py not found)"
fi

# Test 6: Check for Python 3 issues
if [ -f ./check_py3_issues.py ]; then
    run_test "Python 3 issues check" "python3 check_py3_issues.py --generate-report" || ((FAILURES++))
else
    echo "Skipping Python 3 issues check (check_py3_issues.py not found)"
fi

# Test 7: Build verification
if [ -f ./verify_py3_build.py ]; then
    run_test "Build verification" "python3 verify_py3_build.py --skip-test" || ((FAILURES++))
else
    echo "Skipping Build verification test (verify_py3_build.py not found)"
fi

# Print summary
echo ""
echo "Test Summary"
echo "-------------"
if [ $FAILURES -eq 0 ]; then
    echo "All tests PASSED\!"
    exit 0
else
    echo "$FAILURES test(s) FAILED\!"
    exit 1
fi
