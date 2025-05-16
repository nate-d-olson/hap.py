#!/bin/bash
# test_py3_core.sh - Test the core functionality of hap.py with Python 3
#
# This script tests the core functionality of hap.py with vcfeval as the
# comparison engine to ensure it works properly after streamlining.

set -e

# Set directory variables
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BUILD_DIR="${1:-/tmp/happy-py3-core-test}"
LOG_DIR="${SCRIPT_DIR}/build_logs"

# Print colored status messages
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${BLUE}=== Testing hap.py core functionality with Python 3 ===${NC}"

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
    
    # Run the installer script with Python 3
    python3 "${SCRIPT_DIR}/install.py" "${BUILD_DIR}" --no-tests
    
    if [ $? -ne 0 ]; then
        echo -e "${RED}Build failed!${NC}"
        exit 1
    else
        echo -e "${GREEN}Build succeeded!${NC}"
    fi
}

# Test the core functionality with example data
test_core_functionality() {
    echo -e "${BLUE}Testing core functionality with vcfeval engine...${NC}"
    
    # Set up paths
    HAPPY_PATH="${BUILD_DIR}/bin/hap.py"
    EXAMPLE_DIR="${SCRIPT_DIR}/example"
    OUTPUT_DIR="${BUILD_DIR}/test_output"
    
    mkdir -p "${OUTPUT_DIR}"
    
    # Run hap.py with vcfeval engine on example data
    echo -e "${YELLOW}Running hap.py with vcfeval engine...${NC}"
    ${HAPPY_PATH} \
        ${EXAMPLE_DIR}/PG_performance.vcf.gz \
        ${EXAMPLE_DIR}/PG_hc.vcf.gz \
        -r ${EXAMPLE_DIR}/chr21.fa \
        -o ${OUTPUT_DIR}/test.vcfeval \
        --engine vcfeval \
        -f ${EXAMPLE_DIR}/performance.confident.bed.gz \
        --ci-alpha 0.05 \
        --no-decompose \
        --no-leftshift \
        > "${LOG_DIR}/hap.py.vcfeval.log" 2>&1
    
    if [ $? -ne 0 ]; then
        echo -e "${RED}Core functionality test failed!${NC}"
        echo -e "${RED}Check the log file at ${LOG_DIR}/hap.py.vcfeval.log${NC}"
        exit 1
    else
        echo -e "${GREEN}Core functionality test succeeded!${NC}"
    fi
    
    # Check if output files were created
    if [ ! -f "${OUTPUT_DIR}/test.vcfeval.summary.csv" ]; then
        echo -e "${RED}Output file not found: ${OUTPUT_DIR}/test.vcfeval.summary.csv${NC}"
        exit 1
    else
        echo -e "${GREEN}Output files created successfully:${NC}"
        ls -la "${OUTPUT_DIR}/test.vcfeval"*
    fi
}

# Main execution
check_python
build_test_installation
test_core_functionality

echo -e "${GREEN}All tests completed successfully!${NC}"
