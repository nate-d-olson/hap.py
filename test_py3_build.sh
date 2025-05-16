#!/usr/bin/env bash
# This script helps test the Python 3 migration of hap.py
# It creates a Python 3 virtual environment, installs dependencies,
# and attempts a test build using Python 3

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
NC='\033[0m' # No Color

echo -e "${GREEN}=== hap.py Python 3 Build Test ===${NC}"
echo -e "${YELLOW}This script will test building hap.py with Python 3${NC}"

# Create a temporary build directory
BUILD_DIR=$(mktemp -d)
echo -e "${GREEN}Using build directory: ${BUILD_DIR}${NC}"

# Create Python 3 virtual environment
echo -e "${GREEN}Creating Python 3 virtual environment...${NC}"
python3 -m venv "${BUILD_DIR}/venv_py3"
source "${BUILD_DIR}/venv_py3/bin/activate"

# Install dependencies
echo -e "${GREEN}Installing Python 3 dependencies...${NC}"
pip install --upgrade pip
pip install wheel numpy cython pandas scipy

# Check for Python 3 requirements file
if [ -f "happy.requirements.py3.txt" ]; then
    echo -e "${GREEN}Installing from happy.requirements.py3.txt...${NC}"
    pip install -r happy.requirements.py3.txt
else
    echo -e "${YELLOW}No Python 3 requirements file found. Using default requirements.${NC}"
    if [ -f "happy.requirements.txt" ]; then
        echo -e "${YELLOW}Installing from happy.requirements.txt (might not be Python 3 compatible)...${NC}"
        pip install -r happy.requirements.txt
    fi
fi

# Check for Cython module updater
if [ -f "update_cython_modules_py3.py" ]; then
    echo -e "${GREEN}Updating Cython modules for Python 3...${NC}"
    python update_cython_modules_py3.py --src-dir src/python 2>&1 | tee "${BUILD_DIR}/cython_update.log"
fi

# Check for Python 3 install script
if [ -f "install_py3.py" ]; then
    echo -e "${GREEN}Running Python 3 install script...${NC}"
    python install_py3.py "${BUILD_DIR}/build_py3" 2>&1 | tee "${BUILD_DIR}/build_py3.log"
else
    echo -e "${YELLOW}No Python 3 install script found. Attempting to use standard install script.${NC}"
    echo -e "${YELLOW}This may fail if the install script isn't Python 3 compatible.${NC}"
    python install.py "${BUILD_DIR}/build_py3" 2>&1 | tee "${BUILD_DIR}/build_py3.log"
fi

# Test Cython module loading
if [ -f "test_cython_module_py3.py" ]; then
    echo -e "${GREEN}Testing Cython module loading...${NC}"
    python test_cython_module_py3.py --build-dir "${BUILD_DIR}/build_py3" 2>&1 | tee "${BUILD_DIR}/cython_test.log"
fi

# Test hap.py with Python 3
if [ -f "test_happy_py3.sh" ]; then
    echo -e "${GREEN}Testing hap.py with Python 3...${NC}"
    bash test_happy_py3.sh "${BUILD_DIR}/build_py3" 2>&1 | tee "${BUILD_DIR}/happy_test.log"
fi

# Check if build was successful
if [ -d "${BUILD_DIR}/build_py3" ]; then
    echo -e "${GREEN}Build directory was created successfully.${NC}"

    # Test migrated modules
    echo -e "${GREEN}Testing migrated Python 3 modules...${NC}"
    python test_migrated_modules.py --src-dir "${BUILD_DIR}/build_py3/src"

    echo -e "${GREEN}Testing completed. Log files are available at:${NC}"
    echo -e "${YELLOW}${BUILD_DIR}/build_py3.log${NC}"
    echo
    echo -e "${GREEN}To activate this build for testing, run:${NC}"
    echo -e "${YELLOW}source ${BUILD_DIR}/venv_py3/bin/activate${NC}"
    echo -e "${YELLOW}export PYTHONPATH=${BUILD_DIR}/build_py3/lib/python3*${NC}"
else
    echo -e "${RED}Build failed. Please check the log at ${BUILD_DIR}/build_py3.log${NC}"
    exit 1
fi
