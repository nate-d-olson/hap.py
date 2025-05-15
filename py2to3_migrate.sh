#!/bin/bash
# py2to3_migrate.sh - Script to automate the migration of hap.py from Python 2 to Python 3
#
# This script:
# 1. Creates backup copies of original Python 2 files
# 2. Prepares Python 3 versions of key modules
# 3. Helps test the Python 3 compatibility

set -e

# Set directory variables
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
BACKUP_DIR="${SCRIPT_DIR}/py2_backup"
PYTHON_DIR="${SCRIPT_DIR}/src/python"
PY3_SUFFIX=".py3"

# Create backup directory if it doesn't exist
mkdir -p "${BACKUP_DIR}"

# Print colored status messages
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${BLUE}=== hap.py Python 2 to 3 Migration Script ===${NC}"

# Function to backup a file before migrating
backup_file() {
    local file=$1
    local backup="${BACKUP_DIR}/$(basename ${file})"
    
    # Create subdirectories in backup dir if needed
    local rel_path=$(realpath --relative-to="${SCRIPT_DIR}" "$(dirname "${file}")")
    mkdir -p "${BACKUP_DIR}/${rel_path}"
    backup="${BACKUP_DIR}/${rel_path}/$(basename ${file})"
    
    if [ ! -f "${backup}" ]; then
        echo -e "${YELLOW}Backing up: ${file} -> ${backup}${NC}"
        cp -p "${file}" "${backup}"
    else
        echo -e "${YELLOW}Backup already exists for: ${file}${NC}"
    fi
}

# Function to migrate a single file
migrate_file() {
    local file=$1
    local py3_file="${file}${PY3_SUFFIX}"
    
    # Check if Python 3 version exists
    if [ ! -f "${py3_file}" ]; then
        echo -e "${RED}Error: Python 3 version not found: ${py3_file}${NC}"
        return 1
    fi
    
    # Backup original file
    backup_file "${file}"
    
    # Replace original with Python 3 version
    echo -e "${GREEN}Migrating: ${file} (Python 2 -> 3)${NC}"
    cp -p "${py3_file}" "${file}"
    
    return 0
}

# Function to migrate a directory
migrate_directory() {
    local dir=$1
    local count=0
    local errors=0
    
    echo -e "${BLUE}Processing directory: ${dir}${NC}"
    
    # Find all Python 3 versions and migrate their corresponding Python 2 files
    for py3_file in $(find "${dir}" -name "*${PY3_SUFFIX}" -type f); do
        # Get the original file name by removing the .py3 suffix
        original_file="${py3_file%${PY3_SUFFIX}}"
        
        if migrate_file "${original_file}"; then
            count=$((count+1))
        else
            errors=$((errors+1))
        fi
    done
    
    echo -e "${GREEN}Successfully migrated ${count} files in ${dir}${NC}"
    if [ ${errors} -gt 0 ]; then
        echo -e "${RED}Encountered ${errors} errors${NC}"
    fi
}

# Function to set up Cython integration
setup_cython_integration() {
    echo -e "${BLUE}Setting up Cython integration for Python 3...${NC}"
    
    # Create Haplo/cython directory if it doesn't exist
    local cython_dir="${PYTHON_DIR}/Haplo/cython"
    mkdir -p "${cython_dir}"
    
    # Check if we've already created Python 3 versions of the Cython files
    if [ -f "${cython_dir}/_internal.pyx" ]; then
        echo -e "${GREEN}Cython integration files already exist.${NC}"
    else
        # Ensure the directory exists and is initialized as a Python package
        echo -e "${GREEN}Creating Cython integration files...${NC}"
        
        # Create __init__.py if it doesn't exist
        if [ ! -f "${cython_dir}/__init__.py" ]; then
            cat > "${cython_dir}/__init__.py" << EOF
# __init__.py for the Haplo.cython package in Python 3
# This file makes the directory a proper Python package

"""Cython extensions for hap.py"""

# Re-export functions from _internal
try:
    from ._internal import *  # noqa
except ImportError:
    import logging
    logging.warning("Failed to import Cython extension modules. Using fallback implementation.")
    
    # Provide fallback implementations of key functionality
    def get_version():
        return "0.0.0-fallback"
    
    def get_build_time():
        return "unknown"
    
    def test_module():
        return {"version": get_version(), "build_time": get_build_time()}
EOF
        fi
        
        # Create a README.md for documentation
        if [ ! -f "${cython_dir}/README.md" ]; then
            cat > "${cython_dir}/README.md" << EOF
# Haplo Cython Module

This directory contains the Cython modules for the hap.py Python/C++ integration.

## Python 3 Migration Notes

In the Python 3 migration, we've updated the Cython integration to:
1. Use Python 3 string handling (Unicode vs bytes)
2. Support modern Cython language features
3. Provide proper error handling and fallbacks
4. Use Python 3's improved memory management
EOF
        fi
        
        # Add test script for Cython integration
        if [ ! -f "${SCRIPT_DIR}/test_cython_integration.py" ]; then
            echo -e "${GREEN}Creating Cython integration test script...${NC}"
            cat > "${SCRIPT_DIR}/test_cython_integration.py" << EOF
#!/usr/bin/env python3
"""
Test script for Cython integration in hap.py Python 3 migration

This script helps diagnose and verify the integration between Python and C++ 
components in the Python 3 version of hap.py.

Usage:
    python3 test_cython_integration.py [--mock]
"""

import argparse
import os
import sys
import logging
from pathlib import Path

# Configure logging
logging.basicConfig(level=logging.INFO, 
                   format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger('cython-test')

def main():
    """Main entry point"""
    parser = argparse.ArgumentParser(description="Test Cython integration")
    parser.add_argument("--mock", action="store_true", 
                       help="Force use of mock implementation")
    args = parser.parse_args()
    
    # Set environment variable if mock mode is requested
    if args.mock:
        os.environ["HAPLO_USE_MOCK"] = "1"
        logger.info("Using mock implementation")
    
    # Add src/python to Python path if needed
    script_dir = Path(__file__).resolve().parent
    python_dir = script_dir / "src" / "python"
    if python_dir.exists():
        sys.path.insert(0, str(python_dir))
        logger.info(f"Added {python_dir} to Python path")
    
    # Test the import
    try:
        logger.info("Importing Haplo...")
        import Haplo
        logger.info(f"Successfully imported Haplo version: {Haplo.__version__}")
        
        # Try to access the Cython module
        try:
            from Haplo.cython import test_module
            result = test_module()
            logger.info(f"Cython module test result: {result}")
        except (ImportError, AttributeError) as e:
            logger.error(f"Failed to import or use Cython module: {e}")
        
        # Test additional components if needed
        logger.info("Testing specific functionality...")
        if hasattr(Haplo, "quantify"):
            logger.info("Found quantify module")
            
    except ImportError as e:
        logger.error(f"Failed to import Haplo: {e}")
        return 1
    
    return 0

if __name__ == "__main__":
    sys.exit(main())
EOF
            chmod +x "${SCRIPT_DIR}/test_cython_integration.py"
        fi
    fi
    
    echo -e "${GREEN}Cython integration setup complete.${NC}"
}

# Check if we should migrate everything or selective modules
if [ "$1" == "--all" ]; then
    echo -e "${BLUE}Migrating all Python files...${NC}"
    
    # Migrate Tools directory
    migrate_directory "${PYTHON_DIR}/Tools"
    
    # Migrate Haplo directory
    migrate_directory "${PYTHON_DIR}/Haplo"
    
    # Migrate top-level Python files
    migrate_directory "${PYTHON_DIR}"
    
    # Check if we need to update the main scripts
    if [ -f "${SCRIPT_DIR}/install_py3.py" ]; then
        backup_file "${SCRIPT_DIR}/install.py"
        echo -e "${GREEN}Migrating: install.py (Python 2 -> 3)${NC}"
        cp -p "${SCRIPT_DIR}/install_py3.py" "${SCRIPT_DIR}/install.py"
    fi
    
    # Apply CMakeLists updates
    if [ -f "${SCRIPT_DIR}/CMakeLists.txt.py3" ]; then
        backup_file "${SCRIPT_DIR}/CMakeLists.txt"
        echo -e "${GREEN}Migrating: CMakeLists.txt (Python 2 -> 3)${NC}"
        cp -p "${SCRIPT_DIR}/CMakeLists.txt.py3" "${SCRIPT_DIR}/CMakeLists.txt"
    fi
    
    # Apply external dependencies script updates
    if [ -f "${SCRIPT_DIR}/external/make_dependencies_py3.sh" ]; then
        backup_file "${SCRIPT_DIR}/external/make_dependencies.sh"
        echo -e "${GREEN}Migrating: external/make_dependencies.sh (Python 2 -> 3)${NC}"
        cp -p "${SCRIPT_DIR}/external/make_dependencies_py3.sh" "${SCRIPT_DIR}/external/make_dependencies.sh"
    fi
    
    # Apply requirements updates
    if [ -f "${SCRIPT_DIR}/happy.requirements.py3.txt" ]; then
        backup_file "${SCRIPT_DIR}/happy.requirements.txt"
        echo -e "${GREEN}Migrating: happy.requirements.txt (Python 2 -> 3)${NC}"
        cp -p "${SCRIPT_DIR}/happy.requirements.py3.txt" "${SCRIPT_DIR}/happy.requirements.txt"
    fi
    
    # Set up Cython integration
    setup_cython_integration
    
elif [ "$1" == "--core" ]; then
    echo -e "${BLUE}Migrating core modules only...${NC}"
    
    # List of core files to migrate
    core_files=(
        "${PYTHON_DIR}/Tools/__init__.py"
        "${PYTHON_DIR}/Tools/bcftools.py"
        "${PYTHON_DIR}/Tools/fastasize.py"
        "${PYTHON_DIR}/Haplo/__init__.py"
        "${PYTHON_DIR}/Haplo/quantify.py"
        "${PYTHON_DIR}/pre.py"
        "${PYTHON_DIR}/hap.py"
        "${PYTHON_DIR}/qfy.py"
    )
    
    for file in "${core_files[@]}"; do
        migrate_file "${file}" || true
    done
    
    # Apply minimal build system updates
    if [ -f "${SCRIPT_DIR}/happy.requirements.py3.txt" ]; then
        backup_file "${SCRIPT_DIR}/happy.requirements.txt"
        echo -e "${GREEN}Migrating: happy.requirements.txt (Python 2 -> 3)${NC}"
        cp -p "${SCRIPT_DIR}/happy.requirements.py3.txt" "${SCRIPT_DIR}/happy.requirements.txt"
    fi
    
    # Set up the Cython integration
    setup_cython_integration
    
    # Set up Cython integration
    setup_cython_integration
    
else
    # Print usage
    echo -e "${YELLOW}Usage: $0 [--all|--core]${NC}"
    echo "  --all   : Migrate all Python 2 files to Python 3"
    echo "  --core  : Migrate only core modules needed for basic functionality"
    echo ""
    echo "This script helps migrate hap.py from Python 2 to Python 3."
    echo "It creates backups of original files in ${BACKUP_DIR} before making changes."
    exit 1
fi

echo -e "${GREEN}Migration complete!${NC}"
echo -e "${YELLOW}To test the migrated code, run:${NC}"
echo "python3 test_python3_compatibility.py"
echo ""
echo -e "${YELLOW}To run the build and installation with Python 3:${NC}"
echo "python3 install.py [build_dir]"
