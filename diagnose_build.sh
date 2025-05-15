#!/bin/bash
# diagnose_build.sh - Script to diagnose build issues in hap.py
#
# Usage: ./diagnose_build.sh [build_dir]
#
# This script runs the installation with verbose logging, captures
# the output, and helps identify common build errors.

set -e

# Ensure we're in the hap.py directory
cd "$(dirname "$0")"

# Default build directory
BUILD_DIR=${1:-/tmp/happy-build}

# Create logs directory
LOGS_DIR="./build_logs"
mkdir -p "$LOGS_DIR"

echo "===== Starting hap.py build diagnosis ====="
echo "Build directory: $BUILD_DIR"
echo "Logs directory: $LOGS_DIR"

echo "1. Running full build with verbose logging..."
python install.py "$BUILD_DIR" VERBOSE=1 2>&1 | tee "$LOGS_DIR/full_build.log" || true

echo "2. Extracting error information..."
grep -A 10 "error:" "$LOGS_DIR/full_build.log" > "$LOGS_DIR/build_errors.log" || true
grep -A 10 "Error " "$LOGS_DIR/full_build.log" >> "$LOGS_DIR/build_errors.log" || true
grep -A 10 "CMake Error" "$LOGS_DIR/full_build.log" >> "$LOGS_DIR/build_errors.log" || true
grep -A 10 "Failed" "$LOGS_DIR/full_build.log" >> "$LOGS_DIR/build_errors.log" || true

echo "3. Checking external dependencies..."
python install.py --build-externals-only "$BUILD_DIR-ext" 2>&1 | tee "$LOGS_DIR/externals_build.log" || true

echo "4. Analyzing missing dependencies..."
grep -A 5 "Could not find" "$LOGS_DIR/full_build.log" > "$LOGS_DIR/missing_deps.log" || true

echo "5. Checking for C++ compilation issues..."
grep -B 2 -A 5 "error:" "$LOGS_DIR/full_build.log" | grep -E '\.cpp|\.hpp|\.h|\.cc' > "$LOGS_DIR/cpp_errors.log" || true

echo "===== Build diagnosis complete ====="
echo "Check the logs in the $LOGS_DIR directory"
echo ""
echo "Common issues summary:"
echo "---------------------"
wc -l "$LOGS_DIR"/*.log | sort -nr

echo ""
echo "Next steps:"
echo "1. Review the build errors in $LOGS_DIR/build_errors.log"
echo "2. Check for missing dependencies in $LOGS_DIR/missing_deps.log"
echo "3. Fix C++ compilation errors listed in $LOGS_DIR/cpp_errors.log"
