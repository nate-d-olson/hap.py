#!/usr/bin/env bash
set -euo pipefail

# Simplified run_tests.sh - preserve legacy functionality

# Determine script directory
DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$DIR/../.." && pwd)"

# Dry-run helper (set DRY_RUN=1 to skip execution)
DRY_RUN=${DRY_RUN:-0}
_run() {
  if [[ "$DRY_RUN" -eq 1 ]]; then
    echo "DRY_RUN: $*"; return 0
  else
    "$@"
  fi
}

# On macOS, ensure CMake uses the correct SDK
if [[ "$(uname)" == "Darwin" ]]; then
  export EXTRA_CMAKE_OPTS="-DCMAKE_OSX_SYSROOT=$(xcrun --show-sdk-path)"
fi

# Install and build C++ and Python components into local directory
INSTALL_DIR="${INSTALL_DIR:-$PWD/.happy-install}"
echo "Installing into ${INSTALL_DIR}"
_run python3 "$ROOT/install.py" "$INSTALL_DIR"

# Prepend binaries and CLI scripts to PATH
export PATH="$INSTALL_DIR/bin:$PATH"
PYTHON="${PYTHON:-python3}"
HCDIR="$INSTALL_DIR/bin"
export PYTHON HCDIR

# Counters
FAILED_TESTS=""
TEST_COUNT=0
FAIL_COUNT=0

# BOOST unit tests
echo "Running test_haplotypes"
_run "$HCDIR/test_haplotypes"
((TEST_COUNT++))
if [[ $? -ne 0 ]]; then
  echo "[FAILED] test_haplotypes"; FAILED_TESTS+="\n- test_haplotypes"; ((FAIL_COUNT++))
else
  echo "[PASSED] test_haplotypes"
fi

# Other tests via helper scripts
SCRIPTS=(
  run_multimerge_test.sh run_hapenum_test.sh run_hapcmp_test.sh
  run_pathtraversal_test.sh run_fp_accuracy_test.sh run_faulty_variant_test.sh
  run_leftshift_test.sh run_other_vcf_tests.sh run_gvcf_homref_test.sh
  run_chrprefix_test.sh run_decomp_test.sh run_integration_test.sh
  run_scmp_test.sh
)
for script in "${SCRIPTS[@]}"; do
  echo "Running ${script}"
  _run bash "$DIR/$script"
  ((TEST_COUNT++))
  if [[ $? -ne 0 ]]; then
    echo "[FAILED] ${script}"; FAILED_TESTS+="\n- ${script}"; ((FAIL_COUNT++))
  else
    echo "[PASSED] ${script}"
  fi
done

# Python-based fastasize test
echo "Running run_fastasize_test.py"
_run "$PYTHON" "$DIR/run_fastasize_test.py"
((TEST_COUNT++))
if [[ $? -ne 0 ]]; then
  echo "[FAILED] run_fastasize_test.py"; FAILED_TESTS+="\n- run_fastasize_test.py"; ((FAIL_COUNT++))
else
  echo "[PASSED] run_fastasize_test.py"
fi

# Summary
echo -e "\nRan ${TEST_COUNT} tests: ${FAIL_COUNT} failures"
if [[ ${FAIL_COUNT} -ne 0 ]]; then
  echo -e "Failures:${FAILED_TESTS}"; exit 1
fi

echo "All tests passed!"