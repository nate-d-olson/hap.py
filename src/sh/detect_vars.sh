#!/usr/bin/env bash
# detect_vars.sh - simplified environment detection for run_tests.sh

# If HCDIR is already set (e.g., by run_tests.sh), do nothing
if [[ -n "${HCDIR:-}" ]]; then
    return 0
fi

# If hap.py is on PATH, use that for HCDIR
if command -v hap.py &> /dev/null; then
    HCDIR="$(dirname "$(command -v hap.py)")"
    export HCDIR
    return 0
fi

# Otherwise, cannot proceed
echo "Cannot find hap.py binary. Please install package first." >&2
exit 1
