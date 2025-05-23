#!/bin/bash
# Simple wrapper for RTG tools that bypasses architecture detection issues
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
exec java -jar "$SCRIPT_DIR/RTG.jar" "$@"