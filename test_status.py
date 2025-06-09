#!/usr/bin/env python3
"""
Test script to check the current status of the hap.py package
"""
import os
import subprocess
import sys

print("=== HAP.PY STATUS CHECK ===")
print(f"Python version: {sys.version}")
print(f"Current directory: {os.getcwd()}")
print(f"Python path: {sys.path}")

# Check if we can import hap_py
try:
    import hap_py

    print(f"✓ hap_py package found at: {hap_py.__file__}")
    try:
        print(f"✓ hap_py version: {hap_py.__version__}")
    except AttributeError:
        print("⚠ hap_py version not available")
except ImportError as e:
    print(f"✗ Failed to import hap_py: {e}")

# Check if pytest is available
try:
    import pytest

    print(f"✓ pytest available: {pytest.__version__}")
except ImportError:
    print("✗ pytest not available")

# Check RTG tools
rtg_path = "/Users/nolson/hap.py-modern-claude4/hap.py/external/rtg-tools-3.12.1/rtg"
if os.path.exists(rtg_path):
    print(f"✓ RTG tools found at: {rtg_path}")
else:
    print(f"✗ RTG tools not found at: {rtg_path}")

# Check if we can run basic commands
try:
    result = subprocess.run(
        [sys.executable, "--version"], capture_output=True, text=True
    )
    print(f"✓ Python executable works: {result.stdout.strip()}")
except Exception as e:
    print(f"✗ Python executable test failed: {e}")
