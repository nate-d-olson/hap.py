#!/usr/bin/env python3
"""Simple debug script for ftx.py."""

import os
import subprocess
import sys
from pathlib import Path

print("Starting ftx.py debug test...")

# Get the path to the ftx.py script
script_dir = Path(__file__).resolve().parent.parent / "src" / "python"
ftx_script = script_dir / "ftx.py"

print(f"ftx.py script path: {ftx_script}")
print(f"Script exists: {ftx_script.exists()}")

# Run ftx.py with --help
cmd = [sys.executable, str(ftx_script), "--help"]
print(f"Running command: {' '.join(cmd)}")

env = os.environ.copy()
env["HAPLO_USE_MOCK"] = "1"

try:
    result = subprocess.run(
        cmd,
        env=env,
        capture_output=True,
        text=True,
        check=False,
    )

    print(f"Exit code: {result.returncode}")
    print(f"STDOUT: {result.stdout}")
    print(f"STDERR: {result.stderr}")
except Exception as e:
    print(f"Error running ftx.py: {str(e)}")

print("Debug test completed.")
