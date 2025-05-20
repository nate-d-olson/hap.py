#!/usr/bin/env python3
"""
Test script to verify the CLI scripts in source directory.

This script tests each of the command-line tools to ensure they
can be called directly and return appropriate exit codes.
"""

import os
import subprocess
import sys
from pathlib import Path
import shutil


def test_cli_script(script_path, args=None, expected_exit_code=0):
    """Test a CLI script with given arguments and verify the exit code."""
    if args is None:
        args = ["--help"]

    print(f"Testing: {script_path} {' '.join(args)}")

    try:
        result = subprocess.run(
            [sys.executable, script_path] + args,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=False,
        )

        print(f"  Exit code: {result.returncode}")
        if result.stdout:
            print(
                f"  STDOUT: {result.stdout[:100]}..."
                if len(result.stdout) > 100
                else f"  STDOUT: {result.stdout}"
            )

        if result.returncode != expected_exit_code:
            print(
                f"ERROR: {script_path} returned exit code {result.returncode}, expected {expected_exit_code}"
            )
            print(f"STDERR: {result.stderr}")
            return False

        print(
            f"SUCCESS: {script_path} {' '.join(args)} (exit code: {result.returncode})"
        )
        return True
    except Exception as e:
        print(f"ERROR: Failed to run {script_path}: {str(e)}")
        return False


def run_all_tests():
    """Run tests for all CLI scripts in the src/python directory."""
    success = True

    # Base path for scripts
    script_dir = Path(__file__).resolve().parent.parent / "src" / "python"

    # List of scripts to test
    script_names = ["hap.py", "qfy.py", "pre.py", "ftx.py", "ovc.py", "cnx.py"]

    print(f"Checking for scripts in {script_dir}...")
    scripts = []
    for script_name in script_names:
        script_path = script_dir / script_name
        if script_path.exists():
            print(f"Found script: {script_path}")
            scripts.append(script_path)
        else:
            print(f"Script not found: {script_path}")

    if not scripts:
        print("ERROR: No scripts found to test!")
        return False

    # Test each CLI script with --help
    for script_path in scripts:
        if not test_cli_script(script_path):
            success = False

    # Test with invalid arguments (should fail with non-zero exit code)
    if (script_dir / "hap.py").exists():
        if not test_cli_script(
            script_dir / "hap.py", ["--invalid-option"], expected_exit_code=1
        ):
            success = False

    return success


if __name__ == "__main__":
    print("Testing hap.py CLI scripts...")
    success = run_all_tests()

    if success:
        print("\nAll CLI script tests passed!")
        sys.exit(0)
    else:
        print("\nSome CLI script tests failed!")
        sys.exit(1)
