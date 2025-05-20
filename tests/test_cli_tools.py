#!/usr/bin/env python3
"""
Test script to verify the CLI entry points are working correctly.

This script tests each of the command-line tools to ensure they
can be called properly and return appropriate exit codes.
"""

import os
import subprocess
import sys
from pathlib import Path


def test_cli_tool(command, args=None, expected_exit_code=0):
    """Test a CLI tool with given arguments and verify the exit code."""
    if args is None:
        args = ["--help"]

    try:
        result = subprocess.run(
            [command] + args,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=False,
        )

        if result.returncode != expected_exit_code:
            print(
                f"ERROR: {command} returned exit code {result.returncode}, expected {expected_exit_code}"
            )
            print(f"STDERR: {result.stderr}")
            return False

        print(f"SUCCESS: {command} {' '.join(args)} (exit code: {result.returncode})")
        return True
    except Exception as e:
        print(f"ERROR: Failed to run {command}: {str(e)}")
        return False


def run_all_tests():
    """Run tests for all CLI tools."""
    success = True

    # Test each CLI tool with --help
    for command in ["hap", "qfy", "pre", "ftx", "ovc", "cnx"]:
        if not test_cli_tool(command):
            success = False

    # Test each CLI tool with --version
    for command in ["hap", "qfy", "pre", "ftx"]:  # ovc and cnx may not have --version
        if not test_cli_tool(command, ["--version"]):
            success = False

    # Test with invalid arguments (should fail with non-zero exit code)
    if not test_cli_tool("hap", ["--invalid-option"], expected_exit_code=1):
        success = False

    return success


if __name__ == "__main__":
    print("Testing hap.py CLI tools...")
    success = run_all_tests()

    if success:
        print("\nAll CLI tests passed!")
        sys.exit(0)
    else:
        print("\nSome CLI tests failed!")
        sys.exit(1)
