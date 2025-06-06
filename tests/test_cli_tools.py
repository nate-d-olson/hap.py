#!/usr/bin/env python3
"""Test script to verify the CLI entry points are working correctly."""

import os
import subprocess
import sys
from pathlib import Path
import pytest


COMMAND_TO_MODULE = {
    "hap.py": "hap_py.hap",
    "quantify": "hap_py.qfy",
    "preprocess": "hap_py.pre",
}


def _cli_tool(command, args=None, expected_exit_code=0):
    """Test a CLI tool with given arguments and verify the exit code."""
    if args is None:
        args = ["--help"]

    module_name = COMMAND_TO_MODULE.get(command, command)

    try:
        env = os.environ.copy()
        env["PYTHONPATH"] = str(Path(__file__).resolve().parent.parent / "src")
        result = subprocess.run(
            [sys.executable, "-m", module_name] + args,
            capture_output=True,
            text=True,
            check=False,
            env=env,
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


@pytest.mark.parametrize("command", list(COMMAND_TO_MODULE.keys()))
def test_cli_tool_help(command):
    assert _cli_tool(command)


@pytest.mark.parametrize("command", list(COMMAND_TO_MODULE.keys()))
def test_cli_tool_version(command):
    assert _cli_tool(command, ["--version"])


def test_cli_tool_check_deps():
    assert _cli_tool("hap.py", ["--check-deps"])


def test_cli_tool_invalid_option():
    assert _cli_tool("hap.py", ["--invalid-option"], expected_exit_code=1)


def run_all_tests():
    """Run tests for all CLI tools."""
    success = True

    # Test each CLI tool with --help
    for command in COMMAND_TO_MODULE:
        if not _cli_tool(command):
            success = False

    # Test each CLI tool with --version
    for command in COMMAND_TO_MODULE:
        if not _cli_tool(command, ["--version"]):
            success = False

    # Check dependency reporting
    if not _cli_tool("hap.py", ["--check-deps"]):
        success = False

    # Test with invalid arguments (should fail with non-zero exit code)
    if not _cli_tool("hap.py", ["--invalid-option"], expected_exit_code=1):
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
