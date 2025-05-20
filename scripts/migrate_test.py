#!/usr/bin/env python3
# filepath: /Users/nolson/hap.py-update-take2/hap.py/scripts/migrate_test.py
"""
Tool to help migrate shell tests to pytest.

Usage:
    python scripts/migrate_test.py src/sh/run_test_name.sh tests/integration/test_name.py
"""

import os
import re
import sys


def extract_test_info(shell_script_path):
    """Extract information from a shell test script.

    Args:
        shell_script_path: Path to the shell script

    Returns:
        dict with test information
    """
    with open(shell_script_path) as f:
        content = f.read()

    # Extract test name from filename
    filename = os.path.basename(shell_script_path)
    if filename.startswith("run_") and filename.endswith(".sh"):
        test_name = filename[4:-3]  # Remove 'run_' and '.sh'
    else:
        test_name = filename[:-3]  # Just remove '.sh'

    # Try to extract description (first echo line)
    description_match = re.search(r'echo\s+"([^"]+)"', content)
    description = (
        description_match.group(1)
        if description_match
        else f"Test {test_name} functionality."
    )

    # Find command executions
    commands = []
    # Look for typical command patterns in shell scripts
    command_patterns = [
        r"\$\{HCDIR\}/([^\s]+)([^\n]+)",
        r"\$\{PYTHON\}\s+([^\n]+)",
        r"bin/([^\s]+)([^\n]+)",
        r"\$\{DIR\}/../[^/]+/([^\s]+)([^\n]+)",  # Additional pattern
    ]

    for pattern in command_patterns:
        for match in re.finditer(pattern, content):
            cmd = match.group(0).strip()
            if cmd and not cmd.startswith("#") and "|| exit 1" not in cmd:
                commands.append(cmd)

    # Extract input files
    input_files = []
    file_patterns = [
        r"(\$\{DIR\}/[^\s]+\.vcf|\$\{DIR\}/[^\s]+\.bed|\$\{DIR\}/[^\s]+\.fa)",
        r"(\$\{DIR\}/../[^\s]+\.vcf|\$\{DIR\}/../[^\s]+\.bed|\$\{DIR\}/../[^\s]+\.fa)",
        r"(example/[^\s]+\.vcf|example/[^\s]+\.bed|example/[^\s]+\.fa)",
        r"(\$\{[A-Z_]+\}/[^\s]+\.vcf|\$\{[A-Z_]+\}/[^\s]+\.bed|\$\{[A-Z_]+\}/[^\s]+\.fa)",
    ]

    for pattern in file_patterns:
        for match in re.finditer(pattern, content):
            file_path = match.group(1).strip()
            if file_path and not file_path.endswith("> /dev/null"):
                input_files.append(file_path)

    return {
        "name": test_name,
        "description": description,
        "commands": commands,
        "input_files": input_files,
    }


def generate_pytest_template(test_info, output_path):
    """Generate a pytest template from test information.

    Args:
        test_info: Dict with test information
        output_path: Where to write the output file
    """
    template = f'''"""
Integration tests for {test_info["name"]} functionality.
Migrated from src/sh/run_{test_info["name"]}.sh
"""
import os
import subprocess
import pytest
from pathlib import Path
import sys
import tempfile

from tests.utils import (
    get_project_root,
    get_bin_dir,
    get_example_dir,
    run_command,
    get_python_executable,
    find_reference_file
)


@pytest.mark.integration
def test_{test_info["name"]}(tmp_path):
    """Test {test_info["description"]}"""
    # Get paths to required files and tools
    project_root = get_project_root()
    bin_dir = get_bin_dir()
    example_dir = get_example_dir()
    python_exe = get_python_executable()

    # TODO: Replace with actual paths based on the shell script
    # Example input files:
    {chr(10).join([f'# {file}' for file in test_info["input_files"]])}

    # Example commands:
    {chr(10).join([f'# {cmd}' for cmd in test_info["commands"]])}

    # Implement the test based on the shell script logic
    assert True  # Replace with actual assertions
'''

    with open(output_path, "w") as f:
        f.write(template)

    print(f"Generated pytest template at {output_path}")
    print("Please review and complete the implementation based on the shell script.")


def main():
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <shell_script_path> <output_pytest_path>")
        sys.exit(1)

    shell_script_path = sys.argv[1]
    output_path = sys.argv[2]

    if not os.path.exists(shell_script_path):
        print(f"Error: Shell script {shell_script_path} not found.")
        sys.exit(1)

    test_info = extract_test_info(shell_script_path)
    generate_pytest_template(test_info, output_path)


if __name__ == "__main__":
    main()
