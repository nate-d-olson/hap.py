#!/usr/bin/env python3
"""
Update Cython modules for Python 3 compatibility in hap.py

This script creates or updates Cython module files (.pyx, .pxd) to ensure
they're compatible with Python 3. It sets the proper language level directives
and updates string handling to work with Python 3's Unicode strings.

Usage:
    python update_cython_modules_py3.py --src-dir /path/to/src/dir
"""

import argparse
import re
import shutil
import sys
from pathlib import Path


def process_pyx_file(file_path: Path, backup: bool = True) -> bool:
    """Process a .pyx file to ensure Python 3 compatibility

    Args:
        file_path: Path to the .pyx file
        backup: Whether to create a backup of the original file

    Returns:
        True if changes were made, False otherwise
    """
    if not file_path.exists():
        print(f"File {file_path} does not exist")
        return False

    # Read the file
    with open(file_path, encoding="utf-8") as f:
        content = f.read()

    original_content = content

    # Check if language level directive already exists
    lang_level_pattern = re.compile(r"#\s*cython:\s*language_level\s*=\s*\d+")

    if not lang_level_pattern.search(content):
        # Add language level directive at the top of the file
        content = "# cython: language_level=3\n# distutils: language=c++\n" + content
        print(f"Added language level directive to {file_path}")

    # Fix string handling patterns
    string_fixes = [
        # Python 2 style string to bytes comparison
        (
            r"([^\.])cdef\s+char\s*\*\s*([a-zA-Z0-9_]+)\s*=\s*([a-zA-Z0-9_\.]+)",
            r"\1cdef char *\2 = \3.encode('utf-8')",
        ),
        # Python 2 style direct assignment from Python string to char*
        (
            r"([^\.])(\s*)([a-zA-Z0-9_]+)\s*=\s*<char\s*\*>\s*([a-zA-Z0-9_\.]+)",
            r"\1\2\3 = <char *>\4.encode('utf-8')",
        ),
        # String returns
        (
            r"([^\.])return\s+<char\s*\*>\s*([a-zA-Z0-9_\.]+)",
            r"\1cdef bytes py_bytes = \2.encode('utf-8')\n    return <char *>py_bytes",
        ),
        # C string to Python string conversion
        (
            r"([^\.])return\s+([a-zA-Z0-9_]+)",
            r"\1return \2.decode('utf-8') if isinstance(\2, bytes) else \2",
        ),
    ]

    for pattern, replacement in string_fixes:
        content = re.sub(pattern, replacement, content)

    # Check if changes were made
    if content != original_content:
        if backup:
            backup_path = f"{file_path}.py2.bak"
            print(f"Creating backup of {file_path} to {backup_path}")
            shutil.copy2(file_path, backup_path)

        # Write the updated content
        with open(file_path, "w", encoding="utf-8") as f:
            f.write(content)

        return True

    return False


def process_pxd_file(file_path: Path, backup: bool = True) -> bool:
    """Process a .pxd file to ensure Python 3 compatibility

    Args:
        file_path: Path to the .pxd file
        backup: Whether to create a backup of the original file

    Returns:
        True if changes were made, False otherwise
    """
    if not file_path.exists():
        print(f"File {file_path} does not exist")
        return False

    # Read the file
    with open(file_path, encoding="utf-8") as f:
        content = f.read()

    original_content = content

    # Add proper imports if needed
    if "from libcpp.string cimport string" not in content:
        # Add at the top after any cython imports
        import_pos = content.find("cimport")
        if import_pos >= 0:
            # Find the end of the import section
            next_line_pos = content.find("\n", import_pos)
            if next_line_pos >= 0:
                content = (
                    content[: next_line_pos + 1]
                    + "from libcpp.string cimport string\n"
                    + content[next_line_pos + 1 :]
                )
        else:
            # Add at the top
            content = "from libcpp.string cimport string\n" + content

        print(f"Added string import to {file_path}")

    # Fix char* declarations to use string when appropriate
    char_ptr_pattern = re.compile(r"(cdef\s+)char\s*\*\s*([a-zA-Z0-9_]+)")
    content = char_ptr_pattern.sub(r"\1string \2", content)

    # Check if changes were made
    if content != original_content:
        if backup:
            backup_path = f"{file_path}.py2.bak"
            print(f"Creating backup of {file_path} to {backup_path}")
            shutil.copy2(file_path, backup_path)

        # Write the updated content
        with open(file_path, "w", encoding="utf-8") as f:
            f.write(content)

        return True

    return False


def create_mock_implementation(module_path: Path) -> bool:
    """Create a mock implementation for a Cython module

    Args:
        module_path: Path to the Cython module (.pyx file)

    Returns:
        True if mock was created, False otherwise
    """
    if not module_path.exists():
        print(f"File {module_path} does not exist")
        return False

    # Extract module name
    module_name = module_path.stem

    # Create mock file path
    mock_path = module_path.parent / f"mock_{module_name}.py"

    if mock_path.exists():
        print(f"Mock implementation {mock_path} already exists")
        return True

    # Basic mock template
    mock_content = f'''"""
Mock implementation of {module_name} for Python 3 testing without C++ components

This module provides fallback functionality when the Cython module
cannot be imported or compiled.
"""

import logging
import warnings

# Log warning about using mock implementation
warnings.warn(f"Using mock implementation of {module_name}")
logging.warning(f"Using mock implementation of {module_name}")


def test_string_handling():
    """Mock test for string handling"""
    return "Mock implementation working"


def test_basic_functionality():
    """Mock test for basic functionality"""
    return "Mock basic functionality working"


# Add mock implementations of the functions in the Cython module
# TODO: Analyze the original module and add appropriate mock functions
'''

    # Write the mock implementation
    with open(mock_path, "w", encoding="utf-8") as f:
        f.write(mock_content)

    print(f"Created mock implementation {mock_path}")
    return True


def create_module_init(directory: Path) -> bool:
    """Create a Python 3 compatible __init__.py file with fallback mechanism

    Args:
        directory: Path to the directory where the __init__.py file should be created

    Returns:
        True if file was created or updated, False otherwise
    """
    init_path = directory / "__init__.py"

    init_content = '''"""
Python 3 compatible Cython module package with fallback mechanism

This package tries to import the Cython modules, but falls back to
mock implementations if they are not available.
"""

import os
import sys
import warnings
import logging

# Flag to control whether to use mock implementation
USE_MOCK = os.environ.get("HAPLO_USE_MOCK", "0").lower() in ("1", "true", "yes")

try:
    if USE_MOCK:
        # Force using mock implementation
        raise ImportError("Using mock implementation as requested")

    # Try to import the Cython modules
    from ._internal import *
    from .cpp_internal import *

    # Add any other Cython modules here

    # Signal that we're using the real implementation
    USING_MOCK = False

except ImportError as e:
    warnings.warn(f"Failed to import Cython modules: {e}. Using mock implementation.")
    logging.warning(f"Failed to import Cython modules: {e}. Using mock implementation.")

    # Fall back to mock implementation
    from .mock_internal import *
    # Add other mock imports here

    # Signal that we're using the mock implementation
    USING_MOCK = True

# Version information
__version__ = "0.1.0"
'''

    # Write or update the init file
    with open(init_path, "w", encoding="utf-8") as f:
        f.write(init_content)

    print(f"Created/updated {init_path}")
    return True


def update_cython_modules(src_dir: str, dry_run: bool = False) -> None:
    """Update Cython modules to be Python 3 compatible

    Args:
        src_dir: Path to the source directory
        dry_run: If True, don't make any changes
    """
    src_path = Path(src_dir).resolve()

    if not src_path.exists():
        print(f"Source directory {src_path} does not exist")
        sys.exit(1)

    print(f"Processing Cython modules in {src_path}")

    # Find all Cython files
    pyx_files = list(src_path.glob("**/*.pyx"))
    pxd_files = list(src_path.glob("**/*.pxd"))

    print(f"Found {len(pyx_files)} .pyx files and {len(pxd_files)} .pxd files")

    # Process files
    if dry_run:
        print("Dry run - not making any changes")
    else:
        # Keep track of directories with Cython modules
        cython_dirs = set()

        # Process .pyx files
        for pyx_file in pyx_files:
            cython_dirs.add(pyx_file.parent)
            if process_pyx_file(pyx_file):
                print(f"Updated {pyx_file}")
                # Create mock implementation
                create_mock_implementation(pyx_file)

        # Process .pxd files
        for pxd_file in pxd_files:
            cython_dirs.add(pxd_file.parent)
            if process_pxd_file(pxd_file):
                print(f"Updated {pxd_file}")

        # Create __init__.py files for each directory
        for cython_dir in cython_dirs:
            create_module_init(cython_dir)


def main():
    parser = argparse.ArgumentParser(
        description="Update Cython modules for Python 3 compatibility"
    )
    parser.add_argument("--src-dir", required=True, help="Path to source directory")
    parser.add_argument("--dry-run", action="store_true", help="Don't make any changes")

    args = parser.parse_args()

    update_cython_modules(args.src_dir, args.dry_run)


if __name__ == "__main__":
    main()
