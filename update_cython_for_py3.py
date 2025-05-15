#!/usr/bin/env python3
"""
Script to update Cython modules for Python 3 compatibility
"""

import argparse
import os
import re
import sys


def update_cython_file(file_path, apply=False, verbose=False):
    """Update a Cython file for Python 3 compatibility"""
    if not (file_path.endswith(".pyx") or file_path.endswith(".pxd")):
        return False

    try:
        with open(file_path, "r", encoding="utf-8") as f:
            content = f.read()
    except UnicodeDecodeError:
        print(f"Warning: Unable to read {file_path} as UTF-8. Skipping.")
        return False

    # Check if language_level directive is already present
    has_language_level = (
        re.search(r"#\s*cython:\s*language_level\s*=\s*3", content) is not None
    )

    changes_made = False
    new_content = content

    # Add language_level directive if missing
    if not has_language_level:
        if re.match(r"^#\s*distutils:", content):
            # Add after existing distutils directive
            new_content = re.sub(
                r"^(#\s*distutils:.*)$",
                r"\1\n# cython: language_level=3",
                new_content,
                count=1,
                flags=re.MULTILINE,
            )
        else:
            # Add at the top of the file
            new_content = "# cython: language_level=3\n" + new_content
        changes_made = True
        if verbose:
            print(f"  Adding language_level directive to {file_path}")

    # Fix string handling for Python 3
    # Update str/bytes conversions
    str_patterns = [
        (
            r"str\(([^)]+)\)",
            r'\1.decode("utf-8") if isinstance(\1, bytes) else str(\1)',
        ),
        (r"\.encode\(\)", r'.encode("utf-8")'),
        (r"\.decode\(\)", r'.decode("utf-8")'),
        (r"cdef\s+char\s*\*", r"cdef bytes"),  # Replace char* with bytes
        (
            r"string\(\s*<char\*>\s*([^)]+)\)",
            r'string(<char*>(\1 if isinstance(\1, bytes) else \1.encode("utf-8")))',
        ),
    ]

    for pattern, replacement in str_patterns:
        if re.search(pattern, new_content):
            new_content = re.sub(pattern, replacement, new_content)
            changes_made = True
            if verbose:
                print(f"  Fixed string pattern: {pattern} in {file_path}")

    # Fix Python C API functions
    py_api_patterns = [
        (r"PyString_", r"PyBytes_"),
        (r"PyInt_", r"PyLong_"),
    ]

    for pattern, replacement in py_api_patterns:
        if re.search(pattern, new_content):
            new_content = re.sub(pattern, replacement, new_content)
            changes_made = True
            if verbose:
                print(f"  Fixed Python C API pattern: {pattern} in {file_path}")

    if changes_made:
        print(f"Found Python 3 compatibility issues in: {file_path}")

        if apply:
            with open(file_path, "w", encoding="utf-8") as f:
                f.write(new_content)
            print(f"  Applied fixes to {file_path}")
        elif verbose:
            print("  Changes needed:")
            if not has_language_level:
                print("  - Add '# cython: language_level=3' directive")
            for pattern, replacement in str_patterns:
                if re.search(pattern, content):
                    print(f"  - Replace '{pattern}' with '{replacement}'")
            for pattern, replacement in py_api_patterns:
                if re.search(pattern, content):
                    print(f"  - Replace '{pattern}' with '{replacement}'")

    return changes_made


def main():
    parser = argparse.ArgumentParser(
        description="Update Cython modules for Python 3 compatibility"
    )
    parser.add_argument(
        "--path", required=True, help="Path to scan, can be a file or directory"
    )
    parser.add_argument("--apply", action="store_true", help="Apply fixes")
    parser.add_argument(
        "--verbose", "-v", action="store_true", help="Show verbose output"
    )
    args = parser.parse_args()

    path = args.path
    apply = args.apply
    verbose = args.verbose

    files_updated = 0

    if os.path.isfile(path) and (path.endswith(".pyx") or path.endswith(".pxd")):
        if update_cython_file(path, apply, verbose):
            files_updated += 1
    elif os.path.isdir(path):
        for root, _, files in os.walk(path):
            for file in files:
                if file.endswith(".pyx") or file.endswith(".pxd"):
                    file_path = os.path.join(root, file)
                    if update_cython_file(file_path, apply, verbose):
                        files_updated += 1
    else:
        print(f"Error: {path} is not a valid Cython file or directory")
        return 1

    print(
        f"\nSummary: Found {files_updated} Cython files with Python 3 compatibility issues"
    )
    if apply:
        print(f"Applied fixes to {files_updated} files")
    else:
        print(f"Run with --apply to apply fixes")

    return 0


if __name__ == "__main__":
    sys.exit(main())
