#!/usr/bin/env python3
"""
Script to automatically fix remaining Python 3 compatibility issues in hap.py.

This script targets:
1. String/unicode issues in remaining files
2. Exception handling issues (except X, e -> except X as e)
3. Subprocess handling improvements

Usage:
    python3 fix_remaining_py3_issues.py [--apply] [--files file1,file2,...]
"""

import argparse
import logging
import os
import re
import sys

# Configure logging
logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

# List of files with known issues
PROBLEMATIC_FILES = [
    "src/python/hap.py",
    "src/python/cnx.py",
    "src/python/ftx.py",
    "src/python/pre.py",
    "src/python/qfy.py",
    "src/python/Tools/bcftools.py",
    "src/python/Tools/__init__.py",
    "src/python/Tools/vcfextract.py",
    "src/python/Haplo/blocksplit.py",
    "src/python/Haplo/cython_compat.py",
    "src/python/Haplo/quantify.py",
    "src/python/Haplo/cython/__init__.py",
]


def fix_exception_syntax(content):
    """Fix Python 2 style exception handling (except X, e: -> except X as e:)"""
    # Pattern to match Python 2 exception handling
    pattern = r"except\s+([A-Za-z0-9_]+)\s*,\s*([A-Za-z0-9_]+)\s*:"
    replacement = r"except \1 as \2:"

    # Use re.sub to replace all occurrences
    new_content = re.sub(pattern, replacement, content)

    return new_content


def fix_file_open(content):
    """Fix file open calls to include encoding parameter"""
    # Pattern to match file open without encoding (text mode)
    open_text_pattern = r'open\(\s*([^,]+)\s*,\s*[\'"]([rwat])[\'"](\s*\))'

    # Use a function to handle replacement to preserve formatting
    def open_replacement(match):
        return f'open({match.group(1)}, "{match.group(2)}", encoding="utf-8"{match.group(3)}'

    # Replace text mode opens
    new_content = re.sub(open_text_pattern, open_replacement, content)

    # Pattern for read-only opens (common case)
    open_read_pattern = r'open\(\s*([^,]+)\s*,\s*[\'"]r[\'"]\s*\)'

    # Replace read mode opens
    new_content = re.sub(
        open_read_pattern, r'open(\1, "r", encoding="utf-8")', new_content
    )

    return new_content


def fix_subprocess_calls(content):
    """Fix subprocess calls to use text mode in Python 3"""
    # Pattern to match subprocess.Popen without universal_newlines
    popen_pattern = r"subprocess\.Popen\((.*?)(stdout\s*=\s*subprocess\.PIPE|stderr\s*=\s*subprocess\.PIPE)(.*?)(?:,\s*universal_newlines\s*=\s*True)?\)"

    def popen_replacement(match):
        if "universal_newlines" in match.group(0):
            return match.group(0)  # Already has universal_newlines
        else:
            return f"subprocess.Popen({match.group(1)}{match.group(2)}{match.group(3)}, universal_newlines=True)"

    # Replace subprocess.Popen calls
    new_content = re.sub(popen_pattern, popen_replacement, content)

    # Pattern to match subprocess.check_output without universal_newlines
    check_output_pattern = (
        r"subprocess\.check_output\((.*?)(?:,\s*universal_newlines\s*=\s*True)?\)"
    )

    def check_output_replacement(match):
        if "universal_newlines" in match.group(0):
            return match.group(0)  # Already has universal_newlines
        else:
            return f"subprocess.check_output({match.group(1)}, universal_newlines=True)"

    # Replace subprocess.check_output calls
    new_content = re.sub(check_output_pattern, check_output_replacement, new_content)

    return new_content


def fix_string_handling(content):
    """Fix string handling issues for Python 3 compatibility"""
    # Replace u'' string literals (less critical in Python 3.3+)
    unicode_pattern = r'u([\'"])'
    new_content = re.sub(unicode_pattern, r"\1", content)

    # Find potential str() calls that might need explicit encoding handling
    # This is more complex and might need manual review

    def str_replacement(match):
        var_name = match.group(1).strip()
        # Only replace in contexts where binary data is likely
        if (
            "subprocess" in var_name
            or "read(" in var_name
            or "stdout" in var_name
            or "stderr" in var_name
            or "Exception" in var_name
        ):
            return f"{var_name}.decode('utf-8') if isinstance({var_name}, bytes) else str({var_name})"
        return match.group(0)

    # This is a higher-risk transformation, disabled by default
    # new_content = re.sub(str_pattern, str_replacement, new_content)

    return new_content


def fix_python3_issues(file_path, apply=False, verbose=False):
    """Fix Python 3 compatibility issues in the given file"""
    try:
        with open(file_path, encoding="utf-8") as f:
            content = f.read()

        logger.info(f"Processing {file_path}...")

        # Make a copy of the original content to detect changes
        original_content = content

        # Apply fixes
        content = fix_exception_syntax(content)
        content = fix_file_open(content)
        content = fix_subprocess_calls(content)
        content = fix_string_handling(content)

        # Check if content has changed
        if content != original_content:
            if apply:
                with open(file_path, "w", encoding="utf-8") as f:
                    f.write(content)
                logger.info(f"  Fixed Python 3 compatibility issues in {file_path}")
            else:
                logger.info(f"  Found issues in {file_path} (run with --apply to fix)")
            return True
        else:
            if verbose:
                logger.info(f"  No issues found in {file_path}")
            return False

    except Exception as e:
        logger.error(f"Error processing {file_path}: {e}")
        return False


def main():
    """Main entry point"""
    parser = argparse.ArgumentParser(
        description="Fix remaining Python 3 compatibility issues in hap.py"
    )
    parser.add_argument("--apply", action="store_true", help="Apply fixes to files")
    parser.add_argument(
        "--verbose", "-v", action="store_true", help="Show verbose output"
    )
    parser.add_argument(
        "--files",
        help="Comma-separated list of specific files to process (default: all known problematic files)",
    )
    args = parser.parse_args()

    files_to_process = PROBLEMATIC_FILES
    if args.files:
        files_to_process = args.files.split(",")

    fixed_files = 0
    for file_path in files_to_process:
        if os.path.exists(file_path):
            if fix_python3_issues(file_path, args.apply, args.verbose):
                fixed_files += 1
        else:
            logger.warning(f"File not found: {file_path}")

    logger.info(f"\nSummary: Found issues in {fixed_files} files")
    if args.apply:
        logger.info("Applied fixes to all files")
    else:
        logger.info("Run with --apply to fix issues automatically")

    return 0


if __name__ == "__main__":
    sys.exit(main())
