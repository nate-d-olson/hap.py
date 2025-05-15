#!/usr/bin/env python3
"""
Script to specifically find and fix Python 2 style exception syntax (except X, y:)
"""

import argparse
import os
import re
import sys


def find_exception_syntax_issues(file_path):
    """Find Python 2 style exception handling in a file"""
    with open(file_path, "r", encoding="utf-8") as f:
        content = f.read()

    # Pattern to match Python 2 exception handling (except X, e:)
    pattern = r"except\s+([A-Za-z0-9_]+)\s*,\s*([A-Za-z0-9_]+)\s*:"
    issues = []

    for match in re.finditer(pattern, content):
        lineno = content[: match.start()].count("\n") + 1
        exception_class = match.group(1)
        exception_var = match.group(2)
        issues.append(
            {
                "line": lineno,
                "text": match.group(0),
                "replacement": f"except {exception_class} as {exception_var}:",
            }
        )

    return issues, content


def fix_exception_syntax(file_path, apply=False):
    """Fix Python 2 style exception handling in a file"""
    issues, content = find_exception_syntax_issues(file_path)

    if not issues:
        return False

    print(f"File: {file_path}")
    for issue in issues:
        print(f"  Line {issue['line']}: {issue['text']} -> {issue['replacement']}")

    if apply:
        # Replace all occurrences of the pattern
        new_content = re.sub(
            r"except\s+([A-Za-z0-9_]+)\s*,\s*([A-Za-z0-9_]+)\s*:",
            r"except \1 as \2:",
            content,
        )

        # Write back to the file
        with open(file_path, "w", encoding="utf-8") as f:
            f.write(new_content)

        print(f"  Fixed {len(issues)} exception syntax issues in {file_path}")

    return True


def main():
    parser = argparse.ArgumentParser(description="Fix Python 2 style exception syntax")
    parser.add_argument(
        "--path", required=True, help="Path to scan, can be a file or directory"
    )
    parser.add_argument("--apply", action="store_true", help="Apply fixes")
    args = parser.parse_args()

    path = args.path
    apply = args.apply

    if os.path.isfile(path) and path.endswith(".py"):
        fix_exception_syntax(path, apply)
    elif os.path.isdir(path):
        for root, _, files in os.walk(path):
            for file in files:
                if file.endswith(".py"):
                    file_path = os.path.join(root, file)
                    fix_exception_syntax(file_path, apply)
    else:
        print(f"Error: {path} is not a Python file or directory")
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
