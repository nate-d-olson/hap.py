#!/usr/bin/env python3
"""
Script to identify and fix common string/unicode issues in Python 2 to 3 migration
"""

import argparse
import os
import re
import sys


def find_string_unicode_issues(file_path):
    """Find potential string/unicode issues in a file"""
    try:
        with open(file_path, "r", encoding="utf-8") as f:
            content = f.read()
    except UnicodeDecodeError:
        print(f"Warning: Unable to read {file_path} as UTF-8. Skipping.")
        return [], None

    issues = []

    # Check for file open without encoding
    open_pattern = r'open\(\s*([^,]+)\s*,\s*[\'"]r[\'"]\s*\)'
    for match in re.finditer(open_pattern, content):
        lineno = content[: match.start()].count("\n") + 1
        issues.append(
            {
                "type": "file_open",
                "line": lineno,
                "text": match.group(0),
                "replacement": f'open({match.group(1)}, "r", encoding="utf-8")',
            }
        )

    # Check for file open with text mode but no encoding
    open_text_pattern = r'open\(\s*([^,]+)\s*,\s*[\'"]([rwat])[\'"](\s*\))'
    for match in re.finditer(open_text_pattern, content):
        lineno = content[: match.start()].count("\n") + 1
        issues.append(
            {
                "type": "file_open_text",
                "line": lineno,
                "text": match.group(0),
                "replacement": f'open({match.group(1)}, "{match.group(2)}", encoding="utf-8"{match.group(3)}',
            }
        )

    # Check for reading files without decoding
    read_pattern = r"\.read\(\)(\s*\.split|\s*\.strip|\s*\.replace|\s*\+)"
    for match in re.finditer(read_pattern, content):
        lineno = content[: match.start()].count("\n") + 1
        if (
            ".decode(" not in content[max(0, match.start() - 20) : match.start()]
            and "encoding=" not in content[max(0, match.start() - 50) : match.start()]
        ):
            issues.append(
                {
                    "type": "file_read",
                    "line": lineno,
                    "text": "Potential bytes output from read() without decoding",
                    "requires_manual_check": True,
                    "context": content[
                        max(0, match.start() - 50) : min(len(content), match.end() + 50)
                    ],
                }
            )

    # Check for string literal prefixes (less common issue)
    unicode_str_pattern = r'u[\'"](?:[^\'"\\]|\\.)*[\'"]'
    for match in re.finditer(unicode_str_pattern, content):
        lineno = content[: match.start()].count("\n") + 1
        original = match.group(0)
        replacement = original[1:]  # Remove the 'u' prefix
        issues.append(
            {
                "type": "unicode_literal",
                "line": lineno,
                "text": original,
                "replacement": replacement,
            }
        )

    # Check for subprocess calls that might need text mode
    subprocess_patterns = [
        (
            r"subprocess\.Popen\((.*?)stdout\s*=\s*subprocess\.PIPE(.*?)\)",
            r"subprocess.Popen(\1stdout=subprocess.PIPE\2, universal_newlines=True)",
        ),
        (
            r"subprocess\.check_output\((.*?)\)",
            r"subprocess.check_output(\1, universal_newlines=True)",
        ),
    ]

    for pattern, replacement in subprocess_patterns:
        for match in re.finditer(pattern, content):
            if "universal_newlines" not in match.group(0):
                lineno = content[: match.start()].count("\n") + 1
                issues.append(
                    {
                        "type": "subprocess_text_mode",
                        "line": lineno,
                        "text": match.group(0),
                        "replacement": replacement,
                    }
                )

    # Check for str() calls that might need encoding/decoding
    str_pattern = r"str\(([^)]+)\)"
    for match in re.finditer(str_pattern, content):
        lineno = content[: match.start()].count("\n") + 1
        context_before = content[max(0, match.start() - 50) : match.start()]
        context_after = content[match.end() : min(len(content), match.end() + 50)]
        var_name = match.group(1).strip()

        # Consider all str() calls as potentially problematic in Python 3
        # especially those in logging statements or string formatting
        if (
            "%" in context_after
            or "logging" in context_before
            or "bytes" in context_before
            or "encode" in context_before
            or "decode" in context_before
            or "read(" in context_before
            or "subprocess" in context_before
            or "Exception" in var_name
        ):

            # Flag this as a potential issue
            issues.append(
                {
                    "type": "str_call",
                    "line": lineno,
                    "text": match.group(0),  # The actual str() call
                    "requires_manual_check": True,
                    "context": content[
                        max(0, match.start() - 50) : min(len(content), match.end() + 50)
                    ],
                    "replacement_suggestion": f"{var_name}.decode('utf-8') if isinstance({var_name}, bytes) else str({var_name})",
                    "replacement": f"{var_name}.decode('utf-8') if isinstance({var_name}, bytes) else str({var_name})",
                }
            )

    # Check for concatenation of strings that might be bytes in Python 3
    concat_patterns = [
        r"(\w+)\.decode\([^)]*\)\s*\+",
        r"\+\s*(\w+)\.decode\([^)]*\)",
        r"(\w+)\.encode\([^)]*\)\s*\+",
        r"\+\s*(\w+)\.encode\([^)]*\)",
    ]

    for pattern in concat_patterns:
        for match in re.finditer(pattern, content):
            lineno = content[: match.start()].count("\n") + 1
            issues.append(
                {
                    "type": "string_concat",
                    "line": lineno,
                    "text": f"Potential string concatenation issue with encoded/decoded data: {match.group(0)}",
                    "requires_manual_check": True,
                    "context": content[
                        max(0, match.start() - 50) : min(len(content), match.end() + 50)
                    ],
                }
            )

    # Check for string formatting with potential bytes objects
    format_pattern = r"%\s*(\w+)"
    for match in re.finditer(format_pattern, content):
        var_name = match.group(1)
        # Check if variable is likely bytes
        if (
            re.search(r'[\'"]%s[\'"]' % re.escape(var_name), content)
            or f"{var_name}.encode(" in content
            or f"{var_name}.decode(" in content
        ):
            lineno = content[: match.start()].count("\n") + 1
            issues.append(
                {
                    "type": "string_format",
                    "line": lineno,
                    "text": f"Potential string formatting with bytes object: %{var_name}",
                    "requires_manual_check": True,
                    "context": content[
                        max(0, match.start() - 50) : min(len(content), match.end() + 50)
                    ],
                    "replacement_suggestion": f"{var_name}.decode('utf-8') if isinstance({var_name}, bytes) else {var_name}",
                }
            )

    return issues, content


def fix_string_unicode_issues(file_path, apply=False, verbose=False, autofix_all=False):
    """Identify and optionally fix string/unicode issues in a file"""
    issues, content = find_string_unicode_issues(file_path)

    if not issues or content is None:
        return False

    print(f"File: {file_path}")
    for issue in issues:
        if issue.get("requires_manual_check", False):
            if verbose and issue.get("context"):
                print(
                    f"  Line {issue['line']}: {issue['text']} (Requires manual check)"
                )
                print(f"  Context: {issue.get('context', '')}")
                if issue.get("replacement_suggestion"):
                    print(f"  Suggested fix: {issue.get('replacement_suggestion')}")
            else:
                print(
                    f"  Line {issue['line']}: {issue['text']} (Requires manual check)"
                )
        else:
            print(f"  Line {issue['line']}: {issue['text']} -> {issue['replacement']}")

    if apply:
        new_content = content

        # Apply simple replacements and complex replacements if autofix_all is True
        for issue in issues:
            if not issue.get("requires_manual_check", False) or (
                autofix_all and issue.get("replacement")
            ):
                if issue["type"] in [
                    "file_open",
                    "file_open_text",
                    "unicode_literal",
                    "subprocess_text_mode",
                    "str_call",
                ]:
                    # Try to use replacement pattern if available
                    if "replacement" in issue:
                        if verbose:
                            print(
                                f"  Applying fix for line {issue['line']}: {issue['text']} -> {issue['replacement']}"
                            )

                        # Find the text in the content by line number to be more precise
                        lines = new_content.split("\n")
                        if issue["line"] <= len(lines):
                            line = lines[issue["line"] - 1]
                            if issue["text"] in line:
                                # Replace just this instance in this specific line
                                new_line = line.replace(
                                    issue["text"], issue["replacement"]
                                )
                                lines[issue["line"] - 1] = new_line
                                new_content = "\n".join(lines)

        # Write back to the file if changes were made
        if new_content != content:
            with open(file_path, "w", encoding="utf-8") as f:
                f.write(new_content)
            print(f"  Fixed string/unicode issues in {file_path}")
        else:
            print(f"  No fixes applied to {file_path}")

    return True


def main():
    parser = argparse.ArgumentParser(
        description="Fix Python 2 to 3 string/unicode issues"
    )
    parser.add_argument(
        "--path", required=True, help="Path to scan, can be a file or directory"
    )
    parser.add_argument("--apply", action="store_true", help="Apply automatic fixes")
    parser.add_argument(
        "--verbose", "-v", action="store_true", help="Show verbose output"
    )
    parser.add_argument(
        "--autofix-all",
        action="store_true",
        help="Automatically fix all issues, including those that would normally require manual checks",
    )
    args = parser.parse_args()

    path = args.path
    apply = args.apply
    verbose = args.verbose
    autofix_all = args.autofix_all

    files_with_issues = 0
    issues_fixed = 0

    if os.path.isfile(path) and (path.endswith(".py") or path.endswith(".pyx")):
        if fix_string_unicode_issues(path, apply, verbose, autofix_all):
            files_with_issues += 1
    elif os.path.isdir(path):
        for root, _, files in os.walk(path):
            for file in files:
                if file.endswith(".py") or file.endswith(".pyx"):
                    file_path = os.path.join(root, file)
                    if fix_string_unicode_issues(
                        file_path, apply, verbose, autofix_all
                    ):
                        files_with_issues += 1
    else:
        print(f"Error: {path} is not a Python file or directory")
        return 1

    print(f"\nSummary: Found {files_with_issues} files with string/unicode issues")
    if apply:
        print(f"Applied automatic fixes where possible")
        if autofix_all:
            print(
                f"Applied fixes to all detected issues, including those that normally require manual checks"
            )
    else:
        print(f"Run with --apply to fix issues automatically")
        print(f"Run with --autofix-all to automatically fix all issues")

    return 0


if __name__ == "__main__":
    sys.exit(main())
