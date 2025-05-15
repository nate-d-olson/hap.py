#!/usr/bin/env python3
"""
Verify Python 3 compatibility issues in hap.py codebase.

This script analyzes the Python modules to identify common Python 2 to 3 migration issues
that may need to be fixed during the migration process. It also tracks migration progress
and generates a progress report.

Usage:
    python3 check_py3_issues.py [--paths PATHS] [--generate-report]
"""

import argparse
import ast
import os
import re
import sys
from pathlib import Path
from typing import Dict, List, Set, Tuple


class Python2To3Issues:
    """Helper class to find Python 2 to 3 compatibility issues"""

    def __init__(self, verbose=False):
        self.verbose = verbose
        self.issues = {
            "print_statement": [],
            "unicode_issues": [],
            "xrange": [],
            "dict_methods": [],
            "division": [],
            "except_syntax": [],
            "imports": [],
            "cython_integration": [],
        }

        # Common dict method differences between Py2 and Py3
        self.dict_methods_py2 = {
            "iteritems",
            "iterkeys",
            "itervalues",
            "viewitems",
            "viewkeys",
            "viewvalues",
        }

        # Modules renamed or relocated in Python 3
        self.renamed_modules = {
            "ConfigParser": "configparser",
            "Queue": "queue",
            "SocketServer": "socketserver",
            "urllib2": "urllib.request",
            "urlparse": "urllib.parse",
            "StringIO": "io",
            "cStringIO": "io",
        }

    def check_file(self, filepath):
        """Check a single file for Python 2 to 3 issues"""
        try:
            with open(filepath, "r", encoding="utf-8") as f:
                content = f.read()

            # Parse the file with AST to analyze its structure
            try:
                tree = ast.parse(content)
                self._check_ast(tree, filepath)
            except SyntaxError as e:
                print(f"Syntax error in {filepath}: {e}")
                # Continue with regex-based checks even if AST parsing fails

            # Perform regex-based checks
            self._check_with_regex(content, filepath)

            # Special checks for cython files
            if filepath.endswith(".pyx") or filepath.endswith(".pxd"):
                self._check_cython_issues(content, filepath)

        except UnicodeDecodeError:
            print(f"Cannot decode file {filepath} as UTF-8. Skipping.")

    def _check_ast(self, tree, filepath):
        """Check Python 2 to 3 issues using AST"""
        for node in ast.walk(tree):
            # Note: We can't check for print statements using AST in Python 3
            # since ast.Print doesn't exist. We use regex for this instead.

            # Check for except syntax (Python 3 uses ast.Try instead of ast.TryExcept)
            if hasattr(ast, "TryExcept") and isinstance(node, ast.TryExcept):
                for handler in node.handlers:
                    if handler.name is not None and not isinstance(
                        handler.name, ast.Name
                    ):
                        self.issues["except_syntax"].append(
                            {
                                "file": filepath,
                                "line": handler.lineno,
                                "message": "Old exception syntax (except Exception, e)",
                            }
                        )
            elif isinstance(node, ast.Try):
                for handler in node.handlers:
                    # In Python 3, exception variable is in handler.name
                    if (
                        hasattr(handler, "name")
                        and handler.name is not None
                        and not isinstance(handler.name, ast.Name)
                    ):
                        self.issues["except_syntax"].append(
                            {
                                "file": filepath,
                                "line": handler.lineno,
                                "message": "Old exception syntax (except Exception, e)",
                            }
                        )

            # Check for xrange usage
            if (
                isinstance(node, ast.Call)
                and isinstance(node.func, ast.Name)
                and node.func.id == "xrange"
            ):
                self.issues["xrange"].append(
                    {
                        "file": filepath,
                        "line": node.lineno,
                        "message": "Replace xrange with range",
                    }
                )

            # Check for dict methods that have changed
            if isinstance(node, ast.Attribute):
                if (
                    isinstance(node.value, ast.Name)
                    and node.attr in self.dict_methods_py2
                ):
                    self.issues["dict_methods"].append(
                        {
                            "file": filepath,
                            "line": getattr(node, "lineno", "?"),
                            "message": f"Replace {node.attr} with Python 3 equivalent",
                        }
                    )

            # Check for division that might need to be updated
            if isinstance(node, ast.BinOp) and isinstance(node.op, ast.Div):
                if isinstance(node.left, ast.Num) and isinstance(node.right, ast.Num):
                    # Integer division in Python 2 that might need // in Python 3
                    if (
                        hasattr(node.left, "n")
                        and hasattr(node.right, "n")
                        and node.left.n == int(node.left.n)
                        and node.right.n == int(node.right.n)
                    ):
                        self.issues["division"].append(
                            {
                                "file": filepath,
                                "line": node.lineno,
                                "message": "Consider using // for integer division",
                            }
                        )

            # Check for imports of renamed modules
            if isinstance(node, ast.Import):
                for name in node.names:
                    if name.name in self.renamed_modules:
                        self.issues["imports"].append(
                            {
                                "file": filepath,
                                "line": node.lineno,
                                "message": f"Update {name.name} to {self.renamed_modules[name.name]}",
                            }
                        )
            elif isinstance(node, ast.ImportFrom):
                if node.module in self.renamed_modules:
                    self.issues["imports"].append(
                        {
                            "file": filepath,
                            "line": node.lineno,
                            "message": f"Update {node.module} to {self.renamed_modules[node.module]}",
                        }
                    )

    def _check_with_regex(self, content, filepath):
        """Check Python 2 to 3 issues using regex patterns"""
        # Check for Python 2 print statements (not wrapped in parentheses)
        print_pattern = r'(?<!["\'])print\s+[^(]'
        for match in re.finditer(print_pattern, content):
            # Make sure it's not inside a comment or string
            line_start = content.rfind("\n", 0, match.start())
            if line_start == -1:
                line_start = 0
            else:
                line_start += 1  # Skip the newline

            line = content[line_start : match.start()].strip()
            if (
                not line.startswith("#")
                and not line.startswith('"')
                and not line.startswith("'")
            ):
                lineno = content[: match.start()].count("\n") + 1
                self.issues["print_statement"].append(
                    {
                        "file": filepath,
                        "line": lineno,
                        "message": "Python 2 print statement",
                    }
                )

        # Check for unicode/str issues
        unicode_patterns = [
            r"unicode\(",
            r'u[\'"]',
            r"str\(.*\)",
            r"\.encode\(",
            r"\.decode\(",
        ]

        for pattern in unicode_patterns:
            for match in re.finditer(pattern, content):
                lineno = content[: match.start()].count("\n") + 1
                self.issues["unicode_issues"].append(
                    {
                        "file": filepath,
                        "line": lineno,
                        "message": "Potential string/unicode issue",
                    }
                )

    def _check_cython_issues(self, content, filepath):
        """Check for Cython-specific Python 3 issues"""
        # Check for language level directives
        if not re.search(r"#\s*cython:\s*language_level\s*=\s*3", content):
            self.issues["cython_integration"].append(
                {
                    "file": filepath,
                    "line": 0,
                    "message": "Missing language_level=3 directive",
                }
            )

        # Check for string handling in Cython
        if re.search(r"cdef\s+char\s*\*.*=\s*(?!<).*str", content):
            search_results = re.findall(
                r"\n.*cdef\s+char\s*\*.*=\s*(?!<).*str", content
            )
            if search_results:
                lineno = content[: content.find(search_results[0])].count("\n") + 1
                self.issues["cython_integration"].append(
                    {
                        "file": filepath,
                        "line": lineno,
                        "message": "Potential string handling issue with char* and Python 3 str",
                    }
                )

        # Check for Python C API calls that changed in Python 3
        py3_api_changes = [
            r"PyString_",
            r"PyInt_",
            r"Py_TPFLAGS_HAVE_WEAKREFS",
            r"PyClass_",
        ]

        for pattern in py3_api_changes:
            for match in re.finditer(pattern, content):
                lineno = content[: match.start()].count("\n") + 1
                self.issues["cython_integration"].append(
                    {
                        "file": filepath,
                        "line": lineno,
                        "message": f"Python C API change needed: {pattern}",
                    }
                )

    def report(self):
        """Generate a report of all issues found"""
        total_issues = sum(len(issues) for issues in list(self.issues.values()))

        print("\n=== Python 2 to 3 Migration Issues Report ===")
        print(f"Total issues found: {total_issues}\n")

        # Generic function to print issues with a consistent format
        def print_issue_category(category, description):
            issues = self.issues[category]
            if not issues:
                return

            issue_count = len(issues)
            print(f"{description} ({issue_count}):")
            if self.verbose:
                for issue in issues:
                    if isinstance(issue, tuple):
                        # Handle old tuple format (backward compatibility)
                        if len(issue) == 2:
                            file, line = issue
                            message = ""
                        elif len(issue) == 3:
                            file, line, message = issue
                        else:
                            continue
                        print(f"  {file}:{line}{' - ' + message if message else ''}")
                    elif isinstance(issue, dict):
                        # Handle new dict format
                        file = issue.get("file", "unknown")
                        line = issue.get("line", "?")
                        message = issue.get("message", "")
                        print(f"  {file}:{line}{' - ' + message if message else ''}")
            print()

        # Print each category of issues
        print_issue_category(
            "print_statement", "Print Statements: Replace with print() function"
        )
        print_issue_category(
            "unicode_issues",
            "Unicode/String Issues: Check string handling, unicode() vs str()",
        )
        print_issue_category("xrange", "xrange Usage: Replace with range()")
        print_issue_category(
            "dict_methods",
            "Dict Methods: Update iterator methods (iteritems -> items, etc.)",
        )
        print_issue_category(
            "division", "Division Operator: Check integer division behavior (/ vs //)"
        )
        print_issue_category(
            "except_syntax", "Exception Syntax: Update to 'except X as y'"
        )
        print_issue_category("imports", "Renamed Imports: Update module imports")
        print_issue_category(
            "cython_integration",
            "Cython Integration Issues: Check Cython/Python 3 compatibility",
        )

        return total_issues


def generate_migration_report(checker, python_files):
    """Generate a migration progress report"""
    total_files = len(python_files)
    migrated_files = 0
    partially_migrated = 0
    issues_by_module = {}

    # Group issues by module/package
    for issue_type, issues in list(checker.issues.items()):
        for issue in issues:
            # Handle both dict and tuple formats
            if isinstance(issue, dict):
                path = issue.get("file", "")
            elif isinstance(issue, tuple) and len(issue) >= 1:
                path = issue[0]
            else:
                continue

            if path:
                parts = path.split("/")
                module = parts[2] if len(parts) > 2 else path
                if module not in issues_by_module:
                    issues_by_module[module] = 0
                issues_by_module[module] += 1

    # Calculate files with no issues
    for file in python_files:
        file_has_issues = False
        for issue_type, issues in list(checker.issues.items()):
            for issue in issues:
                # Handle both dict and tuple formats
                if isinstance(issue, dict) and issue.get("file") == file:
                    file_has_issues = True
                    break
                elif isinstance(issue, tuple) and len(issue) >= 1 and issue[0] == file:
                    file_has_issues = True
                    break
            if file_has_issues:
                break
        if not file_has_issues:
            migrated_files += 1
        elif file_has_issues:
            partially_migrated += 1

    # Generate report
    print("\n=== Python 3 Migration Progress Report ===")
    print(f"Total Python files: {total_files}")
    print(
        f"Fully migrated files: {migrated_files} ({migrated_files * 100 / total_files:.1f}%)"
    )
    print(
        f"Partially migrated files: {partially_migrated} ({partially_migrated * 100 / total_files:.1f}%)"
    )
    print(
        f"Files with issues: {total_files - migrated_files} ({(total_files - migrated_files) * 100 / total_files:.1f}%)"
    )

    print("\nIssues by module:")
    for module, count in sorted(
        list(issues_by_module.items()), key=lambda x: x[1], reverse=True
    ):
        print(f"  {module}: {count} issues")

    return migrated_files / total_files  # Return migration progress percentage


def main():
    """Main function"""
    parser = argparse.ArgumentParser(
        description="Check Python 2 to 3 compatibility issues in code"
    )
    parser.add_argument(
        "--paths", default="src/python", help="Comma-separated list of paths to scan"
    )
    parser.add_argument(
        "--verbose", "-v", action="store_true", help="Show detailed issue locations"
    )
    parser.add_argument(
        "--generate-report",
        action="store_true",
        help="Generate migration progress report",
    )
    args = parser.parse_args()

    checker = Python2To3Issues(verbose=args.verbose)
    paths = [p.strip() for p in args.paths.split(",")]

    print(f"Scanning paths: {', '.join(paths)}")
    python_files = []

    # Find all Python files in the specified paths
    for path in paths:
        for root, _, files in os.walk(path):
            for file in files:
                if (
                    file.endswith(".py")
                    or file.endswith(".pyx")
                    or file.endswith(".pxd")
                ):
                    python_files.append(os.path.join(root, file))

    print(f"Found {len(python_files)} Python files to scan")

    # Check each file
    for i, file in enumerate(python_files):
        if i % 10 == 0:
            print(f"Scanning file {i+1}/{len(python_files)}: {file}")
        checker.check_file(file)

    # Report results
    total_issues = checker.report()

    # Generate migration progress report if requested
    if args.generate_report:
        progress = generate_migration_report(checker, python_files)
        print(f"\nOverall migration progress: {progress:.1%}")

    if total_issues > 0:
        print(f"\nFound {total_issues} issues to fix for Python 3 compatibility.")
        return 1
    else:
        print("\nNo Python 2 to 3 compatibility issues found!")
        return 0


if __name__ == "__main__":
    sys.exit(main())
