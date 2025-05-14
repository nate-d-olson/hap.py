#!/usr/bin/env python3
"""
Automated Python 3 updates for hap.py codebase

This script automates several common Python 2 to 3 migration tasks:
1. Converting print statements to functions
2. Updating exception syntax from "except X, e" to "except X as e"
3. Replacing xrange with range
4. Updating dict methods (iteritems to items, etc.)
5. Adding proper string encoding/decoding for Python 3

Usage:
    python3 automate_py3_updates.py [--paths PATHS] [--apply]
"""

import argparse
import ast
import os
import re
import sys
from pathlib import Path
from typing import Any, Dict, List, Set, Tuple, Optional


class Python2To3Updater:
    """Helper class for updating Python 2 code to Python 3"""

    def __init__(self, verbose=False, apply=False):
        self.verbose = verbose
        self.apply = apply  # Whether to apply changes or just report them
        self.changes_made = {
            "print_statement": [],
            "except_syntax": [],
            "xrange": [],
            "dict_methods": [],
            "encoding": [],
        }

        # Common dict method replacements
        self.dict_methods_map = {
            "iteritems": "items",
            "iterkeys": "keys",
            "itervalues": "values",
            "viewitems": "items",
            "viewkeys": "keys",
            "viewvalues": "values",
        }

    def update_file(self, filepath: str) -> bool:
        """Update a single file for Python 3 compatibility
        
        Args:
            filepath: Path to the file to update
            
        Returns:
            True if changes were made, False otherwise
        """
        if self.verbose:
            print(f"Processing {filepath}...")
        
        try:
            with open(filepath, "r", encoding="utf-8") as f:
                content = f.read()
                
            # Make a copy of the original content
            original_content = content
            
            # Apply fixes
            content = self._fix_print_statements(content, filepath)
            content = self._fix_except_syntax(content, filepath)
            content = self._fix_xrange(content, filepath)
            content = self._fix_dict_methods(content, filepath)
            content = self._fix_encoding(content, filepath)
            
            # Check if any changes were made
            changes_made = original_content != content
            
            # Write the updated file if changes were made and apply is True
            if changes_made and self.apply:
                # Create backup
                backup_file = f"{filepath}.py2.bak"
                if not os.path.exists(backup_file):
                    with open(backup_file, "w", encoding="utf-8") as f:
                        f.write(original_content)
                    print(f"Created backup: {backup_file}")
                
                # Write updated content
                with open(filepath, "w", encoding="utf-8") as f:
                    f.write(content)
                print(f"Updated {filepath}")
                
            return changes_made
        except UnicodeDecodeError:
            print(f"Error: Cannot decode {filepath} as UTF-8. Skipping.")
            return False
            
    def _fix_print_statements(self, content: str, filepath: str) -> str:
        """Convert Python 2 print statements to Python 3 print functions"""
        # Simple pattern to match most print statements
        # This won't handle all cases perfectly, but it's a good starting point
        pattern = r'(?<!["\'])print\s+([^(].+?)(?=$|\n|\s*#)'
        
        def replace_print(match):
            # Extract the content of the print statement
            print_content = match.group(1).strip()
            # Record the change
            self.changes_made["print_statement"].append(filepath)
            # Return the print function
            return f"print({print_content})"
        
        return re.sub(pattern, replace_print, content)
        
    def _fix_except_syntax(self, content: str, filepath: str) -> str:
        """Update exception syntax from 'except X, e' to 'except X as e'"""
        # Pattern to match old except syntax
        pattern = r'except\s+([A-Za-z0-9_]+)\s*,\s*([A-Za-z0-9_]+):'
        
        def replace_except(match):
            # Extract the exception class and variable
            exception_class = match.group(1)
            exception_var = match.group(2)
            # Record the change
            self.changes_made["except_syntax"].append(filepath)
            # Return the updated syntax
            return f"except {exception_class} as {exception_var}:"
        
        return re.sub(pattern, replace_except, content)
        
    def _fix_xrange(self, content: str, filepath: str) -> str:
        """Replace xrange with range"""
        if "xrange" in content:
            # Make sure it's not in a string or comment
            pattern = r'(?<!["\'])xrange\('
            if re.search(pattern, content):
                self.changes_made["xrange"].append(filepath)
                # Replace xrange with range, being careful about boundaries
                return re.sub(r'(?<!["\'])xrange\(', "range(", content)
        return content
        
    def _fix_dict_methods(self, content: str, filepath: str) -> str:
        """Update dict methods (iteritems to items, etc.)"""
        updated_content = content
        for old_method, new_method in self.dict_methods_map.items():
            # Look for the old method
            pattern = rf'\.{old_method}\(\)'
            if re.search(pattern, updated_content):
                self.changes_made["dict_methods"].append(filepath)
                # Replace with the new method
                updated_content = re.sub(pattern, f".{new_method}()", updated_content)
        return updated_content
                
    def _fix_encoding(self, content: str, filepath: str) -> str:
        """Add proper string encoding/decoding for Python 3"""
        # This is a more complex task that likely requires manual attention
        # Here we just add common imports if they seem necessary
        if (".encode(" in content or ".decode(" in content) and "from __future__ import unicode_literals" not in content:
            # Add essential imports at the top of the file
            self.changes_made["encoding"].append(filepath)
            
            # Find the best place to add imports - after docstring and existing imports
            import_pattern = r'^(import\s+.*|\s*from\s+.*import.*)$'
            
            # Find the last import line
            last_import = None
            for match in re.finditer(import_pattern, content, re.MULTILINE):
                last_import = match.end()
            
            if last_import:
                # Add our imports after the last existing import
                parts = [
                    content[:last_import],
                    "\n\n# Added for Python 3 compatibility",
                    "# For better string encoding/decoding handling",
                    "from __future__ import unicode_literals",
                    content[last_import:]
                ]
                return "\n".join(parts)
                
        return content
                
    def summarize_changes(self):
        """Print a summary of the changes made or would be made"""
        total_changes = sum(len(changes) for changes in self.changes_made.values())
        
        print(f"\n=== Python 3 Update Summary ===")
        print(f"Total changes {'made' if self.apply else 'identified'}: {total_changes}")
        
        print(f"\nPrint statements {'converted' if self.apply else 'to convert'}: {len(self.changes_made['print_statement'])}")
        print(f"Exception syntax {'updated' if self.apply else 'to update'}: {len(self.changes_made['except_syntax'])}")
        print(f"xrange {'replaced' if self.apply else 'to replace'}: {len(self.changes_made['xrange'])}")
        print(f"Dict methods {'updated' if self.apply else 'to update'}: {len(self.changes_made['dict_methods'])}")
        print(f"String encoding {'fixed' if self.apply else 'to fix'}: {len(self.changes_made['encoding'])}")
        
        if not self.apply and total_changes > 0:
            print("\nTo apply these changes, run again with --apply")


def main():
    """Main function"""
    print("Starting Python 3 migration analysis...")
    parser = argparse.ArgumentParser(
        description="Update Python files for Python 3 compatibility"
    )
    parser.add_argument(
        "--paths", default="src/python", help="Comma-separated list of paths to scan"
    )
    parser.add_argument(
        "--apply", action="store_true", help="Apply the changes (default: just report)"
    )
    parser.add_argument(
        "--verbose", "-v", action="store_true", help="Show detailed information"
    )
    args = parser.parse_args()

    updater = Python2To3Updater(verbose=args.verbose, apply=args.apply)
    paths = [p.strip() for p in args.paths.split(",")]

    print(f"Scanning paths: {', '.join(paths)}")
    python_files = []

    # Find all Python files in the specified paths
    for path in paths:
        for root, _, files in os.walk(path):
            for file in files:
                if file.endswith(".py") and not file.startswith("mock_"):
                    python_files.append(os.path.join(root, file))

    print(f"Found {len(python_files)} Python files to process")

    # Process each file
    files_updated = 0
    for i, file in enumerate(python_files):
        if i % 10 == 0 and args.verbose:
            print(f"Processing file {i+1}/{len(python_files)}: {file}")
        if updater.update_file(file):
            files_updated += 1

    print(f"\nProcessed {len(python_files)} files")
    print(f"Files {'updated' if args.apply else 'with issues'}: {files_updated}")
    
    updater.summarize_changes()
    
    if args.apply:
        print("\nNote: Original files have been backed up with .py2.bak extension")
        print("Check the updated files to ensure correctness")


if __name__ == "__main__":
    main()
