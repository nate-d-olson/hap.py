#!/usr/bin/env python3
"""
Fix common Python 2 artifacts in the codebase.

This script identifies and fixes common Python 2 artifacts like:
- print statements -> print functions
- xrange -> range
- Unicode literals imports
- old-style exceptions
- dict.iteritems() -> dict.items()
- StringIO imports
- old-style super() calls
"""

import argparse
import re
from pathlib import Path
from typing import Dict, Tuple


class Python2Fixer:
    """Fix Python 2 artifacts in source files."""

    def __init__(self, project_root: Path):
        """
        Initialize the Python 2 fixer.
        
        Args:
            project_root: Root directory of the project
        """
        self.project_root = project_root
        self.python_files = list(project_root.glob("src/python/**/*.py"))
        self.fixed_files = 0
        self.skipped_files = 0
        self.fixes_applied = {
            "print_statements": 0,
            "xrange": 0,
            "unicode_literals": 0,
            "old_style_except": 0,
            "old_style_super": 0,
            "iteritems": 0,
            "stringio": 0,
            "other": 0,
        }

    def fix_all_files(self, dry_run: bool = False) -> Dict[str, int]:
        """
        Fix Python 2 artifacts in all Python files.
        
        Args:
            dry_run: If True, don't actually make changes, just report
            
        Returns:
            Dictionary with counts of fixes applied
        """
        print(f"🔍 Found {len(self.python_files)} Python files to check")
        
        for file_path in self.python_files:
            self._fix_file(file_path, dry_run)
        
        print(f"\n✅ Fixed {self.fixed_files} files, skipped {self.skipped_files} files")
        print("Applied fixes:")
        for fix_type, count in self.fixes_applied.items():
            if count > 0:
                print(f"  - {fix_type}: {count}")
        
        return self.fixes_applied

    def _fix_file(self, file_path: Path, dry_run: bool = False) -> bool:
        """
        Fix Python 2 artifacts in a single file.
        
        Args:
            file_path: Path to the file to fix
            dry_run: If True, don't actually make changes, just report
            
        Returns:
            True if fixes were applied, False otherwise
        """
        try:
            # Read file content
            with open(file_path, "r", encoding="utf-8") as f:
                content = f.read()
            
            # Apply fixes
            fixed_content, fixes = self._apply_fixes(content)
            
            # If changes were made, write back
            if content != fixed_content:
                fix_list = ", ".join(f"{k} ({v})" for k, v in fixes.items() if v > 0)
                if dry_run:
                    print(f"🔍 Would fix {file_path.relative_to(self.project_root)}: {fix_list}")
                else:
                    with open(file_path, "w", encoding="utf-8") as f:
                        f.write(fixed_content)
                    print(f"✅ Fixed {file_path.relative_to(self.project_root)}: {fix_list}")
                
                self.fixed_files += 1
                for fix_type, count in fixes.items():
                    self.fixes_applied[fix_type] += count
                
                return True
            else:
                self.skipped_files += 1
                return False
                
        except Exception as e:
            print(f"⚠️ Error fixing {file_path}: {e}")
            self.skipped_files += 1
            return False

    def _apply_fixes(self, content: str) -> Tuple[str, Dict[str, int]]:
        """
        Apply fixes to file content.
        
        Args:
            content: File content to fix
            
        Returns:
            Tuple of (fixed content, dictionary of fix counts)
        """
        fixes = {
            "print_statements": 0,
            "xrange": 0,
            "unicode_literals": 0,
            "old_style_except": 0,
            "old_style_super": 0,
            "iteritems": 0,
            "stringio": 0,
            "other": 0,
        }
        
        # Fix print statements
        if "print " in content and not re.search(r"from __future__ import print_function", content):
            # Don't touch multiline prints or complex expressions for safety
            lines = content.split("\n")
            for i, line in enumerate(lines):
                # Skip comments
                if re.match(r'^\s*#', line):
                    continue
                    
                # Match standalone print statements, not already function calls
                if re.search(r'^\s*print\s+[^(]', line):
                    # Simple case: print followed by a string or variable
                    if '>>' not in line and '\\' not in line:
                        # Replace with print function
                        lines[i] = re.sub(r'^\s*print\s+(.*?)$', r'print(\1)', line)
                        fixes["print_statements"] += 1
            
            content = "\n".join(lines)
        
        # Fix xrange -> range
        if "xrange" in content:
            content = re.sub(r'\bxrange\b', 'range', content)
            fixes["xrange"] += len(re.findall(r'\bxrange\b', content))
        
        # Remove unicode_literals import
        if "unicode_literals" in content:
            content = re.sub(
                r'from __future__ import unicode_literals\s*\n',
                '',
                content
            )
            fixes["unicode_literals"] += 1
        
        # Fix old-style except statements
        if re.search(r'except\s+\w+\s*,\s*\w+\s*:', content):
            content = re.sub(
                r'except\s+(\w+)\s*,\s*(\w+)\s*:',
                r'except \1 as \2:',
                content
            )
            fixes["old_style_except"] += len(re.findall(r'except\s+\w+\s*,\s*\w+\s*:', content))
        
        # Fix old-style super calls
        if "super(" in content:
            # This is a complex transformation, so only handle simple cases
            content = re.sub(
                r'super\((\w+),\s*self\)\.(\w+)',
                r'super().\2',
                content
            )
            fixes["old_style_super"] += len(re.findall(r'super\((\w+),\s*self\)\.(\w+)', content))
        
        # Fix dict.iteritems()
        if ".iteritems" in content:
            content = re.sub(r'\.iteritems\(\)', '.items()', content)
            fixes["iteritems"] += len(re.findall(r'\.iteritems\(\)', content))
        
        # Fix StringIO imports
        if "StringIO" in content:
            if "from StringIO import StringIO" in content:
                content = content.replace(
                    "from StringIO import StringIO",
                    "from io import StringIO"
                )
                fixes["stringio"] += 1
            elif "import StringIO" in content:
                content = content.replace(
                    "import StringIO",
                    "import io"
                )
                # Also fix StringIO.StringIO() calls
                content = re.sub(
                    r'StringIO\.StringIO\(',
                    'io.StringIO(',
                    content
                )
                fixes["stringio"] += 1
        
        return content, fixes


def main():
    """Main function."""
    parser = argparse.ArgumentParser(
        description="Fix Python 2 artifacts in the codebase"
    )
    parser.add_argument(
        "--dry-run", action="store_true",
        help="Don't actually make changes, just report what would be done"
    )
    parser.add_argument(
        "--path", type=str, default=None,
        help="Path to the directory or file to process (default: src/python)"
    )
    
    args = parser.parse_args()
    
    # Get project root directory (parent of the script)
    project_root = Path(__file__).parent.parent
    
    # Create and run the fixer
    fixer = Python2Fixer(project_root)
    fixes = fixer.fix_all_files(args.dry_run)
    
    # Print summary
    total_fixes = sum(fixes.values())
    if total_fixes > 0:
        print(f"\n🎉 Applied {total_fixes} fixes to {fixer.fixed_files} files")
        if args.dry_run:
            print("Run without --dry-run to apply these changes")
    else:
        print("\n👍 No Python 2 artifacts found, codebase is clean!")


if __name__ == "__main__":
    main()
