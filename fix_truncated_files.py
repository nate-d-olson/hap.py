#!/usr/bin/env python3
"""
Fix truncated Python files by using a combination of Python 3-specific versions and
original backups, ensuring Python 3 compatibility is maintained.
"""

import logging
import os
import re
import shutil
import sys

# Configure logging
logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

# List of truncated files identified
TRUNCATED_FILES = [
    "src/python/Haplo/happyroc.py",
    "src/python/Somatic/Pisces.py",
    "src/python/Somatic/Strelka.py",
    "src/python/Somatic/Varscan2.py",
    "src/python/Tools/bcftools.py",
    "src/python/Tools/roc.py",
]


def fix_truncated_file(file_path):
    """Fix a truncated file by using Python 3 version if available, or
    adapting the backup with Python 3 compatibility fixes."""

    backup_path = file_path + ".py2.bak"
    py3_path = None

    # Check if a Python 3 version exists
    if os.path.exists(file_path + ".py3"):
        py3_path = file_path + ".py3"
    else:
        # Handle different naming conventions for Python 3 files
        base_name = os.path.basename(file_path)
        dir_name = os.path.dirname(file_path)
        base_without_ext = os.path.splitext(base_name)[0]
        py3_alt_path = os.path.join(dir_name, f"{base_without_ext}_py3.py")

        if os.path.exists(py3_alt_path):
            py3_path = py3_alt_path

    if not os.path.exists(backup_path):
        logging.error(f"Backup file not found: {backup_path}")
        return False

    # Create a backup of the current file
    current_backup = file_path + ".truncated.bak"
    shutil.copy2(file_path, current_backup)
    logging.info(f"Created backup of current truncated file: {current_backup}")

    # Option 1: Use Python 3 specific version if it exists and is not truncated
    if py3_path and os.path.exists(py3_path):
        py3_lines = count_lines(py3_path)
        backup_lines = count_lines(backup_path)

        # If Python 3 version has at least 90% of the backup's content, use it
        if py3_lines >= 0.9 * backup_lines:
            logging.info(f"Using Python 3 specific version: {py3_path}")
            shutil.copy2(py3_path, file_path)
            return True

    # Option 2: Adapt the backup file with Python 3 compatibility fixes
    logging.info(
        f"Adapting backup file with Python 3 compatibility fixes: {backup_path}"
    )
    try:
        # Read backup file
        with open(backup_path, encoding="utf-8") as f:
            content = f.read()

        # Apply Python 3 fixes
        content = fix_python3_compatibility(content)

        # Write to original location
        with open(file_path, "w", encoding="utf-8") as f:
            f.write(content)

        return True

    except Exception as e:
        logging.error(f"Error fixing {file_path}: {e}")
        return False


def fix_python3_compatibility(content):
    """Apply Python 3 compatibility fixes to the content."""

    # Fix shebang line
    content = re.sub(
        r"#!/usr/bin/env python\s*$",
        "#!/usr/bin/env python3",
        content,
        flags=re.MULTILINE,
    )

    # Fix 'except X, e:' to 'except X as e:'
    content = re.sub(
        r"except\s+([A-Za-z0-9_]+),\s+([A-Za-z0-9_]+)(\s*:)",
        r"except \1 as \2\3",
        content,
    )

    # Fix file open operations for text files to include encoding
    content = re.sub(
        r'open\((.*?),\s*["\']r["\']\)', r'open(\1, "r", encoding="utf-8")', content
    )
    content = re.sub(
        r'open\((.*?),\s*["\']w["\']\)', r'open(\1, "w", encoding="utf-8")', content
    )

    # Fix print statements (assuming this has already been done by 2to3 or other tools)

    return content


def count_lines(file_path):
    """Count the number of lines in a file."""
    try:
        with open(file_path, encoding="utf-8") as f:
            return sum(1 for _ in f)
    except UnicodeDecodeError:
        # Try again with latin-1 encoding if utf-8 fails
        try:
            with open(file_path, encoding="latin-1") as f:
                return sum(1 for _ in f)
        except Exception as e:
            logging.error(f"Error reading {file_path} with latin-1 encoding: {e}")
            return -1
    except Exception as e:
        logging.error(f"Error reading {file_path}: {e}")
        return -1


def main():
    """Main function to fix truncated files."""

    fixed_count = 0
    failed_count = 0

    print("Starting to fix truncated files...")

    for file_path in TRUNCATED_FILES:
        print(f"Fixing truncated file: {file_path}")
        logging.info(f"Fixing truncated file: {file_path}")

        if fix_truncated_file(file_path):
            fixed_count += 1
            original_lines = count_lines(file_path)
            backup_lines = count_lines(file_path + ".py2.bak")
            print(
                f"  Fixed: {file_path} (now has {original_lines} lines, backup has {backup_lines} lines)"
            )
            logging.info(
                f"  Fixed: {file_path} (now has {original_lines} lines, backup has {backup_lines} lines)"
            )
        else:
            failed_count += 1
            print(f"  Failed to fix: {file_path}")
            logging.error(f"  Failed to fix: {file_path}")

    print(f"\nFixed {fixed_count} out of {len(TRUNCATED_FILES)} truncated files.")

    if failed_count > 0:
        print(f"Failed to fix {failed_count} files. Check the logs for details.")
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
