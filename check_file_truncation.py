#!/usr/bin/env python3
"""
Check if any Python files were truncated during migration by comparing
line counts with their backup (.py2.bak) versions.
"""

import os
import sys
import glob
import logging

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')

def count_lines(file_path):
    """Count the number of lines in a file."""
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            return sum(1 for _ in f)
    except UnicodeDecodeError:
        # Try again with latin-1 encoding if utf-8 fails
        try:
            with open(file_path, 'r', encoding='latin-1') as f:
                return sum(1 for _ in f)
        except Exception as e:
            logging.error(f"Error reading {file_path} with latin-1 encoding: {e}")
            return -1
    except Exception as e:
        logging.error(f"Error reading {file_path}: {e}")
        return -1

def main():
    """Main function to check for truncated files."""
    # Find all .py2.bak files
    backup_files = []
    for root, _, _ in os.walk('src/python'):
        backup_files.extend(glob.glob(os.path.join(root, '*.py2.bak')))
    
    # For each backup file, compare line count with the original
    truncated_files = []
    total_checked = 0
    
    for backup_file in sorted(backup_files):
        # Get original filename (remove .py2.bak)
        original_file = backup_file[:-8] 
        
        if not os.path.exists(original_file):
            logging.warning(f"Original file not found: {original_file}")
            continue
        
        backup_lines = count_lines(backup_file)
        original_lines = count_lines(original_file)
        
        total_checked += 1
        
        # If the migrated file has significantly fewer lines (>5% less),
        # consider it possibly truncated
        if original_lines > 0 and backup_lines > 0:
            # Calculate the ratio
            ratio = original_lines / backup_lines
            
            if ratio < 0.95:  # More than 5% fewer lines in migrated file
                truncated_files.append({
                    'original': original_file,
                    'backup': backup_file,
                    'original_lines': original_lines,
                    'backup_lines': backup_lines,
                    'ratio': ratio
                })
    
    # Print results
    print(f"\nChecked {total_checked} files for truncation.\n")
    
    if not truncated_files:
        print("No truncated files detected!")
    else:
        print(f"Found {len(truncated_files)} potentially truncated files:")
        print("=" * 80)
        for info in truncated_files:
            print(f"File: {info['original']}")
            print(f"  Original: {info['original_lines']} lines")
            print(f"  Backup:   {info['backup_lines']} lines")
            print(f"  Ratio:    {info['ratio']:.2f}")
            print("-" * 80)

if __name__ == "__main__":
    main()
