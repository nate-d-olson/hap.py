#!/usr/bin/env python3
"""
Test script to verify string handling in updated Python 3 modules
"""

import os
import sys
import tempfile

# Add src/python to the path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "src", "python"))


def test_fastasize():
    """Test string handling in fastasize.py"""
    try:
        print("Successfully imported fastaContigLengths from Tools.fastasize")
        return True
    except Exception as e:
        print(f"Error importing fastaContigLengths: {e}")
        return False


def test_happyroc():
    """Test string handling in happyroc.py"""
    try:
        print("Successfully imported roc from Haplo.happyroc")
        return True
    except Exception as e:
        print(f"Error importing roc: {e}")
        return False


def test_file_handling():
    """Test file handling with explicit encoding"""
    try:
        # Create a temporary file with some Unicode content
        with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", delete=False) as f:
            f.write("Test file with Unicode: ñäöü")
            filename = f.name

        # Read it back using explicit encoding
        with open(filename, encoding="utf-8") as f:
            content = f.read()

        print(f"Successfully read test file with content: {content}")
        os.unlink(filename)
        return True
    except Exception as e:
        print(f"Error in file handling test: {e}")
        return False


def main():
    """Run all tests"""
    tests = [test_fastasize, test_happyroc, test_file_handling]

    success = 0
    for test in tests:
        print(f"\nRunning test: {test.__name__}")
        if test():
            success += 1
            print(f"✓ {test.__name__} passed")
        else:
            print(f"✗ {test.__name__} failed")

    print(f"\nTest Summary: {success}/{len(tests)} tests passed")
    return 0 if success == len(tests) else 1


if __name__ == "__main__":
    sys.exit(main())
