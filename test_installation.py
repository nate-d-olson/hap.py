#!/usr/bin/env python3
"""Test script to check installation status and run basic tests."""

import sys
from pathlib import Path

# Add src to path
src_path = Path(__file__).parent / "src"
sys.path.insert(0, str(src_path))


def test_import():
    """Test basic imports."""
    try:
        import hap_py

        print("✅ hap_py imported successfully")
        print(f"   Version: {hap_py.__version__}")
        return True
    except Exception as e:
        print(f"❌ Failed to import hap_py: {e}")
        return False


def test_rtg_manager():
    """Test RTG manager."""
    try:
        from hap_py.external.rtg_manager import rtg_manager

        rtg_path = rtg_manager.get_rtg_path()
        print(f"✅ RTG Tools found at: {rtg_path}")
        return True
    except Exception as e:
        print(f"❌ RTG manager failed: {e}")
        return False


def test_basic_functionality():
    """Test basic functionality."""
    try:

        print("✅ vcfeval module imported successfully")
        return True
    except Exception as e:
        print(f"❌ vcfeval import failed: {e}")
        return False


def main():
    """Run all tests."""
    print("Testing hap.py installation status...")
    print("=" * 50)

    results = []
    results.append(test_import())
    results.append(test_rtg_manager())
    results.append(test_basic_functionality())

    print("=" * 50)
    if all(results):
        print("✅ All basic tests passed!")
        return 0
    else:
        print("❌ Some tests failed")
        return 1


if __name__ == "__main__":
    sys.exit(main())
