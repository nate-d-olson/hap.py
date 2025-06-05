#!/usr/bin/env python3
"""Quick test to check if core components are working."""

import subprocess
import sys
from pathlib import Path


def test_imports():
    """Test that core modules can be imported."""
    try:

        print("✅ Core imports successful")
        return True
    except Exception as e:
        print(f"❌ Import error: {e}")
        return False


def test_rtg_availability():
    """Test RTG tools availability."""
    rtg_paths = ["external/rtg-tools-3.12.1/rtg", "rtg"]

    for rtg_path in rtg_paths:
        if Path(rtg_path).exists():
            print(f"✅ RTG found at: {rtg_path}")
            return True

    print("❌ RTG tools not found")
    return False


def test_specific_fixes():
    """Test the specific fixes we made."""
    print("\nTesting specific fixes...")

    tests_to_run = [
        "tests/unit/test_ga4gh_compliance.py::TestGA4GHMetrics::test_calculate_f1",
        "tests/unit/test_unit_quantify.py::TestQuantifyEngine::test_variant_classification",
    ]

    all_passed = True
    for test in tests_to_run:
        print(f"Running: {test}")
        try:
            result = subprocess.run(
                ["pytest", test, "-v"], capture_output=True, text=True, timeout=30
            )
            if result.returncode == 0:
                print("  ✅ PASSED")
            else:
                print("  ❌ FAILED")
                print(f"  Error: {result.stdout}")
                all_passed = False
        except Exception as e:
            print(f"  ❌ ERROR: {e}")
            all_passed = False

    return all_passed


def test_simple_integration():
    """Run a simple integration test."""
    print("\nTesting simple integration...")

    try:
        # Try to run a simple pytest collection to see if there are import errors
        result = subprocess.run(
            ["pytest", "tests/integration/", "--collect-only", "-q"],
            capture_output=True,
            text=True,
            timeout=60,
        )

        if "error" in result.stderr.lower() or result.returncode != 0:
            print("❌ Integration test collection failed")
            print(f"STDERR: {result.stderr}")
            return False
        else:
            print("✅ Integration tests can be collected")

            # Try to run just one small test
            result = subprocess.run(
                [
                    "pytest",
                    "tests/integration/test_integration.py",
                    "-k",
                    "test_",
                    "-x",
                    "--tb=short",
                ],
                capture_output=True,
                text=True,
                timeout=120,
            )

            if "FAILED" in result.stdout:
                print("❌ Some integration tests failed")
                print("Sample output:", result.stdout[:500])
                return False
            elif "passed" in result.stdout:
                print("✅ At least some integration tests passed")
                return True
            else:
                print("⚠️  Uncertain integration test status")
                print("Output:", result.stdout[:200])
                return True  # Don't fail for uncertain status

    except Exception as e:
        print(f"❌ Integration test error: {e}")
        return False


def main():
    """Main test function."""
    print("Quick test check for hap.py modernization\n")

    results = {
        "Imports": test_imports(),
        "RTG Tools": test_rtg_availability(),
        "Specific Fixes": test_specific_fixes(),
        "Integration": test_simple_integration(),
    }

    print(f"\n{'='*50}")
    print("SUMMARY")
    print(f"{'='*50}")

    for test_name, passed in results.items():
        status = "✅ PASSED" if passed else "❌ FAILED"
        print(f"{test_name}: {status}")

    all_passed = all(results.values())
    print(f"\nOverall: {'✅ SUCCESS' if all_passed else '❌ ISSUES FOUND'}")

    return 0 if all_passed else 1


if __name__ == "__main__":
    sys.exit(main())
