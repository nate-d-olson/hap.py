#!/usr/bin/env python3
"""Simple test runner to check unit and integration test status."""

import subprocess
import sys
import time


def run_test_suite(test_path, name):
    """Run a test suite and report results."""
    print(f"\n{'='*50}")
    print(f"Running {name}")
    print(f"{'='*50}")

    try:
        start_time = time.time()
        result = subprocess.run(
            ["pytest", test_path, "-q", "--tb=short"],
            capture_output=True,
            text=True,
            timeout=300,
        )
        end_time = time.time()

        print(f"Execution time: {end_time - start_time:.2f} seconds")
        print(f"Return code: {result.returncode}")

        if result.stdout:
            print("\nSTDOUT:")
            print(result.stdout)

        if result.stderr:
            print("\nSTDERR:")
            print(result.stderr)

        return result.returncode == 0

    except subprocess.TimeoutExpired:
        print(f"ERROR: {name} timed out after 300 seconds")
        return False
    except Exception as e:
        print(f"ERROR running {name}: {e}")
        return False


def main():
    """Main test runner."""
    print("Testing specific fixes...")

    # Test the two specific fixes we made
    specific_tests = [
        "tests/unit/test_ga4gh_compliance.py::TestGA4GHMetrics::test_calculate_f1",
        "tests/unit/test_unit_quantify.py::TestQuantifyEngine::test_variant_classification",
    ]

    for test in specific_tests:
        print(f"\nTesting: {test}")
        result = subprocess.run(
            ["pytest", test, "-v"], capture_output=True, text=True, timeout=60
        )
        if result.returncode == 0:
            print("✅ PASSED")
        else:
            print("❌ FAILED")
            print(result.stdout)
            print(result.stderr)

    # Run unit tests
    unit_success = run_test_suite("tests/unit/", "Unit Tests")

    # Run integration tests
    integration_success = run_test_suite("tests/integration/", "Integration Tests")

    print(f"\n{'='*50}")
    print("SUMMARY")
    print(f"{'='*50}")
    print(f"Unit Tests: {'✅ PASSED' if unit_success else '❌ FAILED'}")
    print(f"Integration Tests: {'✅ PASSED' if integration_success else '❌ FAILED'}")

    return 0 if (unit_success and integration_success) else 1


if __name__ == "__main__":
    sys.exit(main())
