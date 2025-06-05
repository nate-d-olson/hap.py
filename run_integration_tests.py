#!/usr/bin/env python3
"""
Script to run integration tests and capture results systematically
"""

import subprocess
import sys
from pathlib import Path


def run_test(test_path):
    """Run a single test and return result"""
    try:
        cmd = [sys.executable, "-m", "pytest", test_path, "-v", "--tb=short"]

        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=300,  # 5 minute timeout per test
        )

        return {
            "test": test_path,
            "returncode": result.returncode,
            "stdout": result.stdout,
            "stderr": result.stderr,
        }
    except subprocess.TimeoutExpired:
        return {
            "test": test_path,
            "returncode": -1,
            "stdout": "",
            "stderr": "Test timed out after 5 minutes",
        }
    except Exception as e:
        return {
            "test": test_path,
            "returncode": -2,
            "stdout": "",
            "stderr": f"Error running test: {str(e)}",
        }


def main():
    """Main function to run tests"""
    # Change to project directory
    project_root = Path(__file__).parent
    print(f"Running from: {project_root}")

    # List of integration tests to run
    integration_tests = [
        "tests/integration/test_fastasize.py::test_fastasize_calculation",
        "tests/integration/test_chrprefix.py::test_numeric_chrs",
        "tests/integration/test_integration.py",
    ]

    results = []

    for test in integration_tests:
        print(f"\n=== Running {test} ===")
        result = run_test(test)
        results.append(result)

        print(f"Return code: {result['returncode']}")
        if result["stdout"]:
            print("STDOUT:")
            print(result["stdout"][:1000])  # First 1000 chars
        if result["stderr"]:
            print("STDERR:")
            print(result["stderr"][:1000])  # First 1000 chars

        if result["returncode"] == 0:
            print("✅ PASSED")
        else:
            print("❌ FAILED")

    # Summary
    print("\n=== SUMMARY ===")
    passed = sum(1 for r in results if r["returncode"] == 0)
    failed = len(results) - passed
    print(f"Passed: {passed}")
    print(f"Failed: {failed}")

    # Write detailed results to file
    with open("integration_test_results.txt", "w") as f:
        for result in results:
            f.write(f"\n=== {result['test']} ===\n")
            f.write(f"Return code: {result['returncode']}\n")
            f.write("STDOUT:\n")
            f.write(result["stdout"])
            f.write("\nSTDERR:\n")
            f.write(result["stderr"])
            f.write("\n" + "=" * 80 + "\n")


if __name__ == "__main__":
    main()
