#!/usr/bin/env python3
"""
Test Status Summary - Quick analysis of hap.py test status
"""

import subprocess
import sys
from pathlib import Path


def run_command_with_timeout(cmd, timeout=60):
    """Run command with timeout and capture output"""
    try:
        result = subprocess.run(
            cmd, shell=True, capture_output=True, text=True, timeout=timeout
        )
        return {
            "success": result.returncode == 0,
            "returncode": result.returncode,
            "stdout": result.stdout,
            "stderr": result.stderr,
            "timed_out": False,
        }
    except subprocess.TimeoutExpired:
        return {
            "success": False,
            "returncode": -1,
            "stdout": "",
            "stderr": f"Command timed out after {timeout} seconds",
            "timed_out": True,
        }
    except Exception as e:
        return {
            "success": False,
            "returncode": -2,
            "stdout": "",
            "stderr": f"Error: {str(e)}",
            "timed_out": False,
        }


def main():
    print("=== HAP.PY TEST STATUS SUMMARY ===\n")

    # Check environment
    print("1. Environment Check:")
    env_result = run_command_with_timeout("python --version && which python")
    if env_result["success"]:
        print(f"   ✅ Python: {env_result['stdout'].strip()}")
    else:
        print(f"   ❌ Python check failed: {env_result['stderr']}")

    # Check package installation
    print("\n2. Package Installation:")
    import_result = run_command_with_timeout(
        "python -c 'import hap_py; print(hap_py.__file__)'"
    )
    if import_result["success"]:
        print("   ✅ hap_py imported successfully")
        print(f"   Location: {import_result['stdout'].strip()}")
    else:
        print(f"   ❌ Import failed: {import_result['stderr']}")

    # Run unit tests with timeout
    print("\n3. Unit Tests (60s timeout):")
    unit_result = run_command_with_timeout(
        "python -m pytest tests/unit/ --tb=no -q", timeout=60
    )
    if unit_result["success"]:
        print("   ✅ Unit tests passed")
        print(f"   Output: {unit_result['stdout'].strip()}")
    else:
        print("   ❌ Unit tests failed or timed out")
        if unit_result["timed_out"]:
            print("   Reason: Timed out after 60 seconds")
        else:
            print(f"   Return code: {unit_result['returncode']}")
            if unit_result["stderr"]:
                print(f"   Error: {unit_result['stderr'][:500]}")

    # Test a simple integration test
    print("\n4. Simple Integration Test (30s timeout):")
    simple_test = run_command_with_timeout(
        "python -m pytest tests/integration/test_fastasize.py -v", timeout=30
    )
    if simple_test["success"]:
        print("   ✅ Simple integration test passed")
    else:
        print("   ❌ Simple integration test failed or timed out")
        if simple_test["timed_out"]:
            print("   Reason: Timed out after 30 seconds")
        else:
            print(f"   Return code: {simple_test['returncode']}")

    # List available tests
    print("\n5. Available Tests:")
    unit_files = list(Path("tests/unit").glob("*.py"))
    integration_files = list(Path("tests/integration").glob("*.py"))
    print(f"   Unit test files: {len(unit_files)}")
    print(f"   Integration test files: {len(integration_files)}")

    # Summary
    print("\n=== SUMMARY ===")
    tests_passed = 0
    tests_total = 4

    if env_result["success"]:
        tests_passed += 1
    if import_result["success"]:
        tests_passed += 1
    if unit_result["success"]:
        tests_passed += 1
    if simple_test["success"]:
        tests_passed += 1

    print(f"Tests passed: {tests_passed}/{tests_total}")

    if tests_passed == tests_total:
        print("🎉 All basic tests are working!")
    elif tests_passed >= 2:
        print("⚠️  Some tests working, but issues with integration tests")
    else:
        print("🚨 Major issues detected - environment or package problems")

    return tests_passed == tests_total


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
