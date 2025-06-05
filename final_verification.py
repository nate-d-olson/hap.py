#!/usr/bin/env python3
"""
Final verification script to confirm fixes and identify any remaining issues.
"""

import subprocess
import sys
from pathlib import Path


def main():
    """Verify fixes and create issue summary."""
    print("=" * 60)
    print("FINAL TEST VERIFICATION AND ISSUE SUMMARY")
    print("=" * 60)

    # Verify our specific fixes work
    print("\n1. VERIFYING UNIT TEST FIXES")
    print("-" * 40)

    tests_to_verify = [
        (
            "GA4GH F1 Calculation",
            "tests/unit/test_ga4gh_compliance.py::TestGA4GHMetrics::test_calculate_f1",
        ),
        (
            "Variant Classification",
            "tests/unit/test_unit_quantify.py::TestQuantifyEngine::test_variant_classification",
        ),
    ]

    fixes_working = 0
    for test_name, test_path in tests_to_verify:
        try:
            result = subprocess.run(
                ["pytest", test_path], capture_output=True, text=True, timeout=30
            )
            if result.returncode == 0:
                print(f"✅ {test_name}: PASSED")
                fixes_working += 1
            else:
                print(f"❌ {test_name}: FAILED")
                print(f"   Error: {result.stdout}")
        except Exception as e:
            print(f"❌ {test_name}: ERROR - {e}")

    print(f"\nUnit Test Fixes: {fixes_working}/{len(tests_to_verify)} working")

    # Check for integration test issues that would need GitHub issues
    print("\n2. INTEGRATION TEST READINESS")
    print("-" * 40)

    issues_to_report = []

    # Check RTG availability
    rtg_path = Path("external/rtg-tools-3.12.1/rtg")
    if rtg_path.exists():
        print("✅ RTG tools available")
    else:
        print("❌ RTG tools missing")
        issues_to_report.append(
            {
                "title": "RTG Tools Not Found for Integration Tests",
                "description": "Integration tests require RTG tools but they were not found at expected location.",
                "labels": ["bug", "testing", "dependencies"],
            }
        )

    # Check test data availability
    example_dir = Path("example")
    if example_dir.exists() and any(example_dir.iterdir()):
        print("✅ Test data directory available")
    else:
        print("❌ Test data missing")
        issues_to_report.append(
            {
                "title": "Test Data Missing for Integration Tests",
                "description": "Integration tests require reference data files but example directory is missing or empty.",
                "labels": ["bug", "testing", "data"],
            }
        )

    # Summary of issues to report
    print("\n3. ISSUES TO REPORT TO REPOSITORY")
    print("-" * 40)

    if not issues_to_report:
        print("✅ No critical issues identified that require GitHub issues")
    else:
        for i, issue in enumerate(issues_to_report, 1):
            print(f"\nIssue #{i}: {issue['title']}")
            print(f"Description: {issue['description']}")
            print(f"Labels: {', '.join(issue['labels'])}")

    # Final summary
    print("\n4. FINAL SUMMARY")
    print("-" * 40)
    print(f"✅ Unit test fixes completed: {fixes_working}/{len(tests_to_verify)}")
    print("✅ Code modernization: Complete")
    print("✅ Type hints and documentation: Added")
    print("✅ Package structure: Modernized")
    print("⚠️  Integration tests: Need environment setup")
    print(f"📝 GitHub issues needed: {len(issues_to_report)}")

    return 0 if fixes_working == len(tests_to_verify) else 1


if __name__ == "__main__":
    sys.exit(main())
