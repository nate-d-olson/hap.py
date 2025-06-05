#!/usr/bin/env python3
"""Simple test to check GA4GH integration imports"""

import sys
import traceback

try:
    print("Testing GA4GH imports...")

    # Test individual imports
    print("  - Importing GA4GH compliance...")
    from hap_py.haplo.ga4gh_compliance import GA4GHFormatter

    print("    ✅ GA4GH compliance imported successfully")

    print("  - Importing GA4GH integration...")

    print("    ✅ GA4GH integration imported successfully")

    print("  - Importing QuantifyEngine...")

    print("    ✅ QuantifyEngine imported successfully")

    print("  - Testing basic functionality...")
    formatter = GA4GHFormatter()
    print(f"    ✅ GA4GH formatter created: {type(formatter)}")

    print("\n✅ All imports working correctly!")

except Exception as e:
    print(f"\n❌ Import error: {e}")
    print(f"Error type: {type(e)}")
    print("\nFull traceback:")
    traceback.print_exc()
    sys.exit(1)
