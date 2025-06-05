#!/usr/bin/env python3
"""
Debug script for MultiSampleQuantifier class.
"""

import sys

sys.path.insert(0, "src")

try:
    from hap_py.haplo.quantify_phase3 import MultiSampleQuantifier

    print("MultiSampleQuantifier methods:")
    msq = MultiSampleQuantifier()
    print(dir(msq))

    print("\nChecking for register_sample method:")
    if hasattr(msq, "register_sample"):
        print("✅ register_sample method exists")
    else:
        print("❌ register_sample method missing")

    print("\nPrinting class source code location:")
    print(MultiSampleQuantifier.__module__)

except Exception as e:
    print(f"Error: {e}")
