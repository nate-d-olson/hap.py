#!/usr/bin/env python3
"""
Test import of quantify_phase3.py
"""

import traceback

print("Starting import test...")

try:
    print("Attempting to import module...")
    import hap_py.haplo.quantify_phase3 as mod

    print("Module imported successfully")
    print(f"Module contents: {dir(mod)}")

    # Try to access the classes
    if hasattr(mod, "RegionBasedQuantifier"):
        print("RegionBasedQuantifier found")
        print(f"RegionBasedQuantifier class: {mod.RegionBasedQuantifier}")
    else:
        print("RegionBasedQuantifier NOT found")

    if hasattr(mod, "MultiSampleQuantifier"):
        print("MultiSampleQuantifier found")
        print(f"MultiSampleQuantifier class: {mod.MultiSampleQuantifier}")
    else:
        print("MultiSampleQuantifier NOT found")

    # Let's try to exec the file directly to see what happens
    print("\nTrying to exec the file directly...")
    exec_ns = {}
    with open("src/hap_py/haplo/quantify_phase3.py") as f:
        file_content = f.read()
    exec(file_content, exec_ns)
    print(f"Exec namespace keys: {list(exec_ns.keys())}")

except Exception as e:
    print(f"Import error: {e}")
    traceback.print_exc()

print("Import test completed.")
