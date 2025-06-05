#!/usr/bin/env python3

import os

print("1. Script starting...")

try:
    print("2. Importing QuantifyEngine...")
    from hap_py.haplo.python_quantify import QuantifyEngine

    print("3. Import successful")

    # Check files
    truth_vcf = "example/performance.vcf.gz"
    query_vcf = "example/PG_performance.vcf.gz"

    print("4. Checking files...")
    print(f"   Truth VCF exists: {os.path.exists(truth_vcf)}")
    print(f"   Query VCF exists: {os.path.exists(query_vcf)}")

    if os.path.exists(truth_vcf) and os.path.exists(query_vcf):
        print("5. Creating QuantifyEngine...")
        engine = QuantifyEngine(
            truth_vcf=truth_vcf,
            query_vcf=query_vcf,
            enable_superlocus_analysis=True,
            enable_multi_sample=True,
        )
        print("6. QuantifyEngine created successfully!")
        print(f"7. Superlocus enabled: {engine.enable_superlocus_analysis}")
        print(f"8. Multi-sample enabled: {engine.enable_multi_sample}")
    else:
        print("5. Skipping QuantifyEngine creation - files not found")

    print("9. Test completed successfully!")

except Exception as e:
    print(f"Error: {e}")
    import traceback

    traceback.print_exc()
