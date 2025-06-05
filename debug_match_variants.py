#!/usr/bin/env python3
"""Debug script to test _match_variants method."""

import sys
import tempfile
import traceback

sys.path.insert(0, "src")

try:
    from hap_py.haplo.python_quantify import QuantifyEngine

    print("✓ QuantifyEngine imported successfully")

    # Create simple test data
    truth_variants = [{"chrom": "chr1", "pos": 100, "ref": "A", "alt": "G"}]
    query_variants = [{"chrom": "chr1", "pos": 100, "ref": "A", "alt": "G"}]

    # Create temporary VCF file
    with tempfile.NamedTemporaryFile(mode="w", suffix=".vcf", delete=False) as tf:
        tf.write(
            "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t100\t.\tA\tG\t60\tPASS\t.\n"
        )
        tf.flush()

        print("✓ Temporary VCF file created")

        engine = QuantifyEngine(truth_vcf=tf.name, query_vcf=tf.name)
        print("✓ QuantifyEngine created successfully")

        engine.truth_variants = truth_variants
        engine.query_variants = query_variants
        print("✓ Variants assigned")

        if hasattr(engine, "_match_variants"):
            print("✓ _match_variants method exists")
            print("Testing _match_variants...")
            engine._match_variants()
            print("✓ _match_variants executed successfully")

            # Check results
            print(f"Truth variants after matching: {engine.truth_variants}")
            print(f"Query variants after matching: {engine.query_variants}")
        else:
            print("✗ _match_variants method does not exist")

except Exception as e:
    print(f"✗ Error: {e}")
    traceback.print_exc()
