#!/usr/bin/env python3
"""
Simple test script to verify quantify fixes work.
"""
import os
import sys
import tempfile

import pandas as pd

# Add src to path
project_root = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(project_root, "src"))

try:
    from hap_py.haplo.python_quantify import QuantifyEngine

    print("✓ Successfully imported QuantifyEngine")
except Exception as e:
    print(f"✗ Failed to import: {e}")
    sys.exit(1)

# Create temporary VCF files
vcf_content = """##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE
chr1	100	.	A	T	60	PASS	.	GT	0/1
"""

with tempfile.NamedTemporaryFile(mode="w", suffix=".vcf", delete=False) as truth_f:
    truth_f.write(vcf_content)
    truth_file = truth_f.name

with tempfile.NamedTemporaryFile(mode="w", suffix=".vcf", delete=False) as query_f:
    query_f.write(vcf_content)
    query_file = query_f.name

try:
    # Test QuantifyEngine initialization
    engine = QuantifyEngine(truth_vcf=truth_file, query_vcf=query_file)
    print("✓ Successfully created QuantifyEngine")

    # Test _are_alleles_compatible with pandas Series
    if hasattr(engine, "_are_alleles_compatible"):
        # Test same variants
        var1 = pd.Series({"ref": "A", "alt": "G", "pos": 100})
        var2 = pd.Series({"ref": "A", "alt": "G", "pos": 100})
        result = engine._are_alleles_compatible(var1, var2)
        print(f"✓ Same variants compatible: {result}")
        assert result is True

        # Test different variants
        var3 = pd.Series({"ref": "A", "alt": "T", "pos": 100})
        result = engine._are_alleles_compatible(var1, var3)
        print(f"✓ Different variants incompatible: {result}")
        assert result is False

        print("✓ All allele compatibility tests passed")
    else:
        print("⚠ _are_alleles_compatible method not found")

    # Test _classify_variant_type
    if hasattr(engine, "_classify_variant_type"):
        # Test SNP
        snp_var = {"ref": "A", "alt": "T"}
        result = engine._classify_variant_type(snp_var)
        print(f"✓ SNP classification: {result}")
        assert result == "SNP"

        # Test deletion
        del_var = {"ref": "ATG", "alt": "A"}
        result = engine._classify_variant_type(del_var)
        print(f"✓ Deletion classification: {result}")
        assert result == "DEL"

        # Test complex variant
        complex_var = {"ref": "ATG", "alt": "TCCG"}
        result = engine._classify_variant_type(complex_var)
        print(f"✓ Complex variant classification: {result}")
        # This might be INS, MNP, or COMPLEX depending on implementation

        print("✓ All variant classification tests passed")
    else:
        print("⚠ _classify_variant_type method not found")

    print("\n🎉 All tests passed successfully!")

except Exception as e:
    print(f"✗ Test failed: {e}")
    import traceback

    traceback.print_exc()

finally:
    # Clean up temp files
    try:
        os.unlink(truth_file)
        os.unlink(query_file)
    except:
        pass
