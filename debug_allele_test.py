#!/usr/bin/env python3
import os
import sys
import tempfile

import pandas as pd

sys.path.append("src")
from hap_py.haplo.python_quantify import QuantifyEngine

# Create temp files
with tempfile.NamedTemporaryFile(mode="w", suffix=".vcf", delete=False) as f:
    f.write(
        '##fileformat=VCFv4.2\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\nchr1\t100\t.\tA\tT\t60\tPASS\t.\tGT\t0/1\n'
    )
    truth_file = f.name

with tempfile.NamedTemporaryFile(mode="w", suffix=".vcf", delete=False) as f:
    f.write(
        '##fileformat=VCFv4.2\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\nchr1\t100\t.\tA\tT\t60\tPASS\t.\tGT\t0/1\n'
    )
    query_file = f.name

try:
    engine = QuantifyEngine(truth_vcf=truth_file, query_vcf=query_file)

    # Test with pandas Series
    var1 = pd.Series({"ref": "A", "alt": "G", "pos": 100})
    var2 = pd.Series({"ref": "A", "alt": "G", "pos": 100})
    result = engine._are_alleles_compatible(var1, var2)
    print(f"Same SNP test: {result}")

    var3 = pd.Series({"ref": "A", "alt": "T", "pos": 100})
    result = engine._are_alleles_compatible(var1, var3)
    print(f"Different SNP test: {result}")

    # Test the classification function
    variant = {"ref": "ATG", "alt": "TCCG"}
    result = engine._classify_variant_type(variant)
    print(f"Complex variant classification: {result}")

finally:
    os.unlink(truth_file)
    os.unlink(query_file)
