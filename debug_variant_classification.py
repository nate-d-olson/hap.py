#!/usr/bin/env python3

import sys
import tempfile
from pathlib import Path

# Add src to sys.path to allow importing hap_py
project_root_path = Path(__file__).resolve().parent
sys.path.insert(0, str(project_root_path / "src"))

from hap_py.haplo.python_quantify import QuantifyEngine


class MockVariant:
    def __init__(self, ref, alts=None):
        self.ref = ref
        self.alts = alts


def test_variant_classification():
    """Debug variant classification."""

    # Create temporary VCF files
    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".vcf", delete=False
    ) as truth_f, tempfile.NamedTemporaryFile(
        mode="w", suffix=".vcf", delete=False
    ) as query_f:

        # Write minimal valid VCF content
        vcf_content = """##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE
chr1	100	.	A	T	60	PASS	.	GT	0/1
"""
        truth_f.write(vcf_content)
        query_f.write(vcf_content)
        truth_f.flush()
        query_f.flush()

        try:
            engine = QuantifyEngine(truth_vcf=truth_f.name, query_vcf=query_f.name)

            # Test the specific case that's failing
            mnp = MockVariant("AT", ["GC"])
            result = engine._get_variant_type(mnp)
            print(
                f"MNP test: ref='AT', alts=['GC'] -> result='{result}' (expected='MNP')"
            )

            # Test _classify_variant_type with dict input
            if hasattr(engine, "_classify_variant_type"):
                mnp_variant = {"ref": "ATG", "alt": "TCC"}
                result2 = engine._classify_variant_type(mnp_variant)
                print(
                    f"_classify_variant_type test: {mnp_variant} -> result='{result2}' (expected='COMPLEX')"
                )

        except Exception as e:
            print(f"Error: {e}")
            import traceback

            traceback.print_exc()
        finally:
            import os

            os.unlink(truth_f.name)
            os.unlink(query_f.name)


if __name__ == "__main__":
    test_variant_classification()
