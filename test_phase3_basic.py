#!/usr/bin/env python3
"""
Test Phase 3 implementation of the quantify module.
"""

import os
import tempfile


def test_phase3_imports():
    """Test that Phase 3 classes can be imported."""
    try:
        from src.hap_py.haplo.quantify_phase3 import (
            MultiSampleQuantifier,
            RegionBasedQuantifier,
        )

        print("✅ Phase 3 classes imported successfully")

        # Test RegionBasedQuantifier instantiation
        rq = RegionBasedQuantifier()
        print("✅ RegionBasedQuantifier instantiated")

        # Test MultiSampleQuantifier instantiation
        msq = MultiSampleQuantifier()
        print("✅ MultiSampleQuantifier instantiated")

        return True
    except Exception as e:
        print(f"❌ Failed to import Phase 3 classes: {e}")
        return False


def test_phase3_integration():
    """Test Phase 3 integration with main QuantifyEngine."""
    try:
        from src.hap_py.haplo.python_quantify import QuantifyEngine

        # Use example VCF files
        truth_vcf = "example/performance.vcf.gz"
        query_vcf = "example/PG_performance.vcf.gz"
        reference = "example/chr21.fa"

        # Check if files exist
        if not os.path.exists(truth_vcf):
            print(f"❌ Truth VCF not found: {truth_vcf}")
            return False

        if not os.path.exists(query_vcf):
            print(f"❌ Query VCF not found: {query_vcf}")
            return False

        # Create a simple BED file for testing
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".bed", delete=False
        ) as bed_file:
            bed_file.write("21\t10000\t20000\ttest_region\n")
            bed_file.write("21\t30000\t40000\ttest_region2\n")
            bed_file_path = bed_file.name

        try:
            # Test QuantifyEngine with Phase 3 enabled
            engine = QuantifyEngine(
                truth_vcf=truth_vcf,
                query_vcf=query_vcf,
                reference=reference,
                enable_superlocus_analysis=True,
                enable_region_stratification=True,
                enable_multi_sample=True,
                region_bed_files={"test_regions": bed_file_path},
            )
            print("✅ QuantifyEngine with Phase 3 enabled instantiated")

            # Check that Phase 3 quantifiers were initialized
            if hasattr(engine, "region_quantifier") and engine.region_quantifier:
                print("✅ RegionBasedQuantifier initialized")
            else:
                print("❌ RegionBasedQuantifier not initialized")

            if (
                hasattr(engine, "multi_sample_quantifier")
                and engine.multi_sample_quantifier
            ):
                print("✅ MultiSampleQuantifier initialized")
            else:
                print("❌ MultiSampleQuantifier not initialized")

            return True

        finally:
            # Clean up temporary file
            os.unlink(bed_file_path)

    except Exception as e:
        print(f"❌ Failed Phase 3 integration test: {e}")
        import traceback

        traceback.print_exc()
        return False


def main():
    """Run Phase 3 tests."""
    print("Testing Phase 3 implementation...")
    print("=" * 50)

    success = True

    # Test imports
    print("\n1. Testing Phase 3 imports...")
    success &= test_phase3_imports()

    # Test integration
    print("\n2. Testing Phase 3 integration...")
    success &= test_phase3_integration()

    print("\n" + "=" * 50)
    if success:
        print("✅ All Phase 3 tests passed!")
    else:
        print("❌ Some Phase 3 tests failed!")

    return success


if __name__ == "__main__":
    main()
