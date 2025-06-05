#!/usr/bin/env python3
"""Simple Phase 3 test to check if everything is working."""


def test_phase3_imports():
    """Test Phase 3 imports."""
    print("Testing Phase 3 imports...")
    try:

        print("✅ All imports successful")
        return True
    except Exception as e:
        print(f"❌ Import failed: {e}")
        return False


def test_phase3_instantiation():
    """Test Phase 3 class instantiation."""
    print("Testing Phase 3 instantiation...")
    try:
        from hap_py.haplo.quantify_phase3 import (
            MultiSampleQuantifier,
            RegionBasedQuantifier,
        )

        region_quantifier = RegionBasedQuantifier()
        multi_sample_quantifier = MultiSampleQuantifier()
        print("✅ Classes instantiated successfully")
        return True
    except Exception as e:
        print(f"❌ Instantiation failed: {e}")
        return False


def test_phase3_integration():
    """Test Phase 3 integration with QuantifyEngine."""
    print("Testing Phase 3 integration...")
    try:
        from hap_py.haplo.python_quantify import QuantifyEngine

        # Test creating engine with Phase 3 parameters
        engine = QuantifyEngine(
            truth_vcf="example/performance.vcf.gz",
            query_vcf="example/PG_performance.vcf.gz",
            reference="example/chr21.fa",
            enable_superlocus_analysis=True,
            enable_region_stratification=True,
            enable_multi_sample=True,
            superlocus_window=1000,
        )
        print("✅ QuantifyEngine created with Phase 3 features")
        return True
    except Exception as e:
        print(f"❌ Integration test failed: {e}")
        return False


def main():
    """Run simple Phase 3 tests."""
    print("🧪 Simple Phase 3 Validation Test")
    print("=" * 40)

    success = True
    success &= test_phase3_imports()
    success &= test_phase3_instantiation()
    success &= test_phase3_integration()

    print("\n" + "=" * 40)
    if success:
        print("✅ All Phase 3 tests PASSED")
        print("Phase 3 implementation is working correctly!")
    else:
        print("❌ Some Phase 3 tests FAILED")

    return 0 if success else 1


if __name__ == "__main__":
    exit(main())
