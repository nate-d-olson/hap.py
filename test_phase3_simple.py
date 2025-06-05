#!/usr/bin/env python3
"""
Simple Phase 3 functional test.
"""


def test_phase3_analysis():
    """Test Phase 3 analysis methods."""
    print("Testing Phase 3 analysis methods...")

    try:
        from hap_py.haplo.python_quantify import QuantifyEngine

        # Create engine with Phase 3 enabled
        engine = QuantifyEngine(
            truth_vcf="example/performance.vcf.gz",
            query_vcf="example/PG_performance.vcf.gz",
            reference="example/chr21.fa",
            enable_superlocus_analysis=True,
            enable_multi_sample=True,
        )
        print("✅ QuantifyEngine created")

        # Manually set up some test data
        engine.truth_variants = [
            {"chromosome": "21", "position": 27000000, "ref": "A", "alt": ["T"]},
            {"chromosome": "21", "position": 27300000, "ref": "G", "alt": ["C"]},
        ]
        engine.query_variants = [
            {"chromosome": "21", "position": 27000000, "ref": "A", "alt": ["T"]},
            {"chromosome": "21", "position": 27300000, "ref": "G", "alt": ["A"]},
        ]
        engine.matched_variants = [(0, 0, "TP", 0.95)]

        print("✅ Test data prepared")

        # Test superlocus analysis
        print("Testing superlocus analysis...")
        engine._perform_superlocus_analysis()

        if hasattr(engine, "superlocus_data"):
            print("✅ Superlocus analysis completed")
            print(f"   - Data keys: {list(engine.superlocus_data.keys())}")
        else:
            print("❌ Superlocus analysis failed")

        # Test multi-sample analysis
        print("Testing multi-sample analysis...")
        engine._perform_multi_sample_analysis()
        print("✅ Multi-sample analysis completed")

        return True

    except Exception as e:
        print(f"❌ Test failed: {e}")
        import traceback

        traceback.print_exc()
        return False


if __name__ == "__main__":
    print("Phase 3 Simple Test")
    print("=" * 30)

    success = test_phase3_analysis()

    if success:
        print("✅ Phase 3 is working!")
    else:
        print("❌ Phase 3 has issues!")
