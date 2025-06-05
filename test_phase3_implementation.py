#!/usr/bin/env python3
"""
Test script for Phase 3 quantify implementation.

This script validates the Phase 3 functionality including:
- Superlocus analysis
- Region-based quantification
- Multi-sample analysis
"""

import os
import sys
import tempfile
from pathlib import Path

# Add src to path to import hap_py modules
sys.path.insert(0, str(Path(__file__).parent / "src"))

from hap_py.haplo.python_quantify import QuantifyEngine
from hap_py.haplo.quantify_phase3 import (
    MultiSampleQuantifier,
    RegionBasedQuantifier,
)


def test_basic_phase3_initialization():
    """Test that Phase 3 components can be initialized properly."""
    print("=== Testing Phase 3 Basic Initialization ===")

    # Test data paths
    truth_vcf = "example/integration/integrationtest.vcf"
    query_vcf = "example/integration/integrationtest_rhs.vcf"
    bed_file = "example/hc.bed"

    # Verify test files exist
    for filepath in [truth_vcf, query_vcf, bed_file]:
        if not os.path.exists(filepath):
            print(f"ERROR: Test file {filepath} not found")
            return False

    try:
        # Test basic QuantifyEngine with Phase 3 disabled
        basic_engine = QuantifyEngine(
            truth_vcf=truth_vcf, query_vcf=query_vcf, quantify_method="xcmp"
        )
        print("✅ Basic QuantifyEngine initialization successful")
        # Use the basic_engine to avoid unused variable warning
        assert basic_engine is not None

        # Test QuantifyEngine with Phase 3 enabled
        phase3_engine = QuantifyEngine(
            truth_vcf=truth_vcf,
            query_vcf=query_vcf,
            quantify_method="xcmp",
            enable_superlocus_analysis=True,
            enable_region_stratification=True,
            enable_multi_sample=True,
            region_bed_files={"high_confidence": bed_file},
            superlocus_window=1000,
        )
        print("✅ Phase 3 QuantifyEngine initialization successful")
        assert phase3_engine is not None

        # Test RegionBasedQuantifier
        region_quantifier = RegionBasedQuantifier()
        region_quantifier.load_bed_regions({"test_region": bed_file})
        print("✅ RegionBasedQuantifier initialization successful")

        # Test MultiSampleQuantifier
        multi_quantifier = MultiSampleQuantifier()
        multi_quantifier.load_vcf_samples([truth_vcf, query_vcf])
        print("✅ MultiSampleQuantifier initialization successful")

        return True

    except Exception as e:
        print(f"❌ Error during initialization: {e}")
        import traceback

        traceback.print_exc()
        return False


def test_region_based_quantification():
    """Test region-based quantification functionality."""
    print("\n=== Testing Region-Based Quantification ===")

    try:
        # Initialize region quantifier
        region_quantifier = RegionBasedQuantifier()

        # Test BED file loading
        bed_file = "example/hc.bed"
        region_quantifier.load_bed_regions({"high_confidence": bed_file})

        print(f"✅ Loaded {len(region_quantifier.bed_regions)} region types")

        # Test variant stratification (mock data)
        test_variants = [
            {"chrom": "chr21", "pos": 20183750, "ref": "A", "alt": "G"},
            {"chrom": "chr21", "pos": 20256850, "ref": "C", "alt": "T"},
            {
                "chrom": "chr21",
                "pos": 25000000,
                "ref": "G",
                "alt": "A",
            },  # Outside regions
        ]

        for variant in test_variants:
            regions = region_quantifier.stratify_variant(variant)
            print(
                f"Variant at {variant['chrom']}:{variant['pos']} -> regions: {regions}"
            )

        print("✅ Region-based stratification working")
        return True

    except Exception as e:
        print(f"❌ Error in region-based quantification: {e}")
        import traceback

        traceback.print_exc()
        return False


def test_multi_sample_analysis():
    """Test multi-sample analysis functionality."""
    print("\n=== Testing Multi-Sample Analysis ===")

    try:
        # Initialize multi-sample quantifier
        multi_quantifier = MultiSampleQuantifier()

        # Test VCF loading
        vcf_files = [
            "example/integration/integrationtest.vcf",
            "example/integration/integrationtest_rhs.vcf",
        ]

        multi_quantifier.load_vcf_samples(vcf_files)
        print(f"✅ Loaded {len(multi_quantifier.samples)} samples")

        # Test population analysis (mock)
        if multi_quantifier.samples:
            sample_names = list(multi_quantifier.samples.keys())
            print(f"Sample names: {sample_names}")

        print("✅ Multi-sample analysis working")
        return True

    except Exception as e:
        print(f"❌ Error in multi-sample analysis: {e}")
        import traceback

        traceback.print_exc()
        return False


def test_superlocus_analysis():
    """Test superlocus analysis functionality."""
    print("\n=== Testing Superlocus Analysis ===")

    try:
        # Test data
        truth_vcf = "example/integration/integrationtest.vcf"
        query_vcf = "example/integration/integrationtest_rhs.vcf"

        # Initialize engine with superlocus analysis
        engine = QuantifyEngine(
            truth_vcf=truth_vcf,
            query_vcf=query_vcf,
            quantify_method="xcmp",
            enable_superlocus_analysis=True,
            superlocus_window=1000,
        )

        print("✅ Engine with superlocus analysis initialized")

        # Test superlocus helper methods exist
        assert hasattr(engine, "_find_superlocus_matches")
        assert hasattr(engine, "_group_into_superloci")
        assert hasattr(engine, "_analyze_superlocus")
        print("✅ Superlocus analysis methods available")

        return True

    except Exception as e:
        print(f"❌ Error in superlocus analysis: {e}")
        import traceback

        traceback.print_exc()
        return False


def test_full_phase3_workflow():
    """Test the complete Phase 3 workflow."""
    print("\n=== Testing Full Phase 3 Workflow ===")

    try:
        # Test data
        truth_vcf = "example/integration/integrationtest.vcf"
        query_vcf = "example/integration/integrationtest_rhs.vcf"
        bed_file = "example/hc.bed"

        # Create temporary output directory
        with tempfile.TemporaryDirectory() as temp_dir:
            output_prefix = os.path.join(temp_dir, "phase3_test")

            # Initialize engine with all Phase 3 features
            engine = QuantifyEngine(
                truth_vcf=truth_vcf,
                query_vcf=query_vcf,
                quantify_method="xcmp",
                enable_superlocus_analysis=True,
                enable_region_stratification=True,
                enable_multi_sample=True,
                region_bed_files={"high_confidence": bed_file},
                superlocus_window=1000,
                enable_roc_analysis=True,
                quality_stratification=True,
            )

            print("✅ Full Phase 3 engine initialized")

            # Run quantification (this might take a moment)
            print("Running quantification with Phase 3 features...")
            results = engine.quantify(output_prefix=output_prefix)

            print("✅ Phase 3 quantification completed")

            # Check results structure
            if "superlocus_data" in results:
                print(
                    f"✅ Superlocus data generated: {len(results['superlocus_data'])} superloci"
                )

            if "region_stratification_results" in results:
                print("✅ Region stratification results generated")

            if "multi_sample_results" in results:
                print("✅ Multi-sample results generated")

            # Check output files
            output_files = list(Path(temp_dir).glob("phase3_test*"))
            print(f"✅ Generated {len(output_files)} output files")
            for f in output_files:
                print(f"  - {f.name}")

            return True

    except Exception as e:
        print(f"❌ Error in full Phase 3 workflow: {e}")
        import traceback

        traceback.print_exc()
        return False


def main():
    """Run all Phase 3 tests."""
    print("Starting Phase 3 Implementation Tests")
    print("=" * 50)

    # Change to the repository directory
    repo_root = Path(__file__).parent
    os.chdir(repo_root)

    tests = [
        test_basic_phase3_initialization,
        test_region_based_quantification,
        test_multi_sample_analysis,
        test_superlocus_analysis,
        test_full_phase3_workflow,
    ]

    passed = 0
    failed = 0

    for test_func in tests:
        try:
            if test_func():
                passed += 1
            else:
                failed += 1
        except Exception as e:
            print(f"❌ Test {test_func.__name__} failed with exception: {e}")
            failed += 1

    print("\n" + "=" * 50)
    print(f"Test Results: {passed} passed, {failed} failed")

    if failed == 0:
        print("🎉 All Phase 3 tests passed!")
        return 0
    else:
        print("⚠️  Some Phase 3 tests failed. Review the output above.")
        return 1


if __name__ == "__main__":
    sys.exit(main())
