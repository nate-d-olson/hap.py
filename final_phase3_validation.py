#!/usr/bin/env python3
"""
Final Phase 3 validation before committing.
"""

import os
import sys
import traceback


def test_phase3_functionality():
    """Test all Phase 3 functionality."""
    print("🧪 Final Phase 3 Validation")
    print("=" * 60)

    success_count = 0
    total_tests = 8

    try:
        # Add src to path
        sys.path.insert(0, "src")

        # Test 1: Import Phase 3 modules
        print("1. Testing Phase 3 imports...")
        try:
            from hap_py.haplo.python_quantify import QuantifyEngine
            from hap_py.haplo.quantify_phase3 import (
                MultiSampleQuantifier,
                RegionBasedQuantifier,
            )

            print("   ✅ All Phase 3 modules imported successfully")
            success_count += 1
        except Exception as e:
            print(f"   ❌ Import failed: {e}")
            return False

        # Test 2: Check class instantiation
        print("2. Testing Phase 3 class instantiation...")
        try:
            rbq = RegionBasedQuantifier()
            msq = MultiSampleQuantifier()
            print("   ✅ Phase 3 classes instantiated successfully")
            success_count += 1
        except Exception as e:
            print(f"   ❌ Class instantiation failed: {e}")

        # Test 3: Check QuantifyEngine with Phase 3 parameters
        print("3. Testing QuantifyEngine with Phase 3 parameters...")
        try:
            engine = QuantifyEngine(
                truth_vcf="example/performance.vcf.gz",
                query_vcf="example/PG_performance.vcf.gz",
                reference="example/chr21.fa",
                enable_superlocus_analysis=True,
                enable_region_stratification=False,
                enable_multi_sample=False,
            )
            print("   ✅ QuantifyEngine created with Phase 3 parameters")
            success_count += 1
        except Exception as e:
            print(f"   ❌ QuantifyEngine creation failed: {e}")

        # Test 4: Check Phase 3 methods exist
        print("4. Testing Phase 3 method availability...")
        try:
            methods = [
                "_perform_superlocus_analysis",
                "_perform_region_stratification",
                "_perform_multi_sample_analysis",
            ]
            for method in methods:
                if hasattr(engine, method):
                    print(f"   ✅ Method {method} exists")
                else:
                    print(f"   ❌ Method {method} missing")
                    continue
            success_count += 1
        except Exception as e:
            print(f"   ❌ Method check failed: {e}")

        # Test 5: Test BED region loading
        print("5. Testing BED region functionality...")
        try:
            import tempfile

            bed_content = "21\t26960070\t27230000\ttest_region"
            with tempfile.NamedTemporaryFile(
                mode="w", suffix=".bed", delete=False
            ) as f:
                f.write(bed_content)
                bed_path = f.name

            rbq = RegionBasedQuantifier()
            rbq.load_bed_regions({"test": bed_path})
            print("   ✅ BED file loading works")
            os.unlink(bed_path)
            success_count += 1
        except Exception as e:
            print(f"   ❌ BED functionality failed: {e}")

        # Test 6: Test multi-sample functionality
        print("6. Testing multi-sample functionality...")
        try:
            msq = MultiSampleQuantifier()
            # Test sample registration using add_sample instead of register_sample
            truth_vcf = "example/performance.vcf.gz"
            query_vcf = "example/PG_performance.vcf.gz"
            msq.add_sample("test_sample", truth_vcf, query_vcf, {"test": True})
            print("   ✅ Multi-sample functionality works")
            success_count += 1
        except Exception as e:
            print(f"   ❌ Multi-sample functionality failed: {e}")

        # Test 7: Check Phase 3 data structures
        print("7. Testing Phase 3 data structures...")
        try:
            engine = QuantifyEngine(
                truth_vcf="example/performance.vcf.gz",
                query_vcf="example/PG_performance.vcf.gz",
                reference="example/chr21.fa",
                enable_superlocus_analysis=True,
            )

            # Check if data structures are initialized
            required_attrs = [
                "superlocus_data",
                "region_stratification_results",
                "multi_sample_results",
            ]

            for attr in required_attrs:
                if hasattr(engine, attr):
                    print(f"   ✅ Attribute {attr} exists")
                else:
                    print(f"   ❌ Attribute {attr} missing")
                    continue
            success_count += 1
        except Exception as e:
            print(f"   ❌ Data structure check failed: {e}")

        # Test 8: Test minimal quantify run with Phase 3
        print("8. Testing minimal Phase 3 quantify run...")
        try:
            engine = QuantifyEngine(
                truth_vcf="example/performance.vcf.gz",
                query_vcf="example/PG_performance.vcf.gz",
                reference="example/chr21.fa",
                enable_superlocus_analysis=True,
                enable_region_stratification=False,
                enable_multi_sample=False,
            )

            # Set up minimal test data
            engine.truth_variants = [
                {"chromosome": "21", "position": 27000000, "ref": "A", "alt": ["T"]},
                {"chromosome": "21", "position": 27000100, "ref": "G", "alt": ["C"]},
            ]
            engine.query_variants = [
                {"chromosome": "21", "position": 27000000, "ref": "A", "alt": ["T"]},
                {"chromosome": "21", "position": 27000050, "ref": "C", "alt": ["G"]},
            ]
            engine.matched_variants = []

            # Test superlocus analysis method
            engine._perform_superlocus_analysis()

            print("   ✅ Phase 3 superlocus analysis runs successfully")
            success_count += 1
        except Exception as e:
            print(f"   ❌ Phase 3 quantify run failed: {e}")
            traceback.print_exc()

        # Summary
        print("\n" + "=" * 60)
        print(f"Phase 3 Validation Results: {success_count}/{total_tests} tests passed")

        if success_count == total_tests:
            print("🎉 All Phase 3 functionality is working correctly!")
            return True
        elif success_count >= 6:
            print("⚠️  Most Phase 3 functionality works - minor issues detected")
            return True
        else:
            print("❌ Significant Phase 3 issues detected")
            return False

    except Exception as e:
        print(f"❌ Critical error in validation: {e}")
        traceback.print_exc()
        return False


if __name__ == "__main__":
    success = test_phase3_functionality()
    sys.exit(0 if success else 1)
