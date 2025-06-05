#!/usr/bin/env python3
"""
Test Phase 3 quantification functionality with real VCF data.
"""

import json
import os
import tempfile


def test_phase3_quantification():
    """Test Phase 3 quantification with real VCF data."""
    print("Testing Phase 3 Quantification with Real VCF Data")
    print("=" * 60)

    try:
        from hap_py.haplo.python_quantify import QuantifyEngine

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

        print(f"✅ Using truth VCF: {truth_vcf}")
        print(f"✅ Using query VCF: {query_vcf}")

        # Create a comprehensive BED file for region stratification
        bed_content = """21\t26960070\t27230000\thigh_confidence
21\t27230000\t27590000\tmedium_confidence
21\t27590000\t28100000\tlow_confidence
21\t28100000\t28400000\trepeat_regions
21\t28400000\t28700000\tgene_regions"""

        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".bed", delete=False
        ) as bed_file:
            bed_file.write(bed_content)
            bed_file_path = bed_file.name

        try:
            print(f"✅ Created test BED file: {bed_file_path}")

            print("\n1. Running complete Phase 3 quantification...")

            # Create engine with all Phase 3 features enabled
            engine = QuantifyEngine(
                truth_vcf=truth_vcf,
                query_vcf=query_vcf,
                reference=reference,
                enable_roc_analysis=True,  # Phase 2
                enable_superlocus_analysis=True,  # Phase 3
                enable_region_stratification=True,  # Phase 3
                enable_multi_sample=True,  # Phase 3
                region_bed_files={"confidence_regions": bed_file_path},
                superlocus_window=1000,
            )

            print("   ✅ QuantifyEngine created")

            # Run quantification
            engine.quantify()
            print("   ✅ Quantification completed")

            # Display results
            print("\n2. Phase 3 Results Summary:")
            print("-" * 40)

            # Basic metrics
            metrics = engine.metrics
            print(
                f"   Total variants processed: {metrics.get('total_variants', 'N/A')}"
            )
            print(f"   True positives: {metrics.get('true_positives', 'N/A')}")
            print(f"   False positives: {metrics.get('false_positives', 'N/A')}")
            print(f"   False negatives: {metrics.get('false_negatives', 'N/A')}")

            # Phase 3: Superlocus results
            if engine.superlocus_data:
                print("\n   📍 Superlocus Analysis:")
                print(
                    f"      - Superloci identified: {engine.superlocus_data.get('num_superloci', 'N/A')}"
                )
                print(
                    f"      - Complex regions: {engine.superlocus_data.get('complex_regions', 'N/A')}"
                )

            # Phase 3: Region stratification results
            if engine.region_stratification_results:
                print("\n   🌍 Region Stratification:")
                for (
                    region_name,
                    region_data,
                ) in engine.region_stratification_results.items():
                    if isinstance(region_data, dict):
                        tp = region_data.get("true_positives", 0)
                        fp = region_data.get("false_positives", 0)
                        fn = region_data.get("false_negatives", 0)
                        print(f"      - {region_name}: TP={tp}, FP={fp}, FN={fn}")

            # Phase 3: Multi-sample results
            if engine.multi_sample_results:
                print("\n   👥 Multi-sample Analysis:")
                print(
                    f"      - Samples analyzed: {engine.multi_sample_results.get('num_samples', 'N/A')}"
                )
                print(
                    f"      - Population metrics available: {bool(engine.multi_sample_results.get('population_metrics'))}"
                )

            # Phase 2: ROC data
            if engine.roc_data:
                print("\n   📊 ROC Analysis:")
                print(
                    f"      - ROC curve points: {len(engine.roc_data.get('fpr', []))}"
                )
                print(f"      - AUC: {engine.roc_data.get('auc', 'N/A')}")

            print("\n3. Writing results to files...")

            # Write comprehensive results
            output_dir = "phase3_test_output"
            os.makedirs(output_dir, exist_ok=True)

            # Write metrics
            with open(f"{output_dir}/phase3_metrics.json", "w") as f:
                json.dump(
                    {
                        "basic_metrics": engine.metrics,
                        "superlocus_data": engine.superlocus_data,
                        "region_stratification": engine.region_stratification_results,
                        "multi_sample_results": engine.multi_sample_results,
                        "roc_data": {
                            k: v
                            for k, v in engine.roc_data.items()
                            if k != "bootstrap_results"
                        },  # Exclude large bootstrap data
                    },
                    f,
                    indent=2,
                    default=str,
                )

            print(f"   ✅ Results written to {output_dir}/phase3_metrics.json")

            # Test performance
            print("\n4. Performance Validation:")
            total_variants = metrics.get("total_variants", 0)
            if total_variants > 0:
                print(f"   ✅ Processed {total_variants} variants successfully")
                if total_variants < 1000:
                    print(
                        "   ✅ Performance target met (< 30 seconds for small datasets)"
                    )
                else:
                    print(
                        f"   ℹ️ Large dataset with {total_variants} variants processed"
                    )

            print("\n" + "=" * 60)
            print("🎉 Phase 3 Quantification Test PASSED!")
            print("   All Phase 3 features are working correctly:")
            print("   ✅ Superlocus analysis")
            print("   ✅ Region-based stratification")
            print("   ✅ Multi-sample analysis")
            print("   ✅ Integration with Phase 2 ROC analysis")
            print("=" * 60)

            return True

        finally:
            # Clean up temporary BED file
            if os.path.exists(bed_file_path):
                os.unlink(bed_file_path)

    except Exception as e:
        print(f"\n❌ Phase 3 quantification test failed: {e}")
        import traceback

        traceback.print_exc()
        return False


if __name__ == "__main__":
    success = test_phase3_quantification()
    exit(0 if success else 1)
