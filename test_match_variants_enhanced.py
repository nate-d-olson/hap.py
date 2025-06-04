#!/usr/bin/env python3
"""
Test script to verify the enhanced _match_variants() implementation.

This script tests the Phase 1 enhancements to the quantify module,
specifically the sophisticated variant matching algorithms.
"""

import logging
import os
import tempfile

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def create_test_vcf(variants, filename):
    """Create a test VCF file with the given variants."""
    vcf_content = """##fileformat=VCFv4.2
##INFO=<ID=BD,Number=1,Type=String,Description="Benchmarking Decision">
##INFO=<ID=BVT,Number=1,Type=String,Description="Benchmarking Variant Type">
##INFO=<ID=QQ,Number=1,Type=Integer,Description="Quality Quantile">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE
"""

    for variant in variants:
        line = f"{variant['chrom']}\t{variant['pos']}\t{variant.get('id', '.')}\t{variant['ref']}\t{variant['alt']}\t{variant.get('qual', 60)}\tPASS\t.\tGT:GQ\t{variant.get('gt', '0/1')}:60\n"
        vcf_content += line

    with open(filename, "w") as f:
        f.write(vcf_content)


def create_test_reference(filename):
    """Create a simple test reference file."""
    ref_content = """>chr1
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
"""
    with open(filename, "w") as f:
        f.write(ref_content)


def test_basic_matching():
    """Test basic exact matching functionality."""
    logger.info("Testing basic exact matching...")

    # Import after ensuring the environment is set up
    import sys

    sys.path.insert(0, "/Users/nolson/hap.py-modern-claude4/hap.py/src")

    from hap_py.haplo.python_quantify import QuantifyEngine

    with tempfile.TemporaryDirectory() as tmpdir:
        # Create test VCF files
        truth_variants = [
            {"chrom": "chr1", "pos": 100, "ref": "A", "alt": "T", "qual": 60},
            {"chrom": "chr1", "pos": 200, "ref": "GG", "alt": "G", "qual": 50},
            {"chrom": "chr1", "pos": 300, "ref": "C", "alt": "CTA", "qual": 40},
        ]

        query_variants = [
            {
                "chrom": "chr1",
                "pos": 100,
                "ref": "A",
                "alt": "T",
                "qual": 55,
            },  # Exact match
            {
                "chrom": "chr1",
                "pos": 200,
                "ref": "GG",
                "alt": "G",
                "qual": 45,
            },  # Exact match
            {
                "chrom": "chr1",
                "pos": 400,
                "ref": "T",
                "alt": "G",
                "qual": 30,
            },  # False positive
        ]

        truth_vcf = os.path.join(tmpdir, "truth.vcf")
        query_vcf = os.path.join(tmpdir, "query.vcf")
        ref_fasta = os.path.join(tmpdir, "ref.fa")

        create_test_vcf(truth_variants, truth_vcf)
        create_test_vcf(query_variants, query_vcf)
        create_test_reference(ref_fasta)

        # Test XCMP method
        logger.info("Testing XCMP method...")
        engine_xcmp = QuantifyEngine(
            truth_vcf=truth_vcf,
            query_vcf=query_vcf,
            reference=ref_fasta,
            quantify_method="xcmp",
        )

        results_xcmp = engine_xcmp.quantify()

        logger.info(
            f"XCMP Results: TP={results_xcmp['metrics']['TP']}, "
            f"FP={results_xcmp['metrics']['FP']}, "
            f"FN={results_xcmp['metrics']['FN']}"
        )

        # Test GA4GH method
        logger.info("Testing GA4GH method...")
        engine_ga4gh = QuantifyEngine(
            truth_vcf=truth_vcf,
            query_vcf=query_vcf,
            reference=ref_fasta,
            quantify_method="ga4gh",
        )

        results_ga4gh = engine_ga4gh.quantify()

        logger.info(
            f"GA4GH Results: TP={results_ga4gh['metrics']['TP']}, "
            f"FP={results_ga4gh['metrics']['FP']}, "
            f"FN={results_ga4gh['metrics']['FN']}"
        )

        # Verify results
        expected_tp = 2  # Two exact matches
        expected_fp = 1  # One false positive in query
        expected_fn = 1  # One false negative in truth

        assert (
            results_xcmp["metrics"]["TP"] == expected_tp
        ), f"XCMP TP mismatch: expected {expected_tp}, got {results_xcmp['metrics']['TP']}"
        assert (
            results_xcmp["metrics"]["FP"] == expected_fp
        ), f"XCMP FP mismatch: expected {expected_fp}, got {results_xcmp['metrics']['FP']}"
        assert (
            results_xcmp["metrics"]["FN"] == expected_fn
        ), f"XCMP FN mismatch: expected {expected_fn}, got {results_xcmp['metrics']['FN']}"

        logger.info("Basic matching test PASSED!")


def test_multiallelic_matching():
    """Test multi-allelic variant matching."""
    logger.info("Testing multi-allelic variant matching...")

    import sys

    sys.path.insert(0, "/Users/nolson/hap.py-modern-claude4/hap.py/src")

    from hap_py.haplo.python_quantify import QuantifyEngine

    with tempfile.TemporaryDirectory() as tmpdir:
        # Create test VCF files with multi-allelic variants
        truth_variants = [
            {
                "chrom": "chr1",
                "pos": 100,
                "ref": "A",
                "alt": "T,G",
                "qual": 60,
            },  # Multi-allelic
        ]

        query_variants = [
            {
                "chrom": "chr1",
                "pos": 100,
                "ref": "A",
                "alt": "T",
                "qual": 55,
            },  # Should match first allele
            {
                "chrom": "chr1",
                "pos": 100,
                "ref": "A",
                "alt": "G",
                "qual": 50,
            },  # Should match second allele
        ]

        truth_vcf = os.path.join(tmpdir, "truth.vcf")
        query_vcf = os.path.join(tmpdir, "query.vcf")
        ref_fasta = os.path.join(tmpdir, "ref.fa")

        create_test_vcf(truth_variants, truth_vcf)
        create_test_vcf(query_variants, query_vcf)
        create_test_reference(ref_fasta)

        # Test with XCMP method
        engine = QuantifyEngine(
            truth_vcf=truth_vcf,
            query_vcf=query_vcf,
            reference=ref_fasta,
            quantify_method="xcmp",
        )

        results = engine.quantify()

        logger.info(
            f"Multi-allelic Results: TP={results['metrics']['TP']}, "
            f"FP={results['metrics']['FP']}, "
            f"FN={results['metrics']['FN']}"
        )

        # Should have good matching performance
        assert (
            results["metrics"]["TP"] > 0
        ), "Should have some true positives from multi-allelic matching"

        logger.info("Multi-allelic matching test PASSED!")


def test_benchmarking_decisions():
    """Test benchmarking decision tracking (BD, BVT, QQ fields)."""
    logger.info("Testing benchmarking decision tracking...")

    import sys

    sys.path.insert(0, "/Users/nolson/hap.py-modern-claude4/hap.py/src")

    from hap_py.haplo.python_quantify import QuantifyEngine

    with tempfile.TemporaryDirectory() as tmpdir:
        truth_variants = [
            {"chrom": "chr1", "pos": 100, "ref": "A", "alt": "T", "qual": 60},
        ]

        query_variants = [
            {"chrom": "chr1", "pos": 100, "ref": "A", "alt": "T", "qual": 55},
        ]

        truth_vcf = os.path.join(tmpdir, "truth.vcf")
        query_vcf = os.path.join(tmpdir, "query.vcf")

        create_test_vcf(truth_variants, truth_vcf)
        create_test_vcf(query_variants, query_vcf)

        engine = QuantifyEngine(
            truth_vcf=truth_vcf,
            query_vcf=query_vcf,
            quantify_method="xcmp",
            output_vtc=True,
        )

        results = engine.quantify()

        # Check that BD, BVT, QQ fields are set
        truth_variant = engine.truth_variants[0]
        query_variant = engine.query_variants[0]

        assert "BD" in truth_variant, "Truth variant should have BD field"
        assert "BVT" in truth_variant, "Truth variant should have BVT field"
        assert "BD" in query_variant, "Query variant should have BD field"
        assert "BVT" in query_variant, "Query variant should have BVT field"
        assert "QQ" in query_variant, "Query variant should have QQ field"

        logger.info(f"Truth BD: {truth_variant['BD']}, BVT: {truth_variant['BVT']}")
        logger.info(
            f"Query BD: {query_variant['BD']}, BVT: {query_variant['BVT']}, QQ: {query_variant['QQ']}"
        )

        logger.info("Benchmarking decision tracking test PASSED!")


def main():
    """Run all tests."""
    logger.info("Starting enhanced _match_variants() tests...")

    try:
        test_basic_matching()
        test_multiallelic_matching()
        test_benchmarking_decisions()

        logger.info("All tests PASSED! ✅")
        logger.info("Enhanced _match_variants() implementation is working correctly.")

    except Exception as e:
        logger.error(f"Test failed: {e}")
        raise


if __name__ == "__main__":
    main()
