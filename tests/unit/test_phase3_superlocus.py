#!/usr/bin/env python3
"""
Unit tests for Phase 3 quantify implementation (Superlocus Analysis).

This module tests the Phase 3 functionality including:
- Superlocus analysis
- Region-based quantification
- Multi-sample analysis
"""

import os
import tempfile

import pytest

from hap_py.haplo.python_quantify import QuantifyEngine
from hap_py.haplo.quantify_phase3 import (
    MultiSampleQuantifier,
    RegionBasedQuantifier,
)


class TestPhase3BasicFunctionality:
    """Test basic Phase 3 component initialization and functionality."""

    def test_region_based_quantifier_initialization(self):
        """Test that RegionBasedQuantifier can be initialized properly."""
        region_quantifier = RegionBasedQuantifier()
        assert region_quantifier is not None
        assert hasattr(region_quantifier, "bed_regions")
        assert hasattr(region_quantifier, "load_bed_regions")

    def test_multi_sample_quantifier_initialization(self):
        """Test that MultiSampleQuantifier can be initialized properly."""
        multi_quantifier = MultiSampleQuantifier()
        assert multi_quantifier is not None
        assert hasattr(multi_quantifier, "samples")
        assert hasattr(multi_quantifier, "load_vcf_samples")

    def test_quantify_engine_phase3_parameters(self):
        """Test that QuantifyEngine accepts Phase 3 parameters."""
        # Create minimal test VCF files
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".vcf", delete=False
        ) as truth_vcf, tempfile.NamedTemporaryFile(
            mode="w", suffix=".vcf", delete=False
        ) as query_vcf:

            # Write minimal VCF content
            vcf_content = """##fileformat=VCFv4.2
##contig=<ID=chr21>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr21	100000	.	A	G	60	PASS	.
"""
            truth_vcf.write(vcf_content)
            query_vcf.write(vcf_content)
            truth_vcf.flush()
            query_vcf.flush()

            try:
                # Test QuantifyEngine with Phase 3 parameters
                engine = QuantifyEngine(
                    truth_vcf=truth_vcf.name,
                    query_vcf=query_vcf.name,
                    quantify_method="xcmp",
                    enable_superlocus_analysis=True,
                    enable_region_stratification=True,
                    enable_multi_sample=True,
                    superlocus_window=1000,
                )

                assert engine is not None
                assert hasattr(engine, "enable_superlocus_analysis")
                assert hasattr(engine, "enable_region_stratification")
                assert hasattr(engine, "enable_multi_sample")

            finally:
                # Clean up temporary files
                os.unlink(truth_vcf.name)
                os.unlink(query_vcf.name)


class TestRegionBasedQuantifier:
    """Test RegionBasedQuantifier functionality."""

    def create_test_bed_file(self):
        """Create a test BED file for region testing."""
        bed_content = """chr21\t26960070\t27230000\thigh_confidence_region
chr21\t27230000\t27590000\tmedium_confidence_region
chr21\t27590000\t28100000\tlow_confidence_region"""

        bed_file = tempfile.NamedTemporaryFile(mode="w", suffix=".bed", delete=False)
        bed_file.write(bed_content)
        bed_file.close()
        return bed_file.name

    def test_bed_region_loading(self):
        """Test loading BED regions from file."""
        bed_file = self.create_test_bed_file()

        try:
            region_quantifier = RegionBasedQuantifier()
            region_quantifier.load_bed_regions({"test_regions": bed_file})

            assert "test_regions" in region_quantifier.bed_regions
            assert len(region_quantifier.bed_regions["test_regions"]) > 0

        finally:
            os.unlink(bed_file)

    def test_variant_region_assignment(self):
        """Test assignment of variants to regions."""
        bed_file = self.create_test_bed_file()

        try:
            region_quantifier = RegionBasedQuantifier()
            region_quantifier.load_bed_regions({"test_regions": bed_file})

            # Test variants within and outside regions
            test_variants = [
                {
                    "chrom": "chr21",
                    "pos": 27000000,
                    "ref": "A",
                    "alt": "G",
                },  # Within region
                {
                    "chrom": "chr21",
                    "pos": 30000000,
                    "ref": "C",
                    "alt": "T",
                },  # Outside region
            ]

            stratified = region_quantifier.stratify_variants(test_variants)
            assert stratified is not None

        finally:
            os.unlink(bed_file)


class TestMultiSampleQuantifier:
    """Test MultiSampleQuantifier functionality."""

    def create_test_vcf(self, content_suffix=""):
        """Create a test VCF file."""
        vcf_content = f"""##fileformat=VCFv4.2
##contig=<ID=chr21>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr21	100000	.	A	G	60	PASS	.
chr21	200000	.	C	T	80	PASS	.{content_suffix}
"""
        vcf_file = tempfile.NamedTemporaryFile(mode="w", suffix=".vcf", delete=False)
        vcf_file.write(vcf_content)
        vcf_file.close()
        return vcf_file.name

    def test_multi_sample_loading(self):
        """Test loading multiple VCF samples."""
        vcf1 = self.create_test_vcf()
        vcf2 = self.create_test_vcf("\nchr21\t300000\t.\tG\tA\t70\tPASS\t.")

        try:
            multi_quantifier = MultiSampleQuantifier()
            multi_quantifier.load_vcf_samples([vcf1, vcf2])

            assert len(multi_quantifier.samples) == 2

        finally:
            os.unlink(vcf1)
            os.unlink(vcf2)

    def test_sample_comparison(self):
        """Test comparison between multiple samples."""
        vcf1 = self.create_test_vcf()
        vcf2 = self.create_test_vcf()

        try:
            multi_quantifier = MultiSampleQuantifier()
            multi_quantifier.load_vcf_samples([vcf1, vcf2])

            # Test sample comparison functionality
            comparison_results = multi_quantifier.compare_samples()
            assert comparison_results is not None

        finally:
            os.unlink(vcf1)
            os.unlink(vcf2)


class TestPhase3Integration:
    """Test integration of Phase 3 components."""

    def test_end_to_end_phase3_workflow(self):
        """Test complete Phase 3 workflow with all components."""
        # Create test files
        bed_content = """chr21\t26960070\t27230000\thigh_confidence"""
        vcf_content = """##fileformat=VCFv4.2
##contig=<ID=chr21>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr21	27000000	.	A	G	60	PASS	.
"""

        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".bed", delete=False
        ) as bed_file, tempfile.NamedTemporaryFile(
            mode="w", suffix=".vcf", delete=False
        ) as truth_vcf, tempfile.NamedTemporaryFile(
            mode="w", suffix=".vcf", delete=False
        ) as query_vcf:

            bed_file.write(bed_content)
            truth_vcf.write(vcf_content)
            query_vcf.write(vcf_content)
            bed_file.flush()
            truth_vcf.flush()
            query_vcf.flush()

            try:
                # Test integrated Phase 3 workflow
                engine = QuantifyEngine(
                    truth_vcf=truth_vcf.name,
                    query_vcf=query_vcf.name,
                    quantify_method="xcmp",
                    enable_superlocus_analysis=True,
                    enable_region_stratification=True,
                    region_bed_files={"high_confidence": bed_file.name},
                )

                # Verify engine configuration
                assert engine.enable_superlocus_analysis is True
                assert engine.enable_region_stratification is True
                assert "high_confidence" in engine.region_bed_files

                # Test that analysis can be initiated (without running full analysis)
                assert hasattr(engine, "run")

            finally:
                os.unlink(bed_file.name)
                os.unlink(truth_vcf.name)
                os.unlink(query_vcf.name)


@pytest.mark.integration
class TestPhase3IntegrationWithExampleData:
    """Integration tests using example data files."""

    def test_with_example_data(self):
        """Test Phase 3 functionality with actual example data if available."""
        # Check if example data is available
        truth_vcf = "example/integration/integrationtest.vcf"
        query_vcf = "example/integration/integrationtest_rhs.vcf"
        bed_file = "example/hc.bed"

        if not all(os.path.exists(f) for f in [truth_vcf, query_vcf, bed_file]):
            pytest.skip("Example data files not available")

        try:
            # Test with actual example data
            engine = QuantifyEngine(
                truth_vcf=truth_vcf,
                query_vcf=query_vcf,
                quantify_method="xcmp",
                enable_superlocus_analysis=True,
                enable_region_stratification=True,
                region_bed_files={"high_confidence": bed_file},
            )

            assert engine is not None

            # Test region loading
            region_quantifier = RegionBasedQuantifier()
            region_quantifier.load_bed_regions({"high_confidence": bed_file})
            assert len(region_quantifier.bed_regions["high_confidence"]) > 0

        except Exception as e:
            pytest.fail(f"Phase 3 integration test failed: {e}")
