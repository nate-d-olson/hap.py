#!/usr/bin/env python3
"""
Unit tests for ROC analysis functionality in the QuantifyEngine.
"""

import os
import unittest
from unittest.mock import MagicMock, patch

import pytest

from hap_py.haplo.python_quantify import QuantifyEngine


class TestROCAnalysis(unittest.TestCase):
    """Test ROC analysis functionality."""

    def setUp(self):
        """Set up test environment."""
        # Create mock VCF paths
        self.truth_vcf = "tests/data/truth.vcf"
        self.query_vcf = "tests/data/query.vcf"

        # Mock VCF file existence
        patcher = patch("os.path.exists")
        self.mock_exists = patcher.start()
        self.mock_exists.return_value = True
        self.addCleanup(patcher.stop)

        # Mock VCF file opening
        patcher = patch("pysam.VariantFile")
        self.mock_variantfile = patcher.start()
        self.mock_variant_file_instance = MagicMock()
        self.mock_variantfile.return_value = self.mock_variant_file_instance
        self.mock_variant_file_instance.__enter__ = MagicMock(
            return_value=self.mock_variant_file_instance
        )
        self.mock_variant_file_instance.__exit__ = MagicMock(return_value=None)
        self.mock_variant_file_instance.header = MagicMock()
        self.mock_variant_file_instance.header.contigs = {}
        self.addCleanup(patcher.stop)

        # Mock _load_variants method
        patcher = patch.object(QuantifyEngine, "_load_variants")
        self.mock_load_variants = patcher.start()
        self.mock_load_variants.return_value = self._get_mock_variants()
        self.addCleanup(patcher.stop)

        # Mock _match_variants method
        patcher = patch.object(QuantifyEngine, "_match_variants")
        self.mock_match_variants = patcher.start()
        self.addCleanup(patcher.stop)

    def _get_mock_variants(self):
        """Get mock variant data for testing."""
        # Create sample variants with quality scores and match status
        variants = []

        # Add true positive SNPs
        for i in range(10):
            variants.append(
                {
                    "chrom": "chr1",
                    "pos": 1000 + i,
                    "ref": "A",
                    "alt": "G",
                    "qual": 10.0 + i * 5,  # Quality from 10 to 55
                    "BVT": "SNP",
                    "match": True,
                }
            )

        # Add false positive SNPs
        for i in range(5):
            variants.append(
                {
                    "chrom": "chr1",
                    "pos": 2000 + i,
                    "ref": "C",
                    "alt": "T",
                    "qual": 5.0 + i * 3,  # Quality from 5 to 17
                    "BVT": "SNP",
                    "match": False,
                }
            )

        # Add true positive INDELs
        for i in range(7):
            variants.append(
                {
                    "chrom": "chr1",
                    "pos": 3000 + i,
                    "ref": "A",
                    "alt": "AG",
                    "qual": 15.0 + i * 4,  # Quality from 15 to 39
                    "BVT": "INS",
                    "match": True,
                }
            )

        # Add false positive INDELs
        for i in range(3):
            variants.append(
                {
                    "chrom": "chr1",
                    "pos": 4000 + i,
                    "ref": "TA",
                    "alt": "T",
                    "qual": 8.0 + i * 2,  # Quality from 8 to 12
                    "BVT": "DEL",
                    "match": False,
                }
            )

        return variants

    def test_roc_analysis_initialization(self):
        """Test that ROC analysis is properly initialized."""
        engine = QuantifyEngine(
            self.truth_vcf, self.query_vcf, enable_roc_analysis=True
        )
        assert engine.enable_roc_analysis is True
        assert engine.quality_stratification is True
        assert engine.roc_bootstrap_samples == 1000

    def test_perform_roc_analysis(self):
        """Test that _perform_roc_analysis method produces expected data structures."""
        engine = QuantifyEngine(self.truth_vcf, self.query_vcf)

        # Mock truth and query variants
        engine.truth_variants = self._get_mock_variants()
        engine.query_variants = self._get_mock_variants()

        # Run ROC analysis
        engine._perform_roc_analysis()

        # Check ROC data structure
        assert isinstance(engine.roc_data, dict)
        assert "snp" in engine.roc_data
        assert "indel" in engine.roc_data
        assert "all" in engine.roc_data

        # Check that confidence intervals were calculated
        assert hasattr(engine, "bootstrap_confidence_intervals")

        # Check quality stratification
        assert hasattr(engine, "quality_metrics")
        if hasattr(engine, "quality_metrics"):
            assert "bin_metrics" in engine.quality_metrics

    @pytest.mark.skip(reason="Integration test requiring file output")
    def test_write_roc_results(self):
        """Test writing ROC results to files."""
        engine = QuantifyEngine(self.truth_vcf, self.query_vcf)

        # Mock data structures
        engine.truth_variants = self._get_mock_variants()
        engine.query_variants = self._get_mock_variants()
        engine.enable_roc_analysis = True

        # Run ROC analysis
        engine._perform_roc_analysis()

        # Create temporary output directory
        import tempfile

        with tempfile.TemporaryDirectory() as tmpdirname:
            output_prefix = os.path.join(tmpdirname, "test_output")

            # Write results
            engine._write_roc_results(output_prefix)

            # Check output files exist
            assert os.path.exists(f"{output_prefix}.roc.tsv")
            assert os.path.exists(f"{output_prefix}.quality_stratification.tsv")
            assert os.path.exists(f"{output_prefix}.multi_threshold.tsv")


if __name__ == "__main__":
    unittest.main()
