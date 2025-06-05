#!/usr/bin/env python3
"""
Integration tests for GA4GH compliance functionality.

This module tests the integration of GA4GH standards with the
QuantifyEngine and ensures that GA4GH-compliant outputs are generated.
"""

import os
from pathlib import Path
from unittest.mock import MagicMock, patch

import pandas as pd
import pysam
import pytest

from hap_py.haplo.ga4gh_compliance import (
    GA4GHDecision,
    GA4GHDecisionDetail,
    GA4GHFormatter,
    GA4GHMetrics,
    GA4GHStratification,
)
from hap_py.haplo.ga4gh_integration import (
    GA4GHIntegration,
    enhance_quantify_engine_with_ga4gh,
)
from hap_py.haplo.python_quantify import QuantifyEngine


class TestGA4GHIntegration:
    """Tests for GA4GH integration with QuantifyEngine."""

    @pytest.fixture
    def mock_quantify_engine(self):
        """Create a mock QuantifyEngine for testing."""
        engine = MagicMock(spec=QuantifyEngine)
        engine.quantify_method = "ga4gh"
        engine.stratification_regions = {"highconf": "highconf.bed"}
        engine.confident_regions = "confident.bed"
        engine.results_by_region = {
            "all": {"TP": 90, "FP": 10, "FN": 10},
            "highconf": {"TP": 85, "FP": 5, "FN": 5},
        }
        return engine

    @pytest.fixture
    def mock_vcf_header(self):
        """Create a mock VCF header."""
        header = MagicMock(spec=pysam.VariantHeader)
        header.add_meta = MagicMock()
        return header

    @pytest.fixture
    def mock_vcf_record(self):
        """Create a mock VCF record."""
        record = MagicMock(spec=pysam.VariantRecord)
        record.samples = [MagicMock()]
        record.samples[0].__getitem__ = MagicMock()
        record.samples[0].__setitem__ = MagicMock()
        record.info = {}
        record.ref = "A"
        record.alts = ["G"]
        record.header = MagicMock()
        record.header.info = {"PctSimilarity": MagicMock()}
        return record

    @pytest.fixture
    def ga4gh_integration(self):
        """Create a GA4GHIntegration instance for testing."""
        return GA4GHIntegration(
            stratification_beds={"highconf": "highconf.bed"},
            confidence_regions="confident.bed",
        )

    def test_init(self):
        """Test GA4GHIntegration initialization."""
        integration = GA4GHIntegration()
        assert isinstance(integration.formatter, GA4GHFormatter)
        assert isinstance(integration.stratification, GA4GHStratification)
        assert isinstance(integration.metrics, GA4GHMetrics)

    def test_prepare_vcf_header(self, ga4gh_integration, mock_vcf_header):
        """Test preparing a GA4GH-compliant VCF header."""
        with patch.object(
            ga4gh_integration.formatter, "format_vcf_header"
        ) as mock_format:
            ga4gh_integration.prepare_vcf_header(mock_vcf_header)
            mock_format.assert_called_once_with(mock_vcf_header)

    def test_transform_match_to_ga4gh(self, ga4gh_integration):
        """Test transforming a match tuple to GA4GH decisions."""
        # Create test DataFrames
        truth_df = pd.DataFrame(
            {
                "chrom": ["chr1"],
                "pos": [1000],
                "ref": ["A"],
                "alt": ["G"],
            }
        )
        query_df = pd.DataFrame(
            {
                "chrom": ["chr1"],
                "pos": [1000],
                "ref": ["A"],
                "alt": ["G"],
            }
        )

        # Test true positive with genotype match
        truth_idx, query_idx, match_type = 0, 0, "gt-match"
        decisions = ga4gh_integration.transform_match_to_ga4gh(
            (truth_idx, query_idx, match_type), truth_df, query_df
        )
        assert decisions[0] == GA4GHDecision.TP
        assert decisions[1] == GA4GHDecision.TP
        assert decisions[2] == GA4GHDecisionDetail.GT_MATCH
        assert decisions[3] == GA4GHDecisionDetail.GT_MATCH

        # Test false negative
        truth_idx, query_idx, match_type = 0, -1, "no-match"
        decisions = ga4gh_integration.transform_match_to_ga4gh(
            (truth_idx, query_idx, match_type), truth_df, query_df
        )
        assert decisions[0] == GA4GHDecision.FN
        assert decisions[1] == GA4GHDecision.UNK

        # Test false positive
        truth_idx, query_idx, match_type = -1, 0, "no-match"
        decisions = ga4gh_integration.transform_match_to_ga4gh(
            (truth_idx, query_idx, match_type), truth_df, query_df
        )
        assert decisions[0] == GA4GHDecision.UNK
        assert decisions[1] == GA4GHDecision.FP

    def test_annotate_vcf_record(self, ga4gh_integration, mock_vcf_record):
        """Test annotating a VCF record with GA4GH fields."""
        with patch.object(
            ga4gh_integration.formatter, "annotate_record"
        ) as mock_annotate:
            with patch.object(
                ga4gh_integration.stratification, "get_region_ids_for_variant"
            ) as mock_get_regions:
                mock_get_regions.return_value = ["highconf"]

                ga4gh_integration.annotate_vcf_record(
                    mock_vcf_record,
                    GA4GHDecision.TP,
                    GA4GHDecision.TP,
                    GA4GHDecisionDetail.GT_MATCH,
                    GA4GHDecisionDetail.GT_MATCH,
                    "chr1",
                    1000,
                )

                mock_get_regions.assert_called_once_with("chr1", 1000, None)
                mock_annotate.assert_called_once()

    def test_create_ga4gh_metrics(self, ga4gh_integration):
        """Test creating GA4GH metrics."""
        with patch.object(
            ga4gh_integration.metrics, "calculate_metrics"
        ) as mock_calculate:
            mock_calculate.return_value = {
                "precision": (0.9, 0.85, 0.95),
                "recall": (0.9, 0.85, 0.95),
                "f1": (0.9, 0.85, 0.95),
                "tp": 90,
                "fp": 10,
                "fn": 10,
            }

            metrics = ga4gh_integration.create_ga4gh_metrics(
                90, 10, 10, "highconf", True
            )

            mock_calculate.assert_called_once_with(90, 10, 10, with_ci=True)
            assert metrics["region"] == "highconf"
            assert metrics["precision"] == (0.9, 0.85, 0.95)

    def test_write_ga4gh_metrics_file(self, ga4gh_integration, tmp_path):
        """Test writing GA4GH metrics to a file."""
        metrics_by_region = {
            "all": {
                "precision": (0.9, 0.85, 0.95),
                "recall": (0.9, 0.85, 0.95),
                "f1": (0.9, 0.85, 0.95),
                "tp": 90,
                "fp": 10,
                "fn": 10,
                "types": {
                    "SNP": {
                        "precision": (0.95, 0.9, 0.99),
                        "recall": (0.95, 0.9, 0.99),
                        "f1": (0.95, 0.9, 0.99),
                        "tp": 80,
                        "fp": 5,
                        "fn": 5,
                    },
                    "INDEL": {
                        "precision": (0.8, 0.75, 0.85),
                        "recall": (0.8, 0.75, 0.85),
                        "f1": (0.8, 0.75, 0.85),
                        "tp": 10,
                        "fp": 5,
                        "fn": 5,
                    },
                },
            },
            "highconf": {
                "precision": (0.95, 0.9, 0.99),
                "recall": (0.95, 0.9, 0.99),
                "f1": (0.95, 0.9, 0.99),
                "tp": 85,
                "fp": 5,
                "fn": 5,
            },
        }

        output_path = tmp_path / "metrics.tsv"
        ga4gh_integration.write_ga4gh_metrics_file(metrics_by_region, output_path)

        assert output_path.exists()
        content = output_path.read_text()
        assert "Region\tType\tTP\tFP\tFN\tPrecision\tPrecision_lower\t" in content
        assert (
            "all\tALL\t90\t10\t10\t0.9\t0.85\t0.95\t0.9\t0.85\t0.95\t0.9\t0.85\t0.95"
            in content.replace("\n", "")
        )

    def test_enhance_quantify_engine(self, mock_quantify_engine):
        """Test enhancing a QuantifyEngine with GA4GH support."""
        enhanced = enhance_quantify_engine_with_ga4gh(mock_quantify_engine)

        assert hasattr(enhanced, "ga4gh")
        assert isinstance(enhanced.ga4gh, GA4GHIntegration)


@pytest.mark.integration
class TestGA4GHQuantificationIntegration:
    """Integration tests for GA4GH quantification."""

    def test_ga4gh_quantify_method(self, tmp_path):
        """Test using GA4GH quantification method."""
        # Skip if no test data available
        if not os.path.exists("example/integration/test1-truth.vcf"):
            pytest.skip("Test data not available")

        # Create a QuantifyEngine with GA4GH method
        engine = QuantifyEngine(
            truth_vcf="example/integration/test1-truth.vcf",
            query_vcf="example/integration/test1-query.vcf",
            output_prefix=str(tmp_path / "output"),
            quantify_method="ga4gh",
        )

        # Enhance the engine with GA4GH support
        enhanced_engine = enhance_quantify_engine_with_ga4gh(engine)

        # Run quantification
        results = enhanced_engine.quantify()

        # Check results
        assert results is not None
        assert "all" in results
        assert results["all"]["TP"] >= 0

        # Check for GA4GH-specific output files
        metrics_file = Path(str(tmp_path / "output") + ".ga4gh.metrics.tsv")
        assert metrics_file.exists() or not os.path.exists(
            "example/integration/test1-truth.vcf"
        ), "GA4GH metrics file should be created if test data exists"
