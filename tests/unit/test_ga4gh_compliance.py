#!/usr/bin/env python3
"""
Tests for GA4GH compliance implementation.
"""

from unittest import mock
from unittest.mock import MagicMock, patch

import pandas as pd
import pytest
from pysam import VariantHeader, VariantRecord

from hap_py.haplo.ga4gh_compliance import (
    GA4GHDecision,
    GA4GHDecisionDetail,
    GA4GHFormatter,
    GA4GHMetrics,
    GA4GHStratification,
    GA4GHVariantType,
)


class TestGA4GHDecisionEnum:
    """Tests for GA4GHDecision enum."""

    def test_decision_values(self):
        """Test enum values match GA4GH specifications."""
        assert GA4GHDecision.TP.value == "TP"
        assert GA4GHDecision.FP.value == "FP"
        assert GA4GHDecision.FN.value == "FN"
        assert GA4GHDecision.N.value == "N"
        assert GA4GHDecision.UNK.value == "UNK"


class TestGA4GHFormatter:
    """Tests for GA4GHFormatter class."""

    @pytest.fixture
    def mock_vcf_header(self):
        """Create a mock VCF header."""
        header = MagicMock(spec=VariantHeader)
        header.add_meta = MagicMock()
        return header

    @pytest.fixture
    def mock_vcf_record(self):
        """Create a mock VCF record."""
        record = MagicMock(spec=VariantRecord)
        record.samples = [MagicMock()]
        record.samples[0].__getitem__ = MagicMock()
        record.samples[0].__setitem__ = MagicMock()
        record.info = {}
        record.ref = "A"
        record.alts = ["G"]
        record.header = MagicMock()
        record.header.info = {"PctSimilarity": MagicMock()}
        return record

    def test_init(self):
        """Test GA4GHFormatter initialization."""
        formatter = GA4GHFormatter()
        assert formatter.add_extra_fields is False

        formatter = GA4GHFormatter(add_extra_fields=True)
        assert formatter.add_extra_fields is True

    def test_format_vcf_header(self, mock_vcf_header):
        """Test adding GA4GH fields to VCF header."""
        formatter = GA4GHFormatter()
        result = formatter.format_vcf_header(mock_vcf_header)

        # Verify FORMAT fields were added
        for field in ["BD", "BK", "QD", "QK"]:
            mock_vcf_header.add_meta.assert_any_call(
                "FORMAT",
                items=[
                    ("ID", field),
                    ("Number", mock.ANY),
                    ("Type", mock.ANY),
                    ("Description", mock.ANY),
                ],
            )

        # Verify INFO fields were added
        for field in ["Regions", "TruthStatus", "QueryStatus", "Subtype"]:
            mock_vcf_header.add_meta.assert_any_call(
                "INFO",
                items=[
                    ("ID", field),
                    ("Number", mock.ANY),
                    ("Type", mock.ANY),
                    ("Description", mock.ANY),
                ],
            )

        # Verify META field for GA4GH standard was added
        mock_vcf_header.add_meta.assert_any_call(
            "META",
            items=[
                ("ID", "ga4gh_standard"),
                ("Version", "1.0"),
                ("Description", "This VCF follows GA4GH benchmarking standards"),
            ],
        )

        assert result == mock_vcf_header

    def test_format_vcf_header_with_extra_fields(self, mock_vcf_header):
        """Test adding extra fields to VCF header."""
        formatter = GA4GHFormatter(add_extra_fields=True)
        result = formatter.format_vcf_header(mock_vcf_header)

        # Verify extra INFO fields were added
        mock_vcf_header.add_meta.assert_any_call(
            "INFO",
            items=[
                ("ID", "PctSimilarity"),
                ("Number", "1"),
                ("Type", "Float"),
                ("Description", "Percent similarity to matched variant"),
            ],
        )
        mock_vcf_header.add_meta.assert_any_call(
            "INFO",
            items=[
                ("ID", "MatchId"),
                ("Number", "1"),
                ("Type", "String"),
                ("Description", "Identifier of the matched variant"),
            ],
        )

        assert result == mock_vcf_header

    def test_annotate_record(self, mock_vcf_record):
        """Test adding GA4GH annotations to VCF record."""
        formatter = GA4GHFormatter()
        result = formatter.annotate_record(
            mock_vcf_record,
            GA4GHDecision.TP,
            GA4GHDecision.TP,
            GA4GHDecisionDetail.GT_MATCH,
            GA4GHDecisionDetail.GT_MATCH,
            regions=["highconf", "exome"],
            variant_subtype=GA4GHVariantType.SNP,
        )

        # Check FORMAT fields
        mock_vcf_record.samples[0].__setitem__.assert_any_call("BD", "TP")
        mock_vcf_record.samples[0].__setitem__.assert_any_call("QD", "TP")
        mock_vcf_record.samples[0].__setitem__.assert_any_call("BK", "gt-match")
        mock_vcf_record.samples[0].__setitem__.assert_any_call("QK", "gt-match")

        # Check INFO fields
        assert result.info["TruthStatus"] == "TP"
        assert result.info["QueryStatus"] == "TP"
        assert result.info["Regions"] == ["highconf", "exome"]
        assert result.info["Subtype"] == "SNP"

    def test_annotate_record_auto_subtype(self, mock_vcf_record):
        """Test auto-detecting variant subtype."""
        formatter = GA4GHFormatter()

        # Test SNP
        mock_vcf_record.ref = "A"
        mock_vcf_record.alts = ["G"]
        result = formatter.annotate_record(
            mock_vcf_record, GA4GHDecision.TP, GA4GHDecision.TP
        )
        assert result.info["Subtype"] == "SNP"

        # Test INDEL
        mock_vcf_record.ref = "AG"
        mock_vcf_record.alts = ["A"]
        result = formatter.annotate_record(
            mock_vcf_record, GA4GHDecision.TP, GA4GHDecision.TP
        )
        assert result.info["Subtype"] == "INDEL"

        # Test COMPLEX
        mock_vcf_record.ref = "AT"
        mock_vcf_record.alts = ["GC"]
        result = formatter.annotate_record(
            mock_vcf_record, GA4GHDecision.TP, GA4GHDecision.TP
        )
        assert result.info["Subtype"] == "COMPLEX"


class TestGA4GHStratification:
    """Tests for GA4GHStratification class."""

    @pytest.fixture
    def mock_pybedtools_import(self):
        """Create mock for pybedtools import."""
        with patch.dict("sys.modules", {"pybedtools": MagicMock()}):
            import sys

            pybedtools = sys.modules["pybedtools"]
            pybedtools.BedTool = MagicMock()
            pybedtools.BedTool.return_value.intersect.return_value = [
                1
            ]  # Non-empty intersection
            yield pybedtools

    def test_init(self):
        """Test GA4GHStratification initialization."""
        stratification = GA4GHStratification()
        assert stratification.regions == {}

        stratification = GA4GHStratification({"highconf": "highconf.bed"})
        assert "highconf" in stratification.regions

    def test_add_region(self, mock_pybedtools_import):
        """Test adding stratification region."""
        stratification = GA4GHStratification()
        stratification.add_region("highconf", "highconf.bed")

        assert "highconf" in stratification.regions
        assert stratification.regions["highconf"]["bed_file"] == "highconf.bed"

    @patch("importlib.import_module")
    def test_add_region_no_pybedtools(self, mock_import):
        """Test adding region when pybedtools is not available."""
        mock_import.side_effect = ImportError("No module named 'pybedtools'")

        stratification = GA4GHStratification()
        stratification.add_region("highconf", "highconf.bed")

        assert "highconf" in stratification.regions
        assert stratification.regions["highconf"]["bed_file"] == "highconf.bed"
        assert stratification.regions["highconf"]["region"] is None

    @patch("importlib.import_module")
    def test_get_region_ids_for_variant(self, mock_import):
        """Test getting regions for a variant."""
        # Mock pybedtools module
        mock_pybedtools = MagicMock()
        mock_bedtool = MagicMock()
        mock_bedtool.intersect.return_value = [1]  # Non-empty intersection
        mock_pybedtools.BedTool.return_value = mock_bedtool
        mock_import.return_value = mock_pybedtools

        stratification = GA4GHStratification()
        stratification.regions = {
            "highconf": {"bed_file": "highconf.bed", "region": mock_bedtool},
            "exome": {"bed_file": "exome.bed", "region": mock_bedtool},
        }

        regions = stratification.get_region_ids_for_variant("chr1", 1000, 1001)
        assert sorted(regions) == sorted(["highconf", "exome"])

    def test_stratify_variants(self):
        """Test stratifying variants by region."""
        # Create a mock for get_region_ids_for_variant
        stratification = GA4GHStratification()
        stratification.get_region_ids_for_variant = MagicMock(return_value=["highconf"])
        stratification.regions = {
            "highconf": {"bed_file": "highconf.bed", "region": None}
        }

        # Create test DataFrame
        df = pd.DataFrame(
            {
                "chrom": ["chr1", "chr1", "chr2"],
                "pos": [1000, 2000, 3000],
                "end": [1001, 2001, 3001],
            }
        )

        result = stratification.stratify_variants(df)

        # Check results
        assert "all" in result
        assert len(result["all"]) == 3  # All variants
        assert "highconf" in result
        assert len(result["highconf"]) == 3  # All variants in highconf region


class TestGA4GHMetrics:
    """Tests for GA4GHMetrics class."""

    def test_init(self):
        """Test GA4GHMetrics initialization."""
        metrics = GA4GHMetrics()
        assert metrics.bootstrap_iterations == 1000
        assert metrics.ci_level == 0.95

        metrics = GA4GHMetrics(bootstrap_iterations=100, ci_level=0.99)
        assert metrics.bootstrap_iterations == 100
        assert metrics.ci_level == 0.99

    def test_calculate_precision(self):
        """Test precision calculation."""
        metrics = GA4GHMetrics()

        # Test basic calculation
        precision = metrics.calculate_precision(90, 10)
        assert precision == 0.9

        # Test zero denominator
        precision = metrics.calculate_precision(0, 0)
        assert precision == 0.0

    def test_calculate_precision_with_ci(self):
        """Test precision calculation with confidence intervals."""
        metrics = GA4GHMetrics()

        precision, lower, upper = metrics.calculate_precision(90, 10, with_ci=True)
        assert precision == 0.9
        assert 0 < lower < precision
        assert precision < upper <= 1.0

    def test_calculate_recall(self):
        """Test recall calculation."""
        metrics = GA4GHMetrics()

        # Test basic calculation
        recall = metrics.calculate_recall(90, 10)
        assert recall == 0.9

        # Test zero denominator
        recall = metrics.calculate_recall(0, 0)
        assert recall == 0.0

    def test_calculate_recall_with_ci(self):
        """Test recall calculation with confidence intervals."""
        metrics = GA4GHMetrics()

        recall, lower, upper = metrics.calculate_recall(90, 10, with_ci=True)
        assert recall == 0.9
        assert 0 < lower < recall
        assert recall < upper <= 1.0

    def test_calculate_f1(self):
        """Test F1 calculation."""
        metrics = GA4GHMetrics()

        # Test basic calculation
        f1 = metrics.calculate_f1(0.9, 0.9)
        assert f1 == 0.9

        # Test zero denominator
        f1 = metrics.calculate_f1(0, 0)
        assert f1 == 0.0

    def test_calculate_f1_with_ci(self):
        """Test F1 calculation with confidence intervals."""
        metrics = GA4GHMetrics()

        f1, lower, upper = metrics.calculate_f1(0.9, 0.9, with_ci=True)
        assert f1 == 0.9
        assert lower <= f1
        assert f1 <= upper

    def test_calculate_metrics(self):
        """Test calculating all metrics together."""
        metrics = GA4GHMetrics()

        result = metrics.calculate_metrics(90, 10, 10)
        assert result["precision"] == 0.9
        assert result["recall"] == 0.9
        assert result["f1"] == 0.9
        assert result["tp"] == 90
        assert result["fp"] == 10
        assert result["fn"] == 10

    def test_calculate_metrics_with_ci(self):
        """Test calculating all metrics with confidence intervals."""
        metrics = GA4GHMetrics()

        result = metrics.calculate_metrics(90, 10, 10, with_ci=True)

        # Check precision with CI
        precision, prec_lower, prec_upper = result["precision"]
        assert precision == 0.9
        assert 0 < prec_lower < precision
        assert precision < prec_upper <= 1.0

        # Check recall with CI
        recall, rec_lower, rec_upper = result["recall"]
        assert recall == 0.9
        assert 0 < rec_lower < recall
        assert recall < rec_upper <= 1.0

        # Check F1 with CI
        f1, f1_lower, f1_upper = result["f1"]
        assert f1 == 0.9
        assert f1_lower <= f1
        assert f1 <= f1_upper


from unittest.mock import patch


@pytest.fixture
def mock_vcf_record():
    """Create a mock VCF record for testing."""
    record = MagicMock()
    record.chrom = "chr1"
    record.pos = 12345
    record.samples = {"TRUTH": {}, "QUERY": {}}
    record.info = {}
    return record


@pytest.fixture
def mock_vcf_header():
    """Create a mock VCF header for testing."""
    header = MagicMock()
    header.formats = {}
    header.info = {}
    header.add_format = MagicMock()
    header.add_info = MagicMock()
    return header
