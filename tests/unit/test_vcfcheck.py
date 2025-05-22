#!/usr/bin/env python3
"""
Test cases for Python VCF checking implementation.

This module tests the functionality of python_vcfcheck.py to ensure
it correctly validates VCF files.
"""

import os
import tempfile
from pathlib import Path

import pytest
from Haplo.python_vcfcheck import VCFChecker


class TestVCFChecker:
    """Test cases for VCFChecker class."""

    @pytest.fixture
    def example_vcf(self):
        """Fixture to provide path to example VCF file."""
        # Use a VCF file from the example directory
        project_root = Path(__file__).parent.parent.parent
        example_dir = project_root / "example"
        vcf_path = example_dir / "hc.vcf.gz"

        if not vcf_path.exists():
            pytest.skip(f"Example VCF file not found: {vcf_path}")

        return str(vcf_path)

    @pytest.fixture
    def example_ref(self):
        """Fixture to provide path to example reference FASTA file."""
        project_root = Path(__file__).parent.parent.parent
        example_dir = project_root / "example"
        ref_path = example_dir / "chr21.fa"

        if not ref_path.exists():
            pytest.skip(f"Example reference file not found: {ref_path}")

        return str(ref_path)

    def test_init(self, example_ref):
        """Test initialization with default parameters."""
        checker = VCFChecker()
        assert checker.reference_path is None
        assert checker.strict is False
        assert checker.apply_filters is False

        # With reference
        checker = VCFChecker(reference_path=example_ref)
        assert checker.reference_path == example_ref
        assert checker.reference is not None

    def test_init_custom_params(self):
        """Test initialization with custom parameters."""
        checker = VCFChecker(strict=True, apply_filters=True)
        assert checker.strict is True
        assert checker.apply_filters is True

    def test_check_file(self, example_vcf):
        """Test checking a VCF file."""
        checker = VCFChecker()

        # Check file without output
        stats = checker.check_file(example_vcf)

        # Verify stats were collected
        assert "total_variants" in stats
        assert stats["total_variants"] > 0

        # Check other stats
        assert "filtered_variants" in stats
        assert "invalid_ref_alleles" in stats
        assert "missing_genotypes" in stats
        assert "invalid_genotypes" in stats
        assert "structural_variants" in stats
        assert "overlapping_variants" in stats

    def test_check_file_with_output(self, example_vcf):
        """Test checking a VCF file with output."""
        checker = VCFChecker()

        # Create temporary file for output
        with tempfile.NamedTemporaryFile(suffix=".txt", delete=False) as temp_file:
            temp_path = temp_file.name

        try:
            # Check file with output
            stats = checker.check_file(example_vcf, temp_path)

            # Verify output file was created
            assert os.path.exists(temp_path)

            # If the file has issues, the file should have content
            if (
                sum(
                    v
                    for k, v in stats.items()
                    if k != "total_variants" and k != "filtered_variants"
                )
                > 0
            ):
                assert os.path.getsize(temp_path) > 0

                # Check file format
                with open(temp_path) as f:
                    lines = f.readlines()
                    # Header line should be present
                    assert lines[0].startswith("CHROM")
        finally:
            # Clean up
            if os.path.exists(temp_path):
                os.unlink(temp_path)

    def test_check_with_reference(self, example_vcf, example_ref):
        """Test checking a VCF file with reference sequence."""
        checker = VCFChecker(reference_path=example_ref)

        # Check file with reference
        stats = checker.check_file(example_vcf)

        # Verify stats were collected
        assert "total_variants" in stats
        assert stats["total_variants"] > 0

        # Since we're checking reference alleles, this should be populated
        assert "invalid_ref_alleles" in stats

    def test_is_filtered(self):
        """Test the _is_filtered method."""

        # This test requires creating a mock VariantRecord
        class MockVariantRecord:
            def __init__(self, filter_values=None):
                self.filter = filter_values or []

        checker = VCFChecker()

        # Test unfiltered record
        unfiltered = MockVariantRecord([])
        assert checker._is_filtered(unfiltered) is False

        # Test PASS filter
        pass_filter = MockVariantRecord(["PASS"])
        assert checker._is_filtered(pass_filter) is False

        # Test filtered record
        filtered = MockVariantRecord(["LowQual"])
        assert checker._is_filtered(filtered) is True

        # Test multiple filters
        multi_filter = MockVariantRecord(["LowQual", "IndelGap"])
        assert checker._is_filtered(multi_filter) is True

    def test_check_header(self):
        """Test header checking functionality."""

        # This requires mocking a VCF header
        class MockHeader:
            def __init__(self, fields=None, samples=None, formats=None):
                self.fields = fields or {}
                self.samples = samples or []
                self.formats = formats or {}

            def __contains__(self, item):
                return item in self.fields

        class MockFormat:
            def __init__(self, name):
                self.name = name

        checker = VCFChecker()

        # Valid header
        valid_header = MockHeader(
            fields={"FILTER": True, "FORMAT": True, "INFO": True},
            samples=["Sample1"],
            formats={"GT": MockFormat("GT")},
        )
        issues = checker._check_header(valid_header)
        assert len(issues) == 0

        # Missing required fields
        missing_fields = MockHeader(
            fields={"FORMAT": True},
            samples=["Sample1"],
            formats={"GT": MockFormat("GT")},
        )
        issues = checker._check_header(missing_fields)
        assert len(issues) > 0
        assert any("Missing required header field" in issue for issue in issues)

        # No samples
        no_samples = MockHeader(
            fields={"FILTER": True, "FORMAT": True, "INFO": True},
            samples=[],
            formats={"GT": MockFormat("GT")},
        )
        issues = checker._check_header(no_samples)
        assert len(issues) > 0
        assert any("No sample columns" in issue for issue in issues)

        # Missing GT format
        no_gt = MockHeader(
            fields={"FILTER": True, "FORMAT": True, "INFO": True},
            samples=["Sample1"],
            formats={"DP": MockFormat("DP")},
        )
        issues = checker._check_header(no_gt)
        assert len(issues) > 0
        assert any("Missing GT format" in issue for issue in issues)

    def test_print_summary(self, example_vcf, capfd):
        """Test the print_summary method."""
        checker = VCFChecker()

        # Check file
        checker.check_file(example_vcf)

        # Print summary
        checker.print_summary()

        # Capture output
        out, err = capfd.readouterr()

        # Verify output contains key information
        assert "VCF Check Summary" in out
        assert "Total variants:" in out
        assert "Total issues:" in out


if __name__ == "__main__":
    pytest.main(["-xvs", __file__])
