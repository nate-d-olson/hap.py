#!/usr/bin/env python3
"""
Test cases for Python quantify implementation.

This module tests the functionality of python_quantify.py to ensure
it correctly quantifies variant calls in VCF files.
"""

import json
import os
import tempfile
from pathlib import Path

import pandas as pd
import pytest
from Haplo.python_quantify import QuantifyEngine


class TestQuantifyEngine:
    """Test cases for QuantifyEngine class."""

    @pytest.fixture
    def example_vcfs(self):
        """Fixture to provide paths to example VCF files."""
        # Use VCF files from the example directory
        project_root = Path(__file__).parent.parent.parent
        example_dir = project_root / "example"

        truth_vcf = example_dir / "chr21.refcalls.vcf.gz"
        query_vcf = example_dir / "chr21.refcalls.vcf.gz"  # Using same file for testing

        if not truth_vcf.exists() or not query_vcf.exists():
            pytest.skip("Example VCF files not found")

        return {
            "truth": str(truth_vcf),
            "query": str(query_vcf),
        }

    @pytest.fixture
    def example_ref(self):
        """Fixture to provide path to example reference FASTA file."""
        project_root = Path(__file__).parent.parent.parent
        example_dir = project_root / "example"
        ref_path = example_dir / "chr21.fa"

        if not ref_path.exists():
            pytest.skip(f"Example reference file not found: {ref_path}")

        return str(ref_path)

    @pytest.fixture
    def temp_output_prefix(self):
        """Fixture to provide a temporary output prefix."""
        with tempfile.TemporaryDirectory() as temp_dir:
            output_prefix = os.path.join(temp_dir, "test_output")
            yield output_prefix

    def test_init(self, example_vcfs):
        """Test initialization with default parameters."""
        engine = QuantifyEngine(
            truth_vcf=example_vcfs["truth"], query_vcf=example_vcfs["query"]
        )

        assert engine.truth_vcf == example_vcfs["truth"]
        assert engine.query_vcf == example_vcfs["query"]
        assert engine.reference is None
        assert engine.regions is None
        assert engine.apply_filters is False
        assert engine.output_vtc is False

    def test_init_custom_params(self, example_vcfs, example_ref):
        """Test initialization with custom parameters."""
        engine = QuantifyEngine(
            truth_vcf=example_vcfs["truth"],
            query_vcf=example_vcfs["query"],
            reference=example_ref,
            apply_filters=True,
            output_vtc=True,
        )

        assert engine.truth_vcf == example_vcfs["truth"]
        assert engine.query_vcf == example_vcfs["query"]
        assert engine.reference == example_ref
        assert engine.regions is None
        assert engine.apply_filters is True
        assert engine.output_vtc is True

    def test_quantify(self, example_vcfs):
        """Test quantifying variants."""
        engine = QuantifyEngine(
            truth_vcf=example_vcfs["truth"], query_vcf=example_vcfs["query"]
        )

        results = engine.quantify()

        # Basic structure checks
        assert "metrics" in results
        assert "stratifications" in results

        # Metrics checks
        metrics = results["metrics"]
        assert "TP" in metrics
        assert "FP" in metrics
        assert "FN" in metrics
        assert "PRECISION" in metrics
        assert "RECALL" in metrics
        assert "F1" in metrics

        # When using same file for truth and query, we expect perfect results
        assert metrics["PRECISION"] == pytest.approx(1.0)
        assert metrics["RECALL"] == pytest.approx(1.0)
        assert metrics["F1"] == pytest.approx(1.0)

        # Stratifications checks
        strats = results["stratifications"]
        assert "variant_type" in strats
        assert "indel_size" in strats
        assert "zygosity" in strats

    def test_write_results(self, example_vcfs, temp_output_prefix):
        """Test writing results to files."""
        engine = QuantifyEngine(
            truth_vcf=example_vcfs["truth"],
            query_vcf=example_vcfs["query"],
            output_vtc=True,
        )

        engine.quantify()
        engine.write_results(temp_output_prefix)

        # Check output files
        metrics_file = f"{temp_output_prefix}.metrics.json"
        summary_file = f"{temp_output_prefix}.summary.tsv"
        truth_vtc_file = f"{temp_output_prefix}.truth.vtc.tsv"
        query_vtc_file = f"{temp_output_prefix}.query.vtc.tsv"

        assert os.path.exists(metrics_file)
        assert os.path.exists(summary_file)
        assert os.path.exists(truth_vtc_file)
        assert os.path.exists(query_vtc_file)

        # Check JSON content
        with open(metrics_file) as f:
            metrics_data = json.load(f)
            assert "metrics" in metrics_data
            assert "stratifications" in metrics_data

        # Check TSV content
        summary_df = pd.read_csv(summary_file, sep="\t")
        assert "Type" in summary_df.columns
        assert "TP" in summary_df.columns
        assert "Precision" in summary_df.columns

        # Check VTC content
        truth_vtc_df = pd.read_csv(truth_vtc_file, sep="\t")
        query_vtc_df = pd.read_csv(query_vtc_file, sep="\t")

        assert "category" in truth_vtc_df.columns
        assert "category" in query_vtc_df.columns

    def test_is_filtered(self):
        """Test the _is_filtered method."""

        # This test requires creating a mock VariantRecord
        class MockVariantRecord:
            def __init__(self, filter_values=None):
                self.filter = filter_values or []

        engine = QuantifyEngine(truth_vcf="dummy.vcf", query_vcf="dummy.vcf")

        # Test unfiltered record
        unfiltered = MockVariantRecord([])
        assert engine._is_filtered(unfiltered) is False

        # Test PASS filter
        pass_filter = MockVariantRecord(["PASS"])
        assert engine._is_filtered(pass_filter) is False

        # Test filtered record
        filtered = MockVariantRecord(["LowQual"])
        assert engine._is_filtered(filtered) is True

        # Test multiple filters
        multi_filter = MockVariantRecord(["LowQual", "IndelGap"])
        assert engine._is_filtered(multi_filter) is True

    def test_get_variant_type(self):
        """Test the _get_variant_type method."""

        # This test requires creating a mock VariantRecord
        class MockVariant:
            def __init__(self, ref, alts=None):
                self.ref = ref
                self.alts = alts

        engine = QuantifyEngine(truth_vcf="dummy.vcf", query_vcf="dummy.vcf")

        # Test SNP
        snp = MockVariant("A", ["G"])
        assert engine._get_variant_type(snp) == "SNP"

        # Test MNP
        mnp = MockVariant("AT", ["GC"])
        assert engine._get_variant_type(mnp) == "MNP"

        # Test insertion
        ins = MockVariant("A", ["ACGT"])
        assert engine._get_variant_type(ins) == "INS"

        # Test deletion
        deletion = MockVariant("ACGT", ["A"])
        assert engine._get_variant_type(deletion) == "DEL"

        # Test reference
        ref = MockVariant("A", None)
        assert engine._get_variant_type(ref) == "REF"

    def test_match_variants(self):
        """Test variant matching functionality."""
        # Setup test variants
        truth_variants = [
            {"chrom": "chr1", "pos": 100, "ref": "A", "alt": "G", "source": "truth"},
            {"chrom": "chr1", "pos": 200, "ref": "AT", "alt": "A", "source": "truth"},
            {"chrom": "chr1", "pos": 300, "ref": "G", "alt": "T", "source": "truth"},
        ]

        query_variants = [
            {"chrom": "chr1", "pos": 100, "ref": "A", "alt": "G", "source": "query"},
            {"chrom": "chr1", "pos": 250, "ref": "C", "alt": "T", "source": "query"},
            {"chrom": "chr1", "pos": 300, "ref": "G", "alt": "C", "source": "query"},
        ]

        engine = QuantifyEngine(truth_vcf="dummy.vcf", query_vcf="dummy.vcf")

        # Set up variants
        engine.truth_variants = truth_variants
        engine.query_variants = query_variants

        # Match variants
        engine._match_variants()

        # Check matching
        assert engine.truth_variants[0]["match"] is True
        assert engine.truth_variants[1]["match"] is False
        assert engine.truth_variants[2]["match"] is False

        assert engine.query_variants[0]["match"] is True
        assert engine.query_variants[1]["match"] is False
        assert engine.query_variants[2]["match"] is False


if __name__ == "__main__":
    pytest.main(["-xvs", __file__])
