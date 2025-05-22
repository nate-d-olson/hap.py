"""
Integration tests for quantify stratification functionality.
Migrated from src/sh/run_quantify_stratification_test.sh
"""

import subprocess
import sys
from pathlib import Path

import pytest

from tests.utils import get_project_root


@pytest.mark.integration
def test_quantify_stratification(tmp_path):
    """Test stratified counting with quantify functionality in hap.py."""
    # Get paths to required files
    project_root = get_project_root()
    example_dir = project_root / "example" / "happy"
    hap_py_script = project_root / "src" / "python" / "hap.py"
    compare_summaries_script = project_root / "src" / "sh" / "compare_summaries.py"
    compare_extended_script = project_root / "src" / "sh" / "compare_extended.py"

    # Input files
    reference = example_dir / "hg38.chr21.fa"
    truth_vcf = example_dir / "PG_NA12878_hg38-chr21.vcf.gz"
    query_vcf = example_dir / "NA12878-GATK3-chr21.vcf.gz"
    confident_regions = example_dir / "PG_Conf_hg38-chr21.bed.gz"
    stratification_file = example_dir / "stratification.tsv"
    expected_summary = example_dir / "expected-stratified.summary.csv"
    expected_extended = example_dir / "expected-stratified.extended.csv"

    # Check that required files exist
    assert truth_vcf.exists(), f"Truth VCF {truth_vcf} not found"
    assert query_vcf.exists(), f"Query VCF {query_vcf} not found"
    assert reference.exists(), f"Reference {reference} not found"
    assert (
        confident_regions.exists()
    ), f"Confident regions file {confident_regions} not found"
    assert (
        stratification_file.exists()
    ), f"Stratification file {stratification_file} not found"
    assert expected_summary.exists(), f"Expected summary {expected_summary} not found"
    assert (
        expected_extended.exists()
    ), f"Expected extended {expected_extended} not found"

    # Output file paths
    output_prefix = tmp_path / "stratified_test"
    output_summary = Path(str(output_prefix) + ".summary.csv")
    output_extended = Path(str(output_prefix) + ".extended.csv")

    # Run hap.py with stratification
    cmd = [
        sys.executable,
        str(hap_py_script),
        "-l",
        "chr21",
        str(truth_vcf),
        str(query_vcf),
        "-f",
        str(confident_regions),
        "-r",
        str(reference),
        "-o",
        str(output_prefix),
        "--stratification",
        str(stratification_file),
        "--force-interactive",
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0, f"hap.py failed with stratification: {result.stderr}"

    # Compare summary files
    compare_summary_cmd = [
        sys.executable,
        str(compare_summaries_script),
        str(output_summary),
        str(expected_summary),
    ]

    compare_summary_result = subprocess.run(
        compare_summary_cmd, capture_output=True, text=True
    )
    assert compare_summary_result.returncode == 0, (
        f"Summary comparison failed: {compare_summary_result.stderr}\n"
        f"Check diff {output_summary} {expected_summary}"
    )

    # Compare extended files
    compare_extended_cmd = [
        sys.executable,
        str(compare_extended_script),
        str(output_extended),
        str(expected_extended),
    ]

    compare_extended_result = subprocess.run(
        compare_extended_cmd, capture_output=True, text=True
    )
    assert compare_extended_result.returncode == 0, (
        f"Extended CSV comparison failed: {compare_extended_result.stderr}\n"
        f"Check diff {output_extended} {expected_extended}"
    )
