"""
Integration tests for chromosome prefix detection.
Migrated from src/sh/run_chrprefix_test.sh
"""

import filecmp
import gzip
import subprocess
import sys
from pathlib import Path

import pytest

from tests.utils import get_project_root


@pytest.mark.integration
def test_numeric_chrs(tmp_path):
    """Test chr prefix detection with numeric chromosomes."""
    # Get paths to required files
    project_root = get_project_root()
    src_data_dir = project_root / "src" / "data" / "numeric_chrs"
    # Using CLI commands instead of script paths

    # Input files
    truth_vcf = src_data_dir / "truth.vcf"
    query_vcf = src_data_dir / "query.vcf"
    fp_bed = src_data_dir / "fp.bed"
    reference = src_data_dir / "test.fa"
    expected_vcf = src_data_dir / "expected.vcf"
    expected_summary = src_data_dir / "expected.summary.csv"

    # Check that required files exist
    assert truth_vcf.exists(), f"Truth VCF {truth_vcf} not found"
    assert query_vcf.exists(), f"Query VCF {query_vcf} not found"
    assert fp_bed.exists(), f"FP BED {fp_bed} not found"
    assert reference.exists(), f"Reference {reference} not found"
    assert expected_vcf.exists(), f"Expected VCF {expected_vcf} not found"
    assert expected_summary.exists(), f"Expected summary {expected_summary} not found"

    # Output file paths
    output_prefix = tmp_path / "numeric_test"
    output_vcf_gz = output_prefix.with_suffix(".vcf.gz")
    output_vcf = output_prefix.with_suffix(".vcf")
    output_summary = Path(str(output_prefix) + ".summary.csv")

    # Run hap.py on numeric chromosome files using CLI command
    cmd = [
        "hap",
        str(truth_vcf),
        str(query_vcf),
        "-f",
        str(fp_bed),
        "-o",
        str(output_prefix),
        "-X",
        "--reference",
        str(reference),
        "-V",
        "--force-interactive",
    ]

    result = subprocess.run(cmd, capture_output=True)
    assert (
        result.returncode == 0
    ), f"hap.py failed with numeric chromosomes: {result.stderr.decode()}"

    # Compare summary files
    compare_cmd = [
        sys.executable,
        str(compare_script),
        str(output_summary),
        str(expected_summary),
    ]

    compare_result = subprocess.run(compare_cmd, capture_output=True)
    assert compare_result.returncode == 0, "Summary comparison failed"

    # Compare VCF files
    with gzip.open(output_vcf_gz, "rt") as f_gz:
        vcf_content = [line for line in f_gz if not line.startswith("#")]

    with open(output_vcf, "w") as f_out:
        f_out.writelines(vcf_content)

    assert filecmp.cmp(output_vcf, expected_vcf), "VCF output differs from expected"


@pytest.mark.integration
def test_chr_prefixed(tmp_path):
    """Test chr prefix detection with 'chr' prefixed chromosomes."""
    # Get paths to required files
    project_root = get_project_root()
    src_data_dir = project_root / "src" / "data" / "numeric_chrs"
    # Using CLI commands instead of script paths

    # Input files
    truth_vcf = src_data_dir / "chrtruth.vcf"
    query_vcf = src_data_dir / "chrquery.vcf"
    fp_bed = src_data_dir / "chrfp.bed"
    reference = src_data_dir / "chrtest.fa"
    expected_vcf = src_data_dir / "chrexpected.vcf"
    expected_summary = src_data_dir / "expected.summary.csv"

    # Check that required files exist
    assert truth_vcf.exists(), f"Truth VCF {truth_vcf} not found"
    assert query_vcf.exists(), f"Query VCF {query_vcf} not found"
    assert fp_bed.exists(), f"FP BED {fp_bed} not found"
    assert reference.exists(), f"Reference {reference} not found"
    assert expected_vcf.exists(), f"Expected VCF {expected_vcf} not found"
    assert expected_summary.exists(), f"Expected summary {expected_summary} not found"

    # Output file paths
    output_prefix = tmp_path / "chr_test"
    output_vcf_gz = output_prefix.with_suffix(".vcf.gz")
    output_vcf = output_prefix.with_suffix(".vcf")
    output_summary = Path(str(output_prefix) + ".summary.csv")

    # Run hap.py on chr-prefixed files using CLI command
    cmd = [
        "hap",
        str(truth_vcf),
        str(query_vcf),
        "-f",
        str(fp_bed),
        "-o",
        str(output_prefix),
        "-X",
        "--reference",
        str(reference),
        "-V",
        "--force-interactive",
    ]

    result = subprocess.run(cmd, capture_output=True)
    assert (
        result.returncode == 0
    ), f"hap.py failed with chr-prefixed chromosomes: {result.stderr.decode()}"

    # Compare summary files
    compare_cmd = [
        sys.executable,
        str(compare_script),
        str(output_summary),
        str(expected_summary),
    ]

    compare_result = subprocess.run(compare_cmd, capture_output=True)
    assert compare_result.returncode == 0, "Summary comparison failed"

    # Compare VCF files
    with gzip.open(output_vcf_gz, "rt") as f_gz:
        vcf_content = [line for line in f_gz if not line.startswith("#")]

    with open(output_vcf, "w") as f_out:
        f_out.writelines(vcf_content)

    assert filecmp.cmp(output_vcf, expected_vcf), "VCF output differs from expected"


@pytest.mark.integration
def test_mixed_chr_prefix(tmp_path):
    """Test chr prefix detection with mixed chromosome naming (chr in truth, numeric in query)."""
    # Get paths to required files
    project_root = get_project_root()
    src_data_dir = project_root / "src" / "data" / "numeric_chrs"
    # Using CLI commands instead of script paths

    # Input files - mixing chr in truth with numeric in query
    truth_vcf = src_data_dir / "chrtruth.vcf"
    query_vcf = src_data_dir / "query.vcf"  # numeric chromosomes
    fp_bed = src_data_dir / "chrfp.bed"
    reference = src_data_dir / "chrtest.fa"
    expected_vcf = src_data_dir / "chrexpected.vcf"
    expected_summary = src_data_dir / "expected.summary.csv"

    # Check that required files exist
    assert truth_vcf.exists(), f"Truth VCF {truth_vcf} not found"
    assert query_vcf.exists(), f"Query VCF {query_vcf} not found"
    assert fp_bed.exists(), f"FP BED {fp_bed} not found"
    assert reference.exists(), f"Reference {reference} not found"
    assert expected_vcf.exists(), f"Expected VCF {expected_vcf} not found"
    assert expected_summary.exists(), f"Expected summary {expected_summary} not found"

    # Output file paths
    output_prefix = tmp_path / "mixed_test"
    output_vcf_gz = output_prefix.with_suffix(".vcf.gz")
    output_vcf = output_prefix.with_suffix(".vcf")
    output_summary = Path(str(output_prefix) + ".summary.csv")

    # Run hap.py with mixed chromosome naming using CLI command
    cmd = [
        "hap",
        str(truth_vcf),
        str(query_vcf),
        "-f",
        str(fp_bed),
        "-o",
        str(output_prefix),
        "-X",
        "--reference",
        str(reference),
        "-V",
        "--force-interactive",
    ]

    result = subprocess.run(cmd, capture_output=True)
    assert (
        result.returncode == 0
    ), f"hap.py failed with mixed chromosome naming: {result.stderr.decode()}"

    # Compare summary files
    compare_cmd = [
        sys.executable,
        str(compare_script),
        str(output_summary),
        str(expected_summary),
    ]

    compare_result = subprocess.run(compare_cmd, capture_output=True)
    assert compare_result.returncode == 0, "Summary comparison failed"

    # Compare VCF files
    with gzip.open(output_vcf_gz, "rt") as f_gz:
        vcf_content = [line for line in f_gz if not line.startswith("#")]

    with open(output_vcf, "w") as f_out:
        f_out.writelines(vcf_content)

    assert filecmp.cmp(output_vcf, expected_vcf), "VCF output differs from expected"
