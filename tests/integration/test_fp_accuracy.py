"""
Integration tests for FP region accuracy in hap.py.
Migrated from src/sh/run_fp_accuracy_test.sh
"""

import gzip
from pathlib import Path

import pytest
from tests.utils import (
    compare_files,
    compare_summary_files,
    get_bin_dir,
    get_project_root,
    get_python_executable,
    run_shell_command,
)


def extract_and_filter_vcf(gz_file: Path, output_file: Path) -> None:
    """Extract a gzipped VCF file and filter out header lines.

    Args:
        gz_file: Path to gzipped VCF file
        output_file: Path to write filtered output
    """
    with gzip.open(gz_file, "rt") as f_in, open(output_file, "w") as f_out:
        for line in f_in:
            if not line.startswith("#"):
                f_out.write(line)


@pytest.mark.integration
def test_fp_region_accuracy(tmp_path):
    """Test if FP regions are processed accurately."""
    # Get paths to required files and tools
    bin_dir = get_bin_dir()
    project_root = get_project_root()
    python_exe = get_python_executable()

    # Set up paths to data files
    data_dir = project_root / "src" / "data" / "fp_region_accuracy"
    truth_vcf = data_dir / "truth.vcf"
    query_vcf = data_dir / "query.vcf"
    fp_bed = data_dir / "fp.bed"
    reference_fa = data_dir / "test.fa"
    expected_summary = data_dir / "expected.summary.csv"
    expected_vcf = data_dir / "expected.vcf"

    # Create output prefix in temporary directory
    output_prefix = tmp_path / "fp_accuracy_output"

    # Run hap.py with FP region
    happy_cmd = [
        python_exe,
        str(bin_dir / "hap.py"),
        str(truth_vcf),
        str(query_vcf),
        "-f",
        str(fp_bed),
        "-o",
        str(output_prefix),
        "-X",
        "--reference",
        str(reference_fa),
        "-l",
        "chrQ",
        "-V",
        "--force-interactive",
    ]

    cmd_str = " ".join(happy_cmd)
    returncode, _, stderr = run_shell_command(cmd_str)
    assert returncode == 0, f"hap.py failed with error: {stderr}"

    # Compare summary file with expected
    assert compare_summary_files(
        output_prefix.with_suffix(".summary.csv"), expected_summary
    ), "Summary files don't match"

    # Extract and filter VCF output
    output_vcf = tmp_path / "fp_accuracy_output.vcf"
    extract_and_filter_vcf(output_prefix.with_suffix(".vcf.gz"), output_vcf)

    # Compare filtered VCF with expected
    assert compare_files(output_vcf, expected_vcf), "Output VCF doesn't match expected"
