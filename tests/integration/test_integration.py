"""
Integration tests for hap.py with different VCF inputs.
Migrated from src/sh/run_integration_test.sh
"""

import os
import shutil
import gzip
import pytest
from pathlib import Path

from tests.utils import (
    check_vcfeval_availability,
    compare_files,
    compare_summary_files,
    find_reference_file,
    get_bin_dir,
    get_example_dir,
    get_python_executable,
    run_shell_command,
)


def gunzip_and_filter_headers(input_gz: Path, output_file: Path) -> None:
    """Extract a gzipped file and filter out header lines.

    Args:
        input_gz: Path to gzipped input file
        output_file: Path to write filtered output
    """
    with gzip.open(input_gz, "rt") as f_in, open(output_file, "w") as f_out:
        for line in f_in:
            if not line.startswith("#"):
                f_out.write(line)


def compare_vcf_files(file1: Path, file2: Path) -> bool:
    """Compare two VCF files ignoring headers.

    Args:
        file1: First VCF file
        file2: Second VCF file

    Returns:
        True if the files match, False otherwise
    """
    # Use grep to filter out header lines (lines starting with #)
    cmd = f"diff -I fileDate -I source_version {file1} {file2}"
    returncode, stdout, stderr = run_shell_command(cmd)
    return returncode == 0


@pytest.mark.integration
@pytest.mark.slow
def test_integration(tmp_path):
    """Test hap.py integration with different input configurations."""
    # Get paths to required files and tools
    bin_dir = get_bin_dir()
    example_dir = get_example_dir()
    python_exe = get_python_executable()

    # Set up paths to example files
    integration_dir = example_dir / "integration"
    lhs_vcf = integration_dir / "integrationtest_lhs.vcf"
    rhs_vcf = integration_dir / "integrationtest_rhs.vcf"
    empty_vcf = integration_dir / "integrationtest_empty.vcf"

    # Locate reference genome
    reference_fa = find_reference_file()
    assert reference_fa is not None, "Reference genome not found"

    # Create temporary output prefix
    output_prefix = tmp_path / "integration_output"

    # Prepare bgzipped and indexed VCF files
    lhs_vcf_gz = tmp_path / "integrationtest_lhs.vcf.gz"
    rhs_vcf_gz = tmp_path / "integrationtest_rhs.vcf.gz"
    empty_vcf_gz = integration_dir / "integrationtest_empty.vcf.gz"

    # Compress and index VCF files if not already compressed
    if not empty_vcf_gz.exists():
        cmd = f"cat {empty_vcf} | bgzip > {empty_vcf_gz}"
        run_shell_command(cmd)
        cmd = f"tabix -f -p vcf {empty_vcf_gz}"
        run_shell_command(cmd)

    # Compress and index the input VCFs
    cmd = f"cat {lhs_vcf} | bgzip > {lhs_vcf_gz}"
    run_shell_command(cmd)
    cmd = f"tabix -f -p vcf {lhs_vcf_gz}"
    run_shell_command(cmd)

    cmd = f"cat {rhs_vcf} | bgzip > {rhs_vcf_gz}"
    run_shell_command(cmd)
    cmd = f"tabix -f -p vcf {rhs_vcf_gz}"
    run_shell_command(cmd)

    # Test multimerge functionality
    multimerge_output = tmp_path / "multimerge_output.vcf"
    multimerge_cmd = [
        str(bin_dir / "multimerge"),
        str(lhs_vcf_gz),
        str(rhs_vcf_gz),
        "-o",
        str(multimerge_output),
        "-r",
        reference_fa,
        "--process-full",
        "1",
    ]

    cmd_str = " ".join(multimerge_cmd)
    returncode, _, stderr = run_shell_command(cmd_str)
    assert returncode == 0, f"multimerge failed with error: {stderr}"

    # Compare to expected merged file
    expected_merged_vcf = integration_dir / "integrationtest_merged.vcf"
    assert compare_vcf_files(
        multimerge_output, expected_merged_vcf
    ), "Merged VCF does not match expected output"

    # Test hap.py with empty truth file
    empty_truth_output = output_prefix.with_suffix(".e0")
    empty_truth_cmd = [
        python_exe,
        str(bin_dir / "hap.py"),
        "-l",
        "chr21",
        str(empty_vcf),
        str(rhs_vcf_gz),
        "-o",
        str(empty_truth_output),
        "-V",
        "-X",
        "--output-vtc",
        "--force-interactive",
    ]

    cmd_str = " ".join(empty_truth_cmd)
    returncode, _, stderr = run_shell_command(cmd_str)
    assert returncode == 0, f"hap.py with empty truth file failed: {stderr}"

    # Test hap.py with empty query file
    empty_query_output = output_prefix.with_suffix(".e1")
    empty_query_cmd = [
        python_exe,
        str(bin_dir / "hap.py"),
        "-l",
        "chr21",
        str(lhs_vcf_gz),
        str(empty_vcf),
        "-o",
        str(empty_query_output),
        "-V",
        "-X",
        "--output-vtc",
        "--force-interactive",
    ]

    cmd_str = " ".join(empty_query_cmd)
    returncode, _, stderr = run_shell_command(cmd_str)
    assert returncode == 0, f"hap.py with empty query file failed: {stderr}"

    # Test standard hap.py
    standard_cmd = [
        python_exe,
        str(bin_dir / "hap.py"),
        "-l",
        "chr21",
        str(lhs_vcf_gz),
        str(rhs_vcf_gz),
        "-o",
        str(output_prefix),
        "-V",
        "-X",
        "--output-vtc",
        "--force-interactive",
    ]

    cmd_str = " ".join(standard_cmd)
    returncode, _, stderr = run_shell_command(cmd_str)
    assert returncode == 0, f"standard hap.py failed: {stderr}"

    # Extract and filter standard output
    standard_output_filtered = tmp_path / "integration_output.vcf"
    gunzip_and_filter_headers(
        output_prefix.with_suffix(".vcf.gz"), standard_output_filtered
    )

    # Test unhappy mode
    unhappy_output = output_prefix.with_suffix(".unhappy")
    unhappy_cmd = [
        python_exe,
        str(bin_dir / "hap.py"),
        "-l",
        "chr21",
        str(lhs_vcf_gz),
        str(rhs_vcf_gz),
        "-o",
        str(unhappy_output),
        "-V",
        "-X",
        "--output-vtc",
        "--force-interactive",
        "--unhappy",
    ]

    cmd_str = " ".join(unhappy_cmd)
    returncode, _, stderr = run_shell_command(cmd_str)
    assert returncode == 0, f"hap.py unhappy mode failed: {stderr}"

    # Extract and filter unhappy output
    unhappy_output_filtered = tmp_path / "integration_output.unhappy.vcf"
    gunzip_and_filter_headers(
        unhappy_output.with_suffix(".vcf.gz"), unhappy_output_filtered
    )

    # Test pass-only mode
    pass_output = output_prefix.with_suffix(".pass")
    pass_cmd = [
        python_exe,
        str(bin_dir / "hap.py"),
        "-l",
        "chr21",
        str(lhs_vcf_gz),
        str(rhs_vcf_gz),
        "-o",
        str(pass_output),
        "-V",
        "-X",
        "--pass-only",
        "--force-interactive",
    ]

    cmd_str = " ".join(pass_cmd)
    returncode, _, stderr = run_shell_command(cmd_str)
    assert returncode == 0, f"hap.py pass-only mode failed: {stderr}"

    # Extract and filter pass-only output
    pass_output_filtered = tmp_path / "integration_output.pass.vcf"
    gunzip_and_filter_headers(pass_output.with_suffix(".vcf.gz"), pass_output_filtered)

    # Compare standard output with expected
    expected_vcf = integration_dir / "integrationtest.vcf"
    assert compare_vcf_files(
        standard_output_filtered, expected_vcf
    ), "Standard output doesn't match expected"

    # Compare unhappy output with expected
    expected_unhappy_vcf = integration_dir / "integrationtest.unhappy.vcf"
    assert compare_vcf_files(
        unhappy_output_filtered, expected_unhappy_vcf
    ), "Unhappy output doesn't match expected"

    # Compare pass-only output with expected
    expected_pass_vcf = integration_dir / "integrationtest.pass.vcf"
    assert compare_vcf_files(
        pass_output_filtered, expected_pass_vcf
    ), "Pass-only output doesn't match expected"

    # Compare summary files
    expected_summary = integration_dir / "integrationtest.summary.csv"
    assert compare_summary_files(
        output_prefix.with_suffix(".summary.csv"), expected_summary
    ), "Summary file doesn't match expected"
