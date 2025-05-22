"""
Integration tests for other VCF functionality.
Migrated from src/sh/run_other_vcf_tests.sh

This includes two tests:
1. Variant filtering test
2. Haploid variants test
"""

import os
import subprocess
from pathlib import Path

import pytest

from tests.utils import (
    compare_summary_files,
    get_bin_dir,
    get_example_dir,
    get_project_root,
    get_python_executable,
    run_command,
)


@pytest.mark.integration
def test_variant_filtering(tmp_path):
    """Test reading and detecting problematic records"""
    # Get paths to required files and tools
    project_root = get_project_root()
    bin_dir = get_bin_dir()
    python_exe = get_python_executable()

    # Define file paths for the test
    data_dir = project_root / "src" / "data"
    per_sample_ft_lhs_vcf = data_dir / "per_sample_ft_lhs.vcf"
    per_sample_ft_rhs_vcf = data_dir / "per_sample_ft_rhs.vcf"
    reference = data_dir / "chrQ.fa"
    expected_summary = data_dir / "per_sample_ft_summary.csv"

    # Define output path
    output_prefix = tmp_path / "filtering_test_out"

    # Determine the path to hap.py script
    hap_py_script = project_root / "src" / "python" / "hap.py"
    if not hap_py_script.exists():
        # Try alternative location
        hap_py_script = bin_dir / "hap.py"

    assert per_sample_ft_lhs_vcf.exists(), f"LHS VCF not found: {per_sample_ft_lhs_vcf}"
    assert per_sample_ft_rhs_vcf.exists(), f"RHS VCF not found: {per_sample_ft_rhs_vcf}"
    assert reference.exists(), f"Reference file not found: {reference}"
    assert expected_summary.exists(), f"Expected summary not found: {expected_summary}"

    # Run hap.py with the same parameters as in the shell script
    cmd = [
        python_exe,
        str(hap_py_script),
        str(per_sample_ft_lhs_vcf),
        str(per_sample_ft_rhs_vcf),
        "-o",
        str(output_prefix),
        "--reference",
        str(reference),
        "--force-interactive",
    ]

    result = run_command(cmd)
    assert (
        result.returncode == 0
    ), f"hap.py failed with output: {result.stdout}\n{result.stderr}"

    # Check summary file
    output_summary = str(output_prefix) + ".summary.csv"
    assert os.path.exists(
        output_summary
    ), f"Output summary not generated: {output_summary}"
    assert compare_summary_files(
        Path(output_summary), expected_summary
    ), f"Summary files differ: {output_summary} vs {expected_summary}"


@pytest.mark.integration
def test_haploid_variants(tmp_path):
    """Test handling of haploid records"""
    # Get paths to required files and tools
    project_root = get_project_root()
    bin_dir = get_bin_dir()
    example_dir = get_example_dir()
    python_exe = get_python_executable()

    # Define file paths for the test
    haploid_dir = example_dir / "haploid"
    truth_vcf = haploid_dir / "truth.vcf"
    query_vcf = haploid_dir / "query.vcf"
    reference = haploid_dir / "test.fa"
    expected_vcf = haploid_dir / "expected.vcf"

    # Define output path
    output_prefix = tmp_path / "haploid_test_out"

    # Determine the path to hap.py script
    hap_py_script = project_root / "src" / "python" / "hap.py"
    if not hap_py_script.exists():
        # Try alternative location
        hap_py_script = bin_dir / "hap.py"

    assert truth_vcf.exists(), f"Truth VCF not found: {truth_vcf}"
    assert query_vcf.exists(), f"Query VCF not found: {query_vcf}"
    assert reference.exists(), f"Reference file not found: {reference}"
    assert expected_vcf.exists(), f"Expected VCF not found: {expected_vcf}"

    # Run hap.py with the same parameters as in the shell script
    cmd = [
        python_exe,
        str(hap_py_script),
        str(truth_vcf),
        str(query_vcf),
        "-o",
        str(output_prefix),
        "-r",
        str(reference),
        "--force-interactive",
    ]

    result = run_command(cmd)
    assert (
        result.returncode == 0
    ), f"hap.py failed with output: {result.stdout}\n{result.stderr}"

    # Check VCF output - need to gunzip first
    output_vcf_gz = str(output_prefix) + ".vcf.gz"
    assert os.path.exists(output_vcf_gz), f"Output VCF not generated: {output_vcf_gz}"

    # Gunzip and compare with expected
    output_vcf = str(output_prefix) + ".vcf"
    gunzip_cmd = ["gunzip", "-c", output_vcf_gz]
    with open(output_vcf, "w") as f:
        result = subprocess.run(gunzip_cmd, stdout=f, check=True)

    # Compare VCF files ignoring header lines (lines starting with #)
    def compare_vcf_files(file1, file2):
        with open(file1) as f1, open(file2) as f2:
            lines1 = [line for line in f1 if not line.startswith("#")]
            lines2 = [line for line in f2 if not line.startswith("#")]
            return lines1 == lines2

    assert compare_vcf_files(
        output_vcf, expected_vcf
    ), f"VCF files differ: {output_vcf} vs {expected_vcf}"
