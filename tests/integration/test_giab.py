"""
Integration tests for GiaB (Genome in a Bottle) functionality.
Migrated from src/sh/run_giab_test.sh
"""

import os
from pathlib import Path

import pytest

from tests.utils import (
    compare_summary_files,
    get_example_dir,
    run_command,
)


@pytest.mark.integration
def test_small_giab_rtg(tmp_path, rtg_executable, reference_file):
    """Test small GiaB/RTG comparison"""
    # Skip test if reference file is not available
    if reference_file is None:
        pytest.skip("Reference file not available for testing")

    # Get paths to required files and tools
    example_dir = get_example_dir()

    # Define file paths for the test
    giab_dir = example_dir / "GiaB"
    nist_vcf = giab_dir / (
        "Complex_2ormoreindels_framerestoring_NIST2.19.ucsccoding.vcf"
    )
    rtg_vcf = giab_dir / ("Complex_2ormoreindels_framerestoring_RTG.ucsccoding.vcf")

    # Define output path
    output_prefix = tmp_path / "small_giab_test_out"

    assert nist_vcf.exists(), f"NIST VCF not found: {nist_vcf}"
    assert rtg_vcf.exists(), f"RTG VCF not found: {rtg_vcf}"

    # Run hap.py with the same parameters as in the shell script
    cmd = [
        "hap.py",
        str(nist_vcf),
        str(rtg_vcf),
        "-r",
        reference_file,  # Add reference file
        "-o",
        str(output_prefix),
        "-X",
        "--force-interactive",
        "--engine-vcfeval-path",
        rtg_executable,  # Add RTG path
    ]

    result = run_command(cmd)
    assert (
        result.returncode == 0
    ), f"hap.py failed with output: {result.stdout}\n{result.stderr}"


@pytest.mark.integration
def test_large_giab_rtg_chr21(tmp_path, rtg_executable, reference_file):
    """Test large GiaB/RTG comparison on chromosome 21"""
    # Skip test if reference file is not available
    if reference_file is None:
        pytest.skip("Reference file not available for testing")

    # Get paths to required files and tools
    example_dir = get_example_dir()

    # Define file paths for the test
    nist_indels_dir = example_dir / "NIST_indels"
    nist_vcf = nist_indels_dir / "Complex_1ormoreindels_NIST2.19.vcf.gz"
    rtg_vcf = nist_indels_dir / "Complex_1ormoreindels_RTG.vcf.gz"
    expected_summary = nist_indels_dir / "expected.summary.21.csv"

    # Define output path
    output_prefix = tmp_path / "large_giab_chr21_test_out"

    assert nist_vcf.exists(), f"NIST VCF not found: {nist_vcf}"
    assert rtg_vcf.exists(), f"RTG VCF not found: {rtg_vcf}"
    assert expected_summary.exists(), f"Expected summary not found: {expected_summary}"

    # Run hap.py with the same parameters as in the shell script
    cmd = [
        "hap.py",
        str(nist_vcf),
        str(rtg_vcf),
        "-r",
        reference_file,  # Add reference file
        "-o",
        str(output_prefix),
        "-l",
        "chr21",
        "-X",
        "--verbose",
        "--force-interactive",
        "--engine-vcfeval-path",
        rtg_executable,  # Add RTG path
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
def test_large_giab_rtg_chr1(tmp_path, rtg_executable, reference_file):
    """Test large GiaB/RTG comparison on chromosome 1"""
    # Skip test if reference file is not available
    if reference_file is None:
        pytest.skip("Reference file not available for testing")

    # Get paths to required files and tools
    example_dir = get_example_dir()

    # Define file paths for the test
    nist_indels_dir = example_dir / "NIST_indels"
    nist_vcf = nist_indels_dir / "Complex_1ormoreindels_NIST2.19.vcf.gz"
    rtg_vcf = nist_indels_dir / "Complex_1ormoreindels_RTG.vcf.gz"
    expected_summary = nist_indels_dir / "expected.summary.1.csv"

    # Define output path
    output_prefix = tmp_path / "large_giab_chr1_test_out"

    assert nist_vcf.exists(), f"NIST VCF not found: {nist_vcf}"
    assert rtg_vcf.exists(), f"RTG VCF not found: {rtg_vcf}"
    assert expected_summary.exists(), f"Expected summary not found: {expected_summary}"

    # Run hap.py with the same parameters as in the shell script
    cmd = [
        "hap.py",
        str(nist_vcf),
        str(rtg_vcf),
        "-r",
        reference_file,  # Add reference file
        "-o",
        str(output_prefix),
        "-l",
        "chr1",
        "-X",
        "--verbose",
        "--force-interactive",
        "--engine-vcfeval-path",
        rtg_executable,  # Add RTG path
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
