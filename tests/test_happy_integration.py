"""
Integration tests for the hap.py CLI.

These tests check that the hap.py CLI can process VCF files correctly
and produce the expected output files.
"""

import os
import subprocess
import sys
from pathlib import Path


def test_happy_basic_cli(sample_vcf_files, sample_reference, tmp_path):
    """Test basic hap.py CLI functionality."""
    truth_vcf, query_vcf = sample_vcf_files
    ref_fasta = sample_reference

    # Create output directory
    output_prefix = str(tmp_path / "test_output")

    # Get the path to the hap.py script
    script_dir = Path(__file__).resolve().parent.parent
    happy_script = script_dir / "src" / "python" / "hap.py"

    # Run hap.py with mock environment to avoid requirement for actual C++ components
    env = os.environ.copy()
    env["HAPLO_USE_MOCK"] = "1"

    # Run the command
    cmd = [
        sys.executable,
        str(happy_script),
        "--force-interactive",  # Avoid SGE requirements
        truth_vcf,
        query_vcf,
        "-r",
        ref_fasta,
        "-o",
        output_prefix,
        "--engine",
        "vcfeval",
    ]

    result = subprocess.run(cmd, env=env, capture_output=True, text=True)

    # Check if command executed successfully
    assert result.returncode == 0, f"hap.py command failed: {result.stderr}"

    # Check if expected output files were created
    expected_files = [f"{output_prefix}.runinfo.json", f"{output_prefix}.summary.csv"]

    for expected_file in expected_files:
        assert os.path.exists(
            expected_file
        ), f"Expected output file {expected_file} not found"


def test_happy_with_bed_file(
    sample_vcf_files, sample_reference, sample_bed_file, tmp_path
):
    """Test hap.py CLI with BED file for stratification."""
    truth_vcf, query_vcf = sample_vcf_files
    ref_fasta = sample_reference
    bed_file = sample_bed_file

    # Create output directory
    output_prefix = str(tmp_path / "test_output_bed")

    # Get the path to the hap.py script
    script_dir = Path(__file__).resolve().parent.parent
    happy_script = script_dir / "src" / "python" / "hap.py"

    # Run hap.py with mock environment
    env = os.environ.copy()
    env["HAPLO_USE_MOCK"] = "1"

    # Run the command with BED file
    cmd = [
        sys.executable,
        str(happy_script),
        "--force-interactive",  # Avoid SGE requirements
        truth_vcf,
        query_vcf,
        "-r",
        ref_fasta,
        "-o",
        output_prefix,
        "-f",
        bed_file,  # Add confident regions BED file
        "--engine",
        "vcfeval",
    ]

    result = subprocess.run(cmd, env=env, capture_output=True, text=True)

    # Check if command executed successfully
    assert result.returncode == 0, f"hap.py command failed: {result.stderr}"

    # Check if expected output files were created
    expected_files = [f"{output_prefix}.runinfo.json", f"{output_prefix}.summary.csv"]

    for expected_file in expected_files:
        assert os.path.exists(
            expected_file
        ), f"Expected output file {expected_file} not found"


def test_error_handling(sample_vcf_files, sample_reference, tmp_path):
    """Test that hap.py handles errors gracefully."""
    truth_vcf, query_vcf = sample_vcf_files

    # Create a non-existent reference file
    nonexistent_ref = str(tmp_path / "nonexistent.fa")

    # Get the path to the hap.py script
    script_dir = Path(__file__).resolve().parent.parent
    happy_script = script_dir / "src" / "python" / "hap.py"

    # Run hap.py with mock environment
    env = os.environ.copy()
    env["HAPLO_USE_MOCK"] = "1"

    # Run the command with non-existent reference
    cmd = [
        sys.executable,
        str(happy_script),
        "--force-interactive",  # Avoid SGE requirements
        truth_vcf,
        query_vcf,
        "-r",
        nonexistent_ref,
        "-o",
        str(tmp_path / "error_output"),
    ]

    result = subprocess.run(cmd, env=env, capture_output=True, text=True)

    # Check that command failed as expected
    assert (
        result.returncode != 0
    ), "hap.py command should have failed with non-existent reference"

    # Check that the error message mentions the reference file
    assert (
        "reference" in result.stderr.lower()
    ), "Error message should mention reference file"
