"""
Integration tests for GiaB (Genome in a Bottle) functionality.
Migrated from src/sh/run_giab_test.sh
"""

import os
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
def test_small_giab_rtg(tmp_path):
    """Test small GiaB/RTG comparison"""
    # Get paths to required files and tools
    project_root = get_project_root()
    bin_dir = get_bin_dir()
    example_dir = get_example_dir()
    python_exe = get_python_executable()

    # Define file paths for the test
    giab_dir = example_dir / "GiaB"
    nist_vcf = giab_dir / "Complex_2ormoreindels_framerestoring_NIST2.19.ucsccoding.vcf"
    rtg_vcf = giab_dir / "Complex_2ormoreindels_framerestoring_RTG.ucsccoding.vcf"

    # Define output path
    output_prefix = tmp_path / "small_giab_test_out"

    # Determine the path to hap.py script
    hap_py_script = project_root / "src" / "python" / "hap.py"
    if not hap_py_script.exists():
        # Try alternative location
        hap_py_script = bin_dir / "hap.py"

    assert nist_vcf.exists(), f"NIST VCF not found: {nist_vcf}"
    assert rtg_vcf.exists(), f"RTG VCF not found: {rtg_vcf}"

    # Run hap.py with the same parameters as in the shell script
    cmd = [
        python_exe,
        str(hap_py_script),
        str(nist_vcf),
        str(rtg_vcf),
        "-o",
        str(output_prefix),
        "-X",
        "--force-interactive",
    ]

    result = run_command(cmd)
    assert (
        result.returncode == 0
    ), f"hap.py failed with output: {result.stdout}\n{result.stderr}"


@pytest.mark.integration
def test_large_giab_rtg_chr21(tmp_path):
    """Test large GiaB/RTG comparison on chromosome 21"""
    # Get paths to required files and tools
    project_root = get_project_root()
    bin_dir = get_bin_dir()
    example_dir = get_example_dir()
    python_exe = get_python_executable()

    # Define file paths for the test
    nist_indels_dir = example_dir / "NIST_indels"
    nist_vcf = nist_indels_dir / "Complex_1ormoreindels_NIST2.19.vcf.gz"
    rtg_vcf = nist_indels_dir / "Complex_1ormoreindels_RTG.vcf.gz"
    expected_summary = nist_indels_dir / "expected.summary.21.csv"

    # Define output path
    output_prefix = tmp_path / "large_giab_chr21_test_out"

    # Determine the path to hap.py script
    hap_py_script = project_root / "src" / "python" / "hap.py"
    if not hap_py_script.exists():
        # Try alternative location
        hap_py_script = bin_dir / "hap.py"

    assert nist_vcf.exists(), f"NIST VCF not found: {nist_vcf}"
    assert rtg_vcf.exists(), f"RTG VCF not found: {rtg_vcf}"
    assert expected_summary.exists(), f"Expected summary not found: {expected_summary}"

    # Run hap.py with the same parameters as in the shell script
    cmd = [
        python_exe,
        str(hap_py_script),
        str(nist_vcf),
        str(rtg_vcf),
        "-o",
        str(output_prefix),
        "-l",
        "chr21",
        "-X",
        "--verbose",
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
def test_large_giab_rtg_chr1(tmp_path):
    """Test large GiaB/RTG comparison on chromosome 1"""
    # Get paths to required files and tools
    project_root = get_project_root()
    bin_dir = get_bin_dir()
    example_dir = get_example_dir()
    python_exe = get_python_executable()

    # Define file paths for the test
    nist_indels_dir = example_dir / "NIST_indels"
    nist_vcf = nist_indels_dir / "Complex_1ormoreindels_NIST2.19.vcf.gz"
    rtg_vcf = nist_indels_dir / "Complex_1ormoreindels_RTG.vcf.gz"
    expected_summary = nist_indels_dir / "expected.summary.1.csv"

    # Define output path
    output_prefix = tmp_path / "large_giab_chr1_test_out"

    # Determine the path to hap.py script
    hap_py_script = project_root / "src" / "python" / "hap.py"
    if not hap_py_script.exists():
        # Try alternative location
        hap_py_script = bin_dir / "hap.py"

    assert nist_vcf.exists(), f"NIST VCF not found: {nist_vcf}"
    assert rtg_vcf.exists(), f"RTG VCF not found: {rtg_vcf}"
    assert expected_summary.exists(), f"Expected summary not found: {expected_summary}"

    # Run hap.py with the same parameters as in the shell script
    cmd = [
        python_exe,
        str(hap_py_script),
        str(nist_vcf),
        str(rtg_vcf),
        "-o",
        str(output_prefix),
        "-l",
        "chr1",
        "-X",
        "--verbose",
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
