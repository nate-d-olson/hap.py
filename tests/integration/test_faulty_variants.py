"""
Integration tests for faulty variant handling.
Migrated from src/sh/run_faulty_variant_test.sh
"""

import filecmp
import gzip
import subprocess
import sys
from pathlib import Path

import pytest


@pytest.mark.integration
def test_faulty_variant_handling(temp_dir):
    """Test handling of faulty variants."""
    # Get paths to reference files
    project_root = Path(__file__).parent.parent.parent
    src_data_dir = project_root / "src" / "data" / "open_indel"
    test_vcf = src_data_dir / "test.vcf"
    test_q_vcf = src_data_dir / "test_q.vcf"
    test_q_failure_vcf = src_data_dir / "test_q_failure.vcf"
    reference = src_data_dir / "test.fa"
    expected_vcf = src_data_dir / "expected.vcf"

    # Define path to hap.py script
    hap_py_script = project_root / "src" / "python" / "hap.py"

    # Output file paths
    output_prefix = Path(temp_dir) / "result"
    output_vcf_gz = output_prefix.with_suffix(".vcf.gz")
    output_vcf = output_prefix.with_suffix(".vcf")

    # Check that required files exist
    assert test_vcf.exists(), f"Test VCF file {test_vcf} not found"
    assert test_q_vcf.exists(), f"Test query VCF file {test_q_vcf} not found"
    assert (
        test_q_failure_vcf.exists()
    ), f"Test query failure VCF file {test_q_failure_vcf} not found"
    assert reference.exists(), f"Reference file {reference} not found"
    assert expected_vcf.exists(), f"Expected output file {expected_vcf} not found"

    # Test 1: hap.py with valid inputs
    cmd = [
        sys.executable,
        str(hap_py_script),
        str(test_vcf),
        str(test_q_vcf),
        "-o",
        str(output_prefix),
        "-X",
        "--reference",
        str(reference),
        "-l",
        "chrQ",
        "-V",
        "--force-interactive",
    ]

    result = subprocess.run(cmd, capture_output=True)
    assert (
        result.returncode == 0
    ), f"hap.py failed with valid inputs: {result.stderr.decode()}"

    # Decompress and compare output to expected
    with gzip.open(output_vcf_gz, "rt") as f_gz:
        vcf_content = [line for line in f_gz if not line.startswith("#")]

    with open(output_vcf, "w") as f_out:
        f_out.writelines(vcf_content)

    assert filecmp.cmp(output_vcf, expected_vcf), "Output differs from expected"

    # Test 2: hap.py with faulty inputs - should fail
    cmd = [
        sys.executable,
        str(hap_py_script),
        str(test_vcf),
        str(test_q_failure_vcf),
        "-o",
        str(output_prefix),
        "-X",
        "--reference",
        str(reference),
        "-l",
        "chrQ",
        "-V",
        "--force-interactive",
    ]

    result = subprocess.run(cmd, capture_output=True)
    assert result.returncode != 0, "hap.py did not fail with faulty inputs"


@pytest.mark.integration
def test_faulty_variant_pre_py(temp_dir):
    """Test handling of faulty variants in pre.py."""
    # Get paths to reference files
    project_root = Path(__file__).parent.parent.parent
    src_data_dir = project_root / "src" / "data"
    faulty_vcf = src_data_dir / "faulty.vcf"
    reference = src_data_dir / "chrQ.fa"

    # Define path to pre.py script
    pre_py_script = project_root / "src" / "python" / "pre.py"

    # Output file path
    output_file = Path(temp_dir) / "result.vcf"

    # Check that required files exist
    assert faulty_vcf.exists(), f"Faulty VCF file {faulty_vcf} not found"
    assert reference.exists(), f"Reference file {reference} not found"

    # Run pre.py with faulty input - should fail
    cmd = [
        sys.executable,
        str(pre_py_script),
        str(faulty_vcf),
        str(output_file),
        "--reference",
        str(reference),
        "-l",
        "chrQ",
        "--decompose",
        "--force-interactive",
    ]

    result = subprocess.run(cmd, capture_output=True)
    assert result.returncode != 0, "pre.py did not fail with faulty inputs"
