"""
Integration tests for variant decomposition functionality.
Migrated from src/sh/run_decomp_test.sh
"""

import os
import subprocess
from pathlib import Path

import pytest

from tests.utils import compare_summary_files, get_example_dir, run_command


@pytest.mark.integration
def test_decomp(tmp_path, rtg_executable, reference_file):
    """Test variant decomposition functionality"""
    # Skip test if reference file is not available
    if reference_file is None:
        pytest.skip("Reference file not available for testing")

    # Get paths to required files and tools
    example_dir = get_example_dir()

    # Define file paths for the test
    decomp_dir = example_dir / "decomp"
    truth_vcf = decomp_dir / "decomp_test.truth.vcf.gz"
    query_vcf = decomp_dir / "decomp_test.query.vcf.gz"
    conf_bed = decomp_dir / "decomp_test.conf.bed.gz"
    expected_summary = decomp_dir / "expected.summary.csv"
    expected_vcf = decomp_dir / "expected.vcf"

    # Define output path
    output_prefix = tmp_path / "decomp_test_out"

    assert truth_vcf.exists(), f"Truth VCF not found: {truth_vcf}"
    assert query_vcf.exists(), f"Query VCF not found: {query_vcf}"
    assert conf_bed.exists(), f"Confident regions BED not found: {conf_bed}"

    # Run hap.py with the same parameters as in the shell script
    cmd = [
        "hap.py",
        str(truth_vcf),
        str(query_vcf),
        "-f",
        str(conf_bed),
        "-r",
        reference_file,  # Add reference file
        "-o",
        str(output_prefix),
        "--preprocess-truth",
        "-X",
        "-V",
        "--force-interactive",
        "--engine-vcfeval-path",
        rtg_executable,  # Use fixture instead of hardcoded path
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
