"""
Integration tests for left-shifting variants.
Migrated from src/sh/run_leftshift_test.sh
"""

import filecmp
import gzip
import os
import subprocess
import sys
from pathlib import Path

import pytest
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from utils import get_project_root


@pytest.mark.integration
def test_leftshift(tmp_path):
    """Test left-shifting functionality in hap.py."""
    # Get paths to required files
    project_root = get_project_root()
    src_data_dir = project_root / "src" / "data" / "leftshifting_example"
    hap_py_script = project_root / "src" / "python" / "hap.py"
    compare_script = project_root / "src" / "sh" / "compare_extended.py"

    # Input files
    truth_vcf = src_data_dir / "truth.vcf"
    query_vcf = src_data_dir / "query.vcf"
    reference = src_data_dir / "ref.fa"
    expected_vcf = src_data_dir / "expected.vcf"
    expected_extended = src_data_dir / "expected.extended.csv"

    # Check that required files exist
    assert truth_vcf.exists(), f"Truth VCF {truth_vcf} not found"
    assert query_vcf.exists(), f"Query VCF {query_vcf} not found"
    assert reference.exists(), f"Reference {reference} not found"
    assert expected_vcf.exists(), f"Expected VCF {expected_vcf} not found"
    assert expected_extended.exists(), (
        f"Expected extended CSV {expected_extended} not found"
    )

    # Output file paths
    output_prefix = tmp_path / "leftshift_test"
    output_vcf_gz = output_prefix.with_suffix(".vcf.gz")
    output_vcf = output_prefix.with_suffix(".vcf")
    output_extended = Path(str(output_prefix) + ".extended.csv")

    # Run hap.py with left-shifting
    cmd = [
        sys.executable,
        str(hap_py_script),
        str(truth_vcf),
        str(query_vcf),
        "-o",
        str(output_prefix),
        "-X",
        "--reference",
        str(reference),
        "-l",
        "chrT",
        "--preprocess-truth",
        "--leftshift",
        "--force-interactive",
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0, f"hap.py failed with left-shifting: {result.stderr}"

    # Compare extended files
    compare_cmd = [
        sys.executable,
        str(compare_script),
        str(output_extended),
        str(expected_extended),
    ]

    compare_result = subprocess.run(compare_cmd, capture_output=True, text=True)
    assert compare_result.returncode == 0, (
        f"Extended CSV comparison failed: {compare_result.stderr}"
    )

    # Compare VCF files - extract non-header lines from the gzipped VCF
    with gzip.open(output_vcf_gz, "rt") as f_gz:
        vcf_content = [line for line in f_gz if not line.startswith("#")]

    with open(output_vcf, "w") as f_out:
        f_out.writelines(vcf_content)

    assert filecmp.cmp(output_vcf, expected_vcf), "VCF output differs from expected"
