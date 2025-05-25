"""
Integration tests for hapcmp functionality.
Migrated from src/sh/run_hapcmp_test.sh
"""

import subprocess
from pathlib import Path
import pytest
from test_utils import get_project_root, get_example_data_dir, get_build_dir


@pytest.mark.integration
@pytest.mark.cpp
def test_hapcmp(example_data_dir, temp_dir):
    """Test hapcmp functionality."""
    # Get paths to reference files
    project_root = get_project_root()
    src_data_dir = project_root / "src" / "data"
    expected_file = src_data_dir / "expected_hapcmp.bed"
    result_file = Path(temp_dir) / "result_hapcmp.bed"

    # Find reference genome (this would typically be from environment variable)
    # For testing purposes, we'll use chr21.fa from example directory
    reference = get_example_data_dir() / "chr21.fa"

    # Define path to hapcmp binary
    hapcmp_bin = get_build_dir() / "bin" / "hapcmp"
    
    # Check that all required paths exist
    if not reference.exists():
        pytest.skip(f"Reference file not found: {reference}")
    if not hapcmp_bin.exists():
        pytest.skip(f"Hapcmp binary not found: {hapcmp_bin}")
    
    # Check for RTG tools dependency
    rtg_path = project_root / "libexec" / "rtg-tools-install"
    if not rtg_path.exists():
        pytest.skip("RTG tools not found. Skipping hapcmp test.")

    # Define input files
    hc_bed = Path(example_data_dir) / "hc.bed"
    hc_vcf = Path(example_data_dir) / "hc.vcf.gz"
    pg_hc_vcf = Path(example_data_dir) / "PG_hc.vcf.gz"

    # Skip test if reference doesn't exist
    if not reference.exists():
        pytest.skip(f"Reference file {reference} not found")

    # Check that the binary exists
    assert hapcmp_bin.exists(), f"hapcmp binary not found at {hapcmp_bin}"

    # Run the hapcmp command
    cmd = [
        str(hapcmp_bin),
        "-r",
        str(reference),
        str(hc_bed),
        str(hc_vcf),
        str(pg_hc_vcf),
        "--progress=0",
        "-n",
        "512",
        "--output-sequences=1",
        "--do-alignment=1",
        "-b",
        str(result_file),
    ]

    subprocess.run(cmd, check=True)

    # Check that the result file exists
    assert result_file.exists(), f"Result file {result_file} not created"

    # Compare result with expected output
    with open(expected_file) as f_expected:
        expected_content = f_expected.readlines()

    with open(result_file) as f_result:
        result_content = f_result.readlines()

    assert result_content == expected_content, "Output does not match expected result"
