"""
Integration tests for path traversal functionality.
Migrated from src/sh/run_pathtraversal_test.sh
"""

import os
import subprocess
from pathlib import Path

import pytest
from tests.utils import (
    compare_summary_files,
    get_project_root,
    get_python_executable,
    run_command,
)


@pytest.mark.integration
def test_pathtraversal(tmp_path):
    """Test hap.py's handling of path traversals."""
    # Get paths to required files and tools
    project_root = get_project_root()
    src_dir = project_root / "src"
    data_dir = src_dir / "data" / "pathtraversal"
    python_exe = get_python_executable()
    happy_script = src_dir / "python" / "hap.py"
    compare_script = src_dir / "sh" / "compare_summaries.py"

    # Input files
    test_vcf = data_dir / "test.vcf"
    test2_vcf = data_dir / "test2.vcf"
    reference_fa = data_dir / "test.fa"
    expected_summary = data_dir / "expected.summary.csv"

    # Output files
    output_prefix = tmp_path / "output"
    output_summary = f"{output_prefix}.summary.csv"

    # Verify input files exist
    assert test_vcf.exists(), f"Test VCF file {test_vcf} not found"
    assert test2_vcf.exists(), f"Test2 VCF file {test2_vcf} not found"
    assert reference_fa.exists(), f"Reference file {reference_fa} not found"
    assert (
        expected_summary.exists()
    ), f"Expected summary file {expected_summary} not found"

    # Run hap.py
    happy_cmd = [
        python_exe,
        str(happy_script),
        str(test_vcf),
        str(test2_vcf),
        "-o",
        str(output_prefix),
        "-X",
        "--reference",
        str(reference_fa),
        "-l",
        "chrQ",
        "--force-interactive",
    ]

    run_command(happy_cmd)

    # Compare summaries
    compare_cmd = [
        python_exe,
        str(compare_script),
        output_summary,
        str(expected_summary),
    ]

    try:
        subprocess.run(compare_cmd, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        pytest.fail(f"Summary comparison failed: {e.stderr}")

    # Alternative direct comparison approach
    assert os.path.exists(
        output_summary
    ), f"Output summary file {output_summary} not found"
    assert compare_summary_files(
        output_summary, expected_summary
    ), "Summary files differ"
