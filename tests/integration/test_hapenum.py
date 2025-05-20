"""
Integration tests for hapenum functionality.
Migrated from src/sh/run_hapenum_test.sh
"""

import filecmp
import subprocess

import pytest
from tests.utils import get_bin_dir, get_project_root, run_command


@pytest.mark.integration
@pytest.mark.cpp
def test_hapenum(tmp_path):
    """Test hapenum's ability to enumerate haplotypes and generate a graph."""
    # Get paths to required files and tools
    project_root = get_project_root()
    bin_dir = get_bin_dir()
    data_dir = project_root / "src" / "data"

    # Input files
    chrq_fa = data_dir / "chrQ.fa"
    refgraph1_vcf = data_dir / "refgraph1.vcf"
    expected_dot = data_dir / "expected_refgraph.dot"

    # Output files
    temp_vcf_gz = tmp_path / "refgraph1.vcf.gz"
    temp_dot = tmp_path / "temp.dot"
    temp_svg = tmp_path / "temp.dot.svg"

    # Make sure the input files exist
    assert chrq_fa.exists(), f"Reference file {chrq_fa} not found"
    assert refgraph1_vcf.exists(), f"VCF file {refgraph1_vcf} not found"
    assert expected_dot.exists(), f"Expected dot file {expected_dot} not found"

    # Compress and index the VCF file
    bgzip_cmd = f"cat {refgraph1_vcf} | bgzip > {temp_vcf_gz}"
    subprocess.run(bgzip_cmd, shell=True, check=True)

    tabix_cmd = f"tabix -p vcf -f {temp_vcf_gz}"
    subprocess.run(tabix_cmd, shell=True, check=True)

    # Run hapenum to generate the dot file
    hapenum_exe = bin_dir / "hapenum"
    hapenum_cmd = [
        str(hapenum_exe),
        "-r",
        str(chrq_fa),
        f"{temp_vcf_gz}:NA12877",
        "--output-dot",
        str(temp_dot),
        "-l",
        "chrQ",
    ]

    run_command(hapenum_cmd)

    # Generate the SVG for visual inspection if needed
    try:
        dot_cmd = ["dot", "-Tsvg", str(temp_dot), "-o", str(temp_svg)]
        run_command(dot_cmd)
    except subprocess.CalledProcessError:
        # Optional - if dot is not installed, still compare dot files
        pass

    # Compare the generated dot file with the expected one
    assert filecmp.cmp(temp_dot, expected_dot), (
        "Generated dot file differs from expected"
    )
