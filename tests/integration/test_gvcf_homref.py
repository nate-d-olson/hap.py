"""
Integration tests for GVCF homref functionality.
Migrated from src/sh/run_gvcf_homref_test.sh
"""

from pathlib import Path

import pytest
from tests.utils import compare_files, get_bin_dir, get_example_dir, run_shell_command


def compress_and_index_vcf(vcf_path: Path) -> Path:
    """Compress and index a VCF file using bgzip and tabix.

    Args:
        vcf_path: Path to the VCF file to compress

    Returns:
        Path to the compressed VCF file
    """
    gz_path = vcf_path.with_suffix(".vcf.gz")

    # Compress the VCF
    cmd = f"cat {vcf_path} | bgzip > {gz_path}"
    returncode, _, stderr = run_shell_command(cmd)
    assert returncode == 0, f"Failed to compress VCF: {stderr}"

    # Index the compressed VCF
    cmd = f"tabix -p vcf {gz_path}"
    returncode, _, stderr = run_shell_command(cmd)
    assert returncode == 0, f"Failed to index VCF: {stderr}"

    return gz_path


@pytest.mark.integration
def test_gvcf_homref(tmp_path):
    """Test multimerge functionality with homref blocks."""
    # Get paths to required files and tools
    bin_dir = get_bin_dir()
    example_dir = get_example_dir()

    # Set up paths
    homref_dir = example_dir / "homref"
    homref_vcf = homref_dir / "homref.vcf"
    homref2_vcf = homref_dir / "homref2.vcf"
    reference_fa = example_dir / "chr21.fa"

    # Compress and index VCF files
    homref_vcf_gz = compress_and_index_vcf(homref_vcf)
    homref2_vcf_gz = compress_and_index_vcf(homref2_vcf)

    # Create temporary output file
    output_vcf = tmp_path / "homref_output.vcf"

    # Run multimerge with homref options
    multimerge_cmd = [
        str(bin_dir / "multimerge"),
        str(homref_vcf_gz),
        str(homref2_vcf_gz),
        "-o",
        str(output_vcf),
        "-r",
        str(reference_fa),
        "--trimalleles=1",
        "--merge-by-location=1",
        "--homref-split=1",
        "--unique-alleles=1",
        "--calls-only=0",
    ]

    cmd_str = " ".join(multimerge_cmd)
    returncode, _, stderr = run_shell_command(cmd_str)
    assert returncode == 0, f"multimerge with homref options failed: {stderr}"

    # Note: The original test commented out the comparison with expected_merge.vcf
    # because of a known issue. We'll skip the comparison as well, but include
    # a comment about it.

    # TODO: This part of the test is disabled in the original shell script
    # because of an issue in VariantReader logic. Once fixed, uncomment this.
    #
    # expected_vcf = homref_dir / "expected_merge.vcf"
    # assert compare_files(
    #     output_vcf, expected_vcf, ignore_comments=True
    # ), "Homref merge output doesn't match expected"


@pytest.mark.integration
def test_gvcf_homref_with_variants(tmp_path):
    """Test multimerge functionality with homref blocks and variants."""
    # Get paths to required files and tools
    bin_dir = get_bin_dir()
    example_dir = get_example_dir()

    # Set up paths
    callsonly_dir = example_dir / "callsonly"
    call_merge_vcf = callsonly_dir / "call_merge.vcf"
    reference_fa = example_dir / "chr21.fa"

    # Compress and index VCF files
    call_merge_vcf_gz = compress_and_index_vcf(call_merge_vcf)

    # Create temporary output file
    output_vcf = tmp_path / "homref_variants_output.vcf"

    # Run multimerge with homref and variants options
    multimerge_cmd = [
        str(bin_dir / "multimerge"),
        f"{call_merge_vcf_gz}:*",
        "-o",
        str(output_vcf),
        "-r",
        str(reference_fa),
        "--process-full=1",
        "--process-formats=1",
    ]

    cmd_str = " ".join(multimerge_cmd)
    returncode, _, stderr = run_shell_command(cmd_str)
    assert returncode == 0, f"multimerge with variants failed: {stderr}"

    # Compare with expected output
    expected_vcf = callsonly_dir / "expected_callsonly.vcf"
    assert compare_files(
        output_vcf, expected_vcf
    ), "Homref+variants output doesn't match expected"
