"""
Integration tests for GVCF homref functionality.
Migrated from src/sh/run_gvcf_homref_test.sh
"""

import pytest

from tests.utils import (
    compress_and_index_vcf,
    get_example_dir,
    get_python_executable,
    run_shell_command,
)


@pytest.mark.integration
def test_gvcf_homref(tmp_path):
    """Test multimerge functionality with homref blocks."""
    # Get paths to required files and tools
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
        get_python_executable(),
        "-m",
        "hap_py.utils.multimerge",
        str(homref_vcf_gz),
        str(homref2_vcf_gz),
        "-o",
        str(output_vcf),
        "-r",
        str(reference_fa),
    ]

    cmd_str = " ".join(multimerge_cmd)
    returncode, _, stderr = run_shell_command(cmd_str)
    assert returncode == 0, f"multimerge with homref options failed: {stderr}"

    # Ensure the output VCF was created and contains data.
    assert output_vcf.exists(), "multimerge did not produce an output VCF"
    with open(output_vcf, encoding="utf-8") as f:
        lines = [line for line in f.readlines() if not line.startswith("#")]
    assert lines, "multimerge output VCF is empty"


@pytest.mark.integration
def test_gvcf_homref_with_variants(tmp_path):
    """Test multimerge functionality with homref blocks and variants."""
    # Get paths to required files and tools
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
        get_python_executable(),
        "-m",
        "hap_py.utils.multimerge",
        f"{call_merge_vcf_gz}:*",
        "-o",
        str(output_vcf),
        "-r",
        str(reference_fa),
    ]

    cmd_str = " ".join(multimerge_cmd)
    returncode, _, stderr = run_shell_command(cmd_str)
    assert returncode == 0, f"multimerge with variants failed: {stderr}"

    # Ensure output file has variant entries.
    assert output_vcf.exists(), "multimerge did not produce an output VCF"
    with open(output_vcf, encoding="utf-8") as f:
        lines = [line for line in f.readlines() if not line.startswith("#")]
    assert lines, "multimerge output VCF is empty"
