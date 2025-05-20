"""
Integration tests for performance and consistency between simplecmp and hapcmp.
Migrated from src/sh/run_performance_test.sh
"""

import re

import pytest

from tests.utils import find_reference_file, get_bin_dir, get_example_dir, run_command


@pytest.mark.integration
@pytest.mark.cpp
def test_performance_vcf(tmp_path):
    """Test performance and consistency between simplecmp and hapcmp on VCF"""
    # Get paths to required files and tools
    bin_dir = get_bin_dir()
    example_dir = get_example_dir()
    reference = find_reference_file()

    # Define file paths for the test
    pg_vcf = example_dir / "PG_performance.vcf.gz"
    performance_vcf = example_dir / "performance.vcf.gz"

    # Define output paths
    output_vcf = tmp_path / "performance.vcf"
    output_bed = tmp_path / "performance.bed"

    # Check that xcmp binary exists
    xcmp_path = bin_dir / "xcmp"
    assert xcmp_path.exists(), f"Cannot find xcmp binary at {xcmp_path}"

    # Check that input files exist
    assert pg_vcf.exists(), f"PG_performance.vcf.gz not found at {pg_vcf}"
    assert (
        performance_vcf.exists()
    ), f"performance.vcf.gz not found at {performance_vcf}"
    assert (
        reference
    ), "Reference file not found. Set HGREF environment variable or use example reference."

    # Run xcmp with the same parameters as in the shell script
    cmd = [
        str(xcmp_path),
        str(pg_vcf),
        str(performance_vcf),
        "--always-hapcmp",
        "1",
        "-e",
        str(output_bed),
        "-o",
        str(output_vcf),
        "-r",
        reference,
        "-f",
        "0",
        "-n",
        "256",
    ]

    result = run_command(cmd)
    assert (
        result.returncode == 0
    ), f"xcmp failed with output: {result.stdout}\n{result.stderr}"

    # Check that output files were created
    assert output_bed.exists(), f"Output BED not generated: {output_bed}"

    # Parse the BED file to count different types of matches and mismatches
    with open(output_bed) as f:
        content = f.read()

    simple_mis_hc_match = len(re.findall(r"simple:match", content))
    simple_mis_hc_mis = len(re.findall(r"mismatch", content))
    simple_match_hc_mis = len(re.findall(r"suspicious", content))

    # Print results for debugging
    print(
        f"Mismatches in both: '{simple_mis_hc_mis}' (== cases not rescued through haplotype comparison)"
    )
    print(
        f"Matches in hapcmp but not simplecmp: '{simple_mis_hc_match}' (== cases rescued through haplotype comparison)"
    )
    print(
        f"Suspicious matches: '{simple_match_hc_mis}' (if this is >0, things are not good)"
    )

    # Check if there are no suspicious matches (which would be bad)
    assert (
        simple_match_hc_mis == 0
    ), "Consistency check failed: there are blocks for which simplecmp reports a match, but hapcmp reports a mismatch."


@pytest.mark.integration
@pytest.mark.cpp
def test_performance_gvcf(tmp_path):
    """Test performance and consistency between simplecmp and hapcmp on GVCF"""
    # Get paths to required files and tools
    bin_dir = get_bin_dir()
    example_dir = get_example_dir()
    reference = find_reference_file()

    # Define file paths for the test
    pg_vcf = example_dir / "PG_performance.vcf.gz"
    performance_vcf = example_dir / "performance.vcf.gz"

    # Define output paths
    output_vcf = tmp_path / "performance.vcf"
    output_bed = tmp_path / "performance.bed"

    # Check that xcmp binary exists
    xcmp_path = bin_dir / "xcmp"
    assert xcmp_path.exists(), f"Cannot find xcmp binary at {xcmp_path}"

    # Check that input files exist
    assert pg_vcf.exists(), f"PG_performance.vcf.gz not found at {pg_vcf}"
    assert (
        performance_vcf.exists()
    ), f"performance.vcf.gz not found at {performance_vcf}"
    assert (
        reference
    ), "Reference file not found. Set HGREF environment variable or use example reference."

    # Run xcmp with the same parameters as in the shell script (GVCF test)
    cmd = [
        str(xcmp_path),
        str(pg_vcf),
        str(performance_vcf),
        "-r",
        reference,
        "-f",
        "0",
        "-n",
        "256",
        "--always-hapcmp",
        "1",
        "-e",
        str(output_bed),
        "-o",
        str(output_vcf),
    ]

    result = run_command(cmd)
    assert (
        result.returncode == 0
    ), f"xcmp failed with output: {result.stdout}\n{result.stderr}"

    # Check that output files were created
    assert output_bed.exists(), f"Output BED not generated: {output_bed}"

    # Parse the BED file to count different types of matches and mismatches
    with open(output_bed) as f:
        content = f.read()

    simple_mis_hc_match = len(re.findall(r"simple:match", content))
    simple_mis_hc_mis = len(re.findall(r"mismatch", content))
    simple_match_hc_mis = len(re.findall(r"suspicious", content))

    # Print results for debugging
    print(
        f"Mismatches in both: '{simple_mis_hc_mis}' (== cases not rescued through haplotype comparison)"
    )
    print(
        f"Matches in hapcmp but not simplecmp: '{simple_mis_hc_match}' (== cases rescued through haplotype comparison)"
    )
    print(
        f"Suspicious matches: '{simple_match_hc_mis}' (if this is >0, things are not good)"
    )

    # Check if there are no suspicious matches (which would be bad)
    assert (
        simple_match_hc_mis == 0
    ), "Consistency check failed: there are blocks for which simplecmp reports a match, but hapcmp reports a mismatch."


@pytest.mark.integration
@pytest.mark.cpp
def test_performance_gvcf_comparison(tmp_path):
    """Test performance and consistency between simplecmp and hapcmp on GVCF"""
    # Get paths to required files and tools
    bin_dir = get_bin_dir()
    example_dir = get_example_dir()
    reference = find_reference_file()

    # Define file paths for the test
    pg_vcf = example_dir / "PG_performance.vcf.gz"
    performance_vcf = example_dir / "performance.vcf.gz"

    # Define output paths
    output_vcf = tmp_path / "performance_gvcf.vcf"
    output_bed = tmp_path / "performance_gvcf.bed"

    # Check that xcmp binary exists
    xcmp_path = bin_dir / "xcmp"
    assert xcmp_path.exists(), f"Cannot find xcmp binary at {xcmp_path}"

    # Check that input files exist
    assert pg_vcf.exists(), f"PG_performance.vcf.gz not found at {pg_vcf}"
    assert (
        performance_vcf.exists()
    ), f"performance.vcf.gz not found at {performance_vcf}"
    assert (
        reference
    ), "Reference file not found. Set HGREF environment variable or use example reference."

    # Run xcmp with the same parameters as in the shell script (GVCF test)
    cmd = [
        str(xcmp_path),
        str(pg_vcf),
        str(performance_vcf),
        "-r",
        reference,
        "-f",
        "0",
        "-n",
        "256",
        "--always-hapcmp",
        "1",
        "-e",
        str(output_bed),
        "-o",
        str(output_vcf),
    ]

    result = run_command(cmd)
    assert (
        result.returncode == 0
    ), f"xcmp failed with output: {result.stdout}\n{result.stderr}"

    # Check that output files were created
    assert output_bed.exists(), f"Output BED not generated: {output_bed}"

    # Parse the BED file to count different types of matches and mismatches
    with open(output_bed) as f:
        content = f.read()

    simple_mis_hc_match = len(re.findall(r"simple:match", content))
    simple_mis_hc_mis = len(re.findall(r"mismatch", content))
    simple_match_hc_mis = len(re.findall(r"suspicious", content))

    # Print results for debugging
    print(
        f"Mismatches in both: '{simple_mis_hc_mis}' (== cases not rescued through haplotype comparison)"
    )
    print(
        f"Matches in hapcmp but not simplecmp: '{simple_mis_hc_match}' (== cases rescued through haplotype comparison)"
    )
    print(
        f"Suspicious matches: '{simple_match_hc_mis}' (if this is >0, things are not good)"
    )

    # Check if there are no suspicious matches (which would be bad)
    assert (
        simple_match_hc_mis == 0
    ), "Consistency check failed: there are blocks for which simplecmp reports a match, but hapcmp reports a mismatch."
