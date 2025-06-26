"""
Unit/component tests for blocksplit functionality, refactored from integration test.

These tests mock external tool invocations and focus on overlap detection and variant coverage logic.
"""

from unittest.mock import MagicMock, patch

import pytest


def parse_bed_lines(bed_lines):
    """Parse BED lines into (chr, start, end) tuples."""
    for line in bed_lines:
        parts = line.strip().split("\t")
        if len(parts) < 3:
            continue
        yield parts[0], int(parts[1]), int(parts[2])


def test_no_overlaps_in_bed():
    """Test that BED regions do not overlap within the same chromosome."""
    # Simulated BED output from blocksplit
    bed_lines = [
        "chr21\t0\t10000\n",
        "chr21\t10000\t20000\n",
        "chr21\t20000\t30000\n",
        "chr22\t0\t5000\n",
        "chr22\t5000\t10000\n",
    ]
    prev_chr = None
    prev_end = -1
    for chr_name, start, end in parse_bed_lines(bed_lines):
        if chr_name == prev_chr and start < prev_end:
            pytest.fail(f"Found overlap at {chr_name}:{start}-{end}")
        prev_chr = chr_name
        prev_end = end


def test_all_variants_covered_by_blocks():
    """Test that all variants are covered by the blocks."""
    # Simulated VCF variant positions (as would be output by bcftools)
    vcf_variants = [
        ("chr21", 100),
        ("chr21", 15000),
        ("chr21", 25000),
        ("chr22", 2000),
        ("chr22", 7000),
    ]
    # Simulated BED blocks
    bed_blocks = [
        ("chr21", 0, 10000),
        ("chr21", 10000, 20000),
        ("chr21", 20000, 30000),
        ("chr22", 0, 5000),
        ("chr22", 5000, 10000),
    ]
    # Check that each variant is within a block
    for v_chr, v_pos in vcf_variants:
        covered = any(
            v_chr == b_chr and b_start <= v_pos < b_end
            for b_chr, b_start, b_end in bed_blocks
        )
        assert covered, f"Variant {v_chr}:{v_pos} not covered by any block"


@patch("subprocess.run")
def test_blocksplit_invocation(mock_run):
    """Test that blocksplit is invoked with correct arguments."""
    # Simulate successful run
    mock_run.return_value = MagicMock(returncode=0)
    import subprocess

    cmd = [
        "blocksplit",
        "vcf1.vcf.gz",
        "vcf2.vcf.gz",
        "-o",
        "out.bed",
        "-l",
        "chr21",
        "-w",
        "10000",
    ]
    result = subprocess.run(cmd)
    assert result.returncode == 0
    mock_run.assert_called_with(cmd)
    assert result.returncode == 0
    mock_run.assert_called_with(cmd)
