#!/usr/bin/env python3
"""
Unit tests for the Python hapcmp implementation.

This module tests the functionality of python_hapcmp.py to ensure
it correctly performs haplotype comparison.
"""

import sys  # Moved to top
from pathlib import Path  # Moved to top

# Add project root to sys.path to allow importing Haplo
project_root_path = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(project_root_path))

import os
import tempfile

import pytest
from Haplo.python_hapcmp import (
    HaploComparator,
    HaplotypeBlock,
)


@pytest.fixture
def reference_path():  # Reverted rename
    """Path to a small reference FASTA file for testing."""
    project_root = Path(__file__).resolve().parent.parent.parent
    return str(project_root / "example" / "chr21.fa")


@pytest.fixture
def example_vcf_paths():  # Reverted rename
    """Paths to example VCF files for testing."""
    project_root = Path(__file__).resolve().parent.parent.parent
    return (
        str(project_root / "example" / "hc.vcf.gz"),
        str(project_root / "example" / "PG_hc.vcf.gz"),
    )


@pytest.fixture
def example_bed_path():  # Reverted rename
    """Path to an example BED file defining regions."""
    project_root = Path(__file__).resolve().parent.parent.parent
    return str(project_root / "example" / "hc.bed")


@pytest.fixture
def temp_output_path():  # Reverted rename
    """Create a temporary file for output."""
    # pylint: disable=consider-using-with
    tmp_file = tempfile.NamedTemporaryFile(suffix=".bed", delete=False)
    tmp_path = tmp_file.name
    tmp_file.close()
    yield tmp_path
    # Clean up after the test
    if os.path.exists(tmp_path):
        os.unlink(tmp_path)


def test_hapcmp_initialization(reference_path):
    """Test initializing the HaploComparator class."""
    haplo_comparator = HaploComparator(reference_path)
    assert haplo_comparator.reference_path == reference_path
    assert haplo_comparator.max_haplotypes == 4096  # Default value
    assert not haplo_comparator.do_alignment  # Default value
    assert not haplo_comparator.output_sequences  # Default value
    assert not haplo_comparator.apply_filters  # Default value


def test_haplotype_block_initialization():
    """Test initializing a HaplotypeBlock."""
    block = HaplotypeBlock("chr1", 1000, 2000)
    assert block.chrom == "chr1"
    assert block.start == 1000
    assert block.end == 2000
    assert block.variants1 == []
    assert block.variants2 == []


@pytest.mark.parametrize("max_haplotypes", [1, 10, 100, 4096])
def test_set_max_haplotypes(reference_path, max_haplotypes):
    """Test setting maximum haplotypes via constructor."""
    haplo_comparator = HaploComparator(reference_path, max_haplotypes=max_haplotypes)
    assert haplo_comparator.max_haplotypes == max_haplotypes


def test_set_do_alignments(reference_path):
    """Test setting alignment option via constructor."""
    # Default is False
    haplo_comparator_default = HaploComparator(reference_path)
    assert not haplo_comparator_default.do_alignment

    # Set to True
    haplo_comparator_true = HaploComparator(reference_path, do_alignment=True)
    assert haplo_comparator_true.do_alignment


@pytest.mark.integration
def test_basic_comparison(reference_path, example_vcf_paths, temp_output_path):
    """Basic integration test for haplotype comparison."""
    vcf1, vcf2 = example_vcf_paths

    haplo_comparator = HaploComparator(reference_path, max_haplotypes=512)

    region_chrom = "chr21"
    region_start = 10000000
    region_end = 10001000

    temp_region_bed_path = ""
    # pylint: disable=consider-using-with
    tmp_bed_file = tempfile.NamedTemporaryFile(
        mode="w", suffix=".bed", delete=False, encoding="utf-8"
    )
    tmp_bed_file.write(f"{region_chrom}\\t{region_start}\\t{region_end}\\n")
    temp_region_bed_path = tmp_bed_file.name
    tmp_bed_file.close()

    try:
        stats = haplo_comparator.compare_files(
            file1=vcf1,
            sample1="",
            file2=vcf2,
            sample2="",
            regions_file=temp_region_bed_path,
            output_bed=temp_output_path,
        )

        assert stats["blocks_total"] == 1
        assert stats["blocks_processed"] == 1
        assert (
            stats["blocks_match"] + stats["blocks_mismatch"] + stats["blocks_error"]
        ) == 1

        with open(temp_output_path, encoding="utf-8") as f_out:
            lines = f_out.readlines()
            assert len(lines) == 2  # Header + 1 data line
            assert lines[0].strip() == "#CHROM\\tSTART\\tEND\\tSTATUS"
            parts = lines[1].strip().split("\\t")
            assert parts[0] == region_chrom
            assert int(parts[1]) == region_start
            assert int(parts[2]) == region_end
            assert parts[3] in ["match", "mismatch", "error"]

    finally:
        if os.path.exists(temp_region_bed_path):
            os.unlink(temp_region_bed_path)
