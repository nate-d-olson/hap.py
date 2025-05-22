#!/usr/bin/env python3
"""
Unit tests for the Python preprocess implementation.
"""

import os
import sys
import tempfile
from pathlib import Path

import pysam
import pytest

# Add src to sys.path to allow importing hap_py
project_root_path = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(project_root_path / "src"))

from hap_py.haplo.python_preprocess import DecomposeLevel, PreprocessEngine


@pytest.fixture
def reference_path():
    """Path to a small reference FASTA file for testing."""
    repo_root = Path(__file__).resolve().parent.parent.parent
    path = repo_root / "example" / "example.fa"
    assert path.exists(), f"Reference FASTA not found at {path}"
    # Ensure the .fai index exists or can be created by pysam
    if not (path.with_suffix(path.suffix + ".fai")).exists():
        try:
            pysam.faidx(str(path))
        except pysam.SamtoolsError as e:
            # If indexing fails, it might indicate a problem with the FASTA file itself
            # or permissions. For now, we'll raise an error to make it visible.
            raise FileNotFoundError(
                f"Failed to find or create FASTA index for {path}: {e}"
            )
    return str(path)


@pytest.fixture
def example_vcf_path():
    """Path to an example VCF file for testing."""
    repo_root = Path(__file__).resolve().parent.parent.parent
    path = repo_root / "example" / "example.vcf.gz"
    assert path.exists(), f"Example VCF not found at {path}"
    # Ensure the .tbi index exists or can be created by pysam
    if not (path.with_suffix(path.suffix + ".tbi")).exists():
        try:
            pysam.tabix_index(str(path), preset="vcf")
        except pysam.SamtoolsError as e:
            raise FileNotFoundError(
                f"Failed to find or create VCF index for {path}: {e}"
            )
    return str(path)


@pytest.fixture
def temp_output_path():
    """Create a temporary file for output."""
    fd, path = tempfile.mkstemp(suffix=".vcf.gz")
    os.close(fd)
    yield path
    # Clean up
    if os.path.exists(path):
        os.remove(path)
    if os.path.exists(path + ".tbi"):
        os.remove(path + ".tbi")


def test_preprocess_initialization(reference_path, example_vcf_path):
    """Test initializing the PreprocessEngine."""
    # Ensure reference_path and example_vcf_path are valid before initializing
    assert os.path.exists(reference_path), f"Reference FASTA missing: {reference_path}"
    assert os.path.exists(example_vcf_path), f"Example VCF missing: {example_vcf_path}"

    engine = PreprocessEngine(
        input_vcf=example_vcf_path, reference_fasta=reference_path
    )

    assert engine.input_vcf == example_vcf_path
    assert engine.reference_fasta == reference_path
    assert engine.decompose_level == DecomposeLevel.CONSERVATIVE
    assert engine.left_shift is True
    assert engine.haploid_x is False


def test_normalize_variant():
    """Test variant normalization."""
    # Create dummy VCF and FASTA files for initialization
    with tempfile.NamedTemporaryFile(
        suffix=".vcf", delete=False
    ) as dummy_vcf_file, tempfile.NamedTemporaryFile(
        suffix=".fa", delete=False
    ) as dummy_fa_file:
        dummy_vcf_path = dummy_vcf_file.name
        dummy_fa_path = dummy_fa_file.name
        # Write minimal valid content to FASTA to avoid pysam errors
        # and ensure it's indexable by pysam.faidx
        dummy_fa_file.write(b">chr1\n")
        dummy_fa_file.write(b"A" * 60 + b"\n")  # Standard FASTA line length
        dummy_fa_file.write(b"C" * 60 + b"\n")
        dummy_fa_file.flush()
        pysam.faidx(dummy_fa_path)  # Create .fai index

    engine = PreprocessEngine(input_vcf=dummy_vcf_path, reference_fasta=dummy_fa_path)

    # Test normalization of SNV - should not change
    pos, ref, alt = engine.normalize_variant("chr1", 100, "A", "G")
    assert pos == 100
    assert ref == "A"
    assert alt == "G"

    # Test normalization with common prefix and suffix
    pos, ref, alt = engine.normalize_variant("chr1", 100, "ATCG", "ATTG")
    assert pos == 101
    assert ref == "TC"
    assert alt == "TT"

    # Test normalization with only common prefix
    pos, ref, alt = engine.normalize_variant("chr1", 100, "ATCG", "ATTT")
    assert pos == 102
    assert ref == "CG"
    assert alt == "TT"

    # Test normalization with only common suffix
    pos, ref, alt = engine.normalize_variant("chr1", 100, "ATCG", "TTCG")
    assert pos == 100
    assert ref == "AT"
    assert alt == "TT"

    # Clean up dummy files
    os.remove(dummy_vcf_path)
    os.remove(dummy_fa_path)
    if os.path.exists(dummy_fa_path + ".fai"):
        os.remove(dummy_fa_path + ".fai")


def test_process_file(reference_path, example_vcf_path, temp_output_path):
    """Test processing a VCF file."""
    engine = PreprocessEngine(
        input_vcf=example_vcf_path,
        reference_fasta=reference_path,
        output_vcf=temp_output_path,
    )

    output_file = engine.process()
    assert os.path.exists(output_file)

    # Check that the output file is a valid VCF
    vcf = pysam.VariantFile(output_file)
    assert vcf is not None

    # Check some basic stats
    assert engine.stats["total_variants"] > 0


def test_decompose_variants(reference_path, example_vcf_path, temp_output_path):
    """Test decomposition of multi-allelic variants."""
    # First count how many multi-allelic variants are in the example file
    vcf = pysam.VariantFile(example_vcf_path)
    multi_allelic_count = 0
    total_variants = 0

    for record in vcf:
        total_variants += 1
        if len(record.alts) > 1:
            multi_allelic_count += 1

    # Now process with decomposition
    engine = PreprocessEngine(
        input_vcf=example_vcf_path,
        reference_fasta=reference_path,
        output_vcf=temp_output_path,
        decompose_level=DecomposeLevel.AGGRESSIVE,
    )

    output_file = engine.process()

    # Check the output file
    vcf_out = pysam.VariantFile(output_file)
    processed_variants = 0
    for _ in vcf_out:
        processed_variants += 1

    # Should have more variants than the input if there were multi-allelic variants
    if multi_allelic_count > 0:
        assert processed_variants > total_variants
        assert engine.stats["decomposed_variants"] > 0


def test_left_shift_variants(
    reference_path,
):  # reference_path fixture will be used by PreprocessEngine
    """Test left-shifting of variants."""
    # Create a simple VCF with variants that need left-shifting
    with tempfile.NamedTemporaryFile(suffix=".vcf", delete=False) as f:
        f.write(
            b"""##fileformat=VCFv4.2
##reference=file://test.fa
##contig=<ID=chr1,length=1000000>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample1
chr1\t100\t.\tGTCG\tGTT\t.\tPASS\t.\tGT\t0/1
chr1\t200\t.\tCT\tCTT\t.\tPASS\t.\tGT\t0/1
chr1\t300\t.\tATCG\tAG\t.\tPASS\t.\tGT\t0/1
"""
        )
        input_vcf = f.name

    # Create a simple reference
    with tempfile.NamedTemporaryFile(suffix=".fa", delete=False) as f:
        f.write(
            b">chr1\n"
            + b"AGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCG\n"
            + b"TCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCG\n"
            + b"TCGTCGTCGTCGTCGCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTATCG\n"
            + b"ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\n"
        )
        f.flush()
        test_reference = f.name
        pysam.faidx(test_reference)  # Ensure .fai index is created

    # Process with left-shifting
    with tempfile.NamedTemporaryFile(suffix=".vcf.gz", delete=False) as f:
        output_vcf = f.name

    engine = PreprocessEngine(
        input_vcf=input_vcf,
        reference_fasta=test_reference,
        output_vcf=output_vcf,
        left_shift=True,
    )

    output_file = engine.process()

    # Check the output file
    vcf_out = pysam.VariantFile(output_file)
    for record in vcf_out:
        # Check if the record has been left-shifted
        if "LEFTSHIFTED" in record.info:
            assert engine.stats["left_shifted_variants"] > 0

    # Clean up
    os.remove(input_vcf)
    os.remove(test_reference)
    os.remove(output_vcf)
    if os.path.exists(output_vcf + ".tbi"):
        os.remove(output_vcf + ".tbi")


def test_haploid_x_handling(reference_path):  # reference_path fixture will be used
    """Test handling of haploid X chromosome."""
    # Create a simple VCF with X chromosome variants
    with tempfile.NamedTemporaryFile(suffix=".vcf", delete=False) as f:
        f.write(
            b"""##fileformat=VCFv4.2
##reference=file://test.fa
##contig=<ID=chrX,length=1000000>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample1
chrX\t100\t.\tA\tG\t.\tPASS\t.\tGT\t0/1
chrX\t200\t.\tC\tT\t.\tPASS\t.\tGT\t0/1
chrX\t300\t.\tG\tC\t.\tPASS\t.\tGT\t1/0
"""
        )
        input_vcf = f.name

    # Process with haploid X
    with tempfile.NamedTemporaryFile(suffix=".vcf.gz", delete=False) as f:
        output_vcf = f.name

    engine = PreprocessEngine(
        input_vcf=input_vcf,
        reference_fasta=reference_path,
        output_vcf=output_vcf,
        haploid_x=True,
    )

    output_file = engine.process()

    # Check the output file
    vcf_out = pysam.VariantFile(output_file)
    for record in vcf_out:
        # All genotypes should be homozygous
        for sample in record.samples:
            gt = record.samples[sample]["GT"]
            assert gt[0] == gt[1], "Genotype should be homozygous due to haploid X"

    assert engine.stats["haploid_adjusted_variants"] == 3

    # Clean up
    os.remove(input_vcf)
    os.remove(output_vcf)
    if os.path.exists(output_vcf + ".tbi"):
        os.remove(output_vcf + ".tbi")


def test_region_filtering(reference_path, example_vcf_path, temp_output_path):
    """Test filtering by region."""
    # First count total variants
    vcf = pysam.VariantFile(example_vcf_path)
    total_variants = sum(1 for _ in vcf)

    # Create a regions file with a subset of the genome
    with tempfile.NamedTemporaryFile(suffix=".bed", delete=False) as f:
        # Get first contig from the VCF
        vcf = pysam.VariantFile(example_vcf_path)
        first_contig = None
        for record in vcf:
            first_contig = record.contig
            break

        if first_contig:
            f.write(f"{first_contig}\t1\t1000000\n".encode())
        regions_file = f.name

    # Process with region filtering
    engine = PreprocessEngine(
        input_vcf=example_vcf_path,
        reference_fasta=reference_path,
        output_vcf=temp_output_path,
        regions=regions_file,
    )

    output_file = engine.process()

    # Count variants in the output
    vcf_out = pysam.VariantFile(output_file)
    filtered_variants = sum(1 for _ in vcf_out)

    # Should have filtered some variants if the region doesn't cover everything
    if first_contig:
        assert filtered_variants <= total_variants

    # Clean up
    os.remove(regions_file)


def test_pass_only_filtering(reference_path):  # reference_path fixture will be used
    """Test filtering for PASS variants only."""
    # Create a simple VCF with a mix of PASS and non-PASS variants
    with tempfile.NamedTemporaryFile(suffix=".vcf", delete=False) as f:
        f.write(
            b"""##fileformat=VCFv4.2
##reference=file://test.fa
##contig=<ID=chr1,length=1000000>
##FILTER=<ID=LowQual,Description="Low quality">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample1
chr1\t100\t.\tA\tG\t30\tPASS\t.\tGT\t0/1
chr1\t200\t.\tC\tT\t10\tLowQual\t.\tGT\t0/1
chr1\t300\t.\tG\tC\t40\tPASS\t.\tGT\t0/1
"""
        )
        input_vcf = f.name

    # Process with pass-only filtering
    with tempfile.NamedTemporaryFile(suffix=".vcf.gz", delete=False) as f:
        output_vcf = f.name

    engine = PreprocessEngine(
        input_vcf=input_vcf,
        reference_fasta=reference_path,
        output_vcf=output_vcf,
        pass_only=True,
    )

    output_file = engine.process()

    # Check the output file
    vcf_out = pysam.VariantFile(output_file)
    count = 0
    for record in vcf_out:
        count += 1
        assert "PASS" in record.filter or len(record.filter) == 0

    assert count == 2  # Only the PASS variants should be included
    assert engine.stats["filtered_variants"] == 1

    # Clean up
    os.remove(input_vcf)
    os.remove(output_vcf)
    if os.path.exists(output_vcf + ".tbi"):
        os.remove(output_vcf + ".tbi")
