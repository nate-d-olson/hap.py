#!/usr/bin/env python3
"""
Test cases for Python sequence utilities implementation.

This module tests the functionality of sequence_utils.py to ensure
it correctly performs sequence manipulations.
"""

import os
import tempfile

import pytest
from Haplo.sequence_utils import FastaReader, SequenceUtils


class TestSequenceUtils:
    """Test cases for SequenceUtils class."""

    def test_complement_sequence(self):
        """Test complementing DNA sequences."""
        # Test string input
        assert SequenceUtils.complement_sequence("ACGT") == "TGCA"
        assert SequenceUtils.complement_sequence("NNATGC") == "NNTACG"

        # Test bytes input
        assert SequenceUtils.complement_sequence(b"ACGT") == "TGCA"

        # Test empty input
        assert SequenceUtils.complement_sequence("") == ""

        # Test with non-standard bases
        assert SequenceUtils.complement_sequence("ACGTRYMKWSBDHVN") == "TGCAYRKMWSVHDB"

    def test_reverse_complement(self):
        """Test reverse complementing DNA sequences."""
        # Test string input
        assert SequenceUtils.reverse_complement("ACGT") == "ACGT"  # Palindromic
        assert SequenceUtils.reverse_complement("AATGC") == "GCATT"

        # Test bytes input
        assert SequenceUtils.reverse_complement(b"AATGC") == "GCATT"

        # Test empty input
        assert SequenceUtils.reverse_complement("") == ""

        # Test with non-standard bases
        assert SequenceUtils.reverse_complement("ACGTN") == "NACGT"

    def test_normalize_sequence(self):
        """Test sequence normalization."""
        # Test uppercase conversion
        assert SequenceUtils.normalize_sequence("acgt") == "ACGT"

        # Test bytes input
        assert SequenceUtils.normalize_sequence(b"acgt") == "ACGT"

        # Test invalid character replacement
        assert SequenceUtils.normalize_sequence("ACGTXYZ") == "ACGTNNN"

        # Test empty input
        assert SequenceUtils.normalize_sequence("") == ""


class TestFastaReader:
    """Test cases for FastaReader class."""

    @pytest.fixture
    def example_fasta(self):
        """Fixture to provide a temporary FASTA file."""
        # Create a temporary FASTA file for testing
        with tempfile.NamedTemporaryFile(suffix=".fa", delete=False) as temp_file:
            temp_file.write(b">chr1\n")
            temp_file.write(b"ACGTACGTACGTACGT\n")
            temp_file.write(b">chr2\n")
            temp_file.write(b"NNNNACGTACGTNNNN\n")
            temp_path = temp_file.name

        try:
            yield temp_path
        finally:
            # Clean up
            if os.path.exists(temp_path):
                os.unlink(temp_path)

            # Also clean up index file if created
            index_path = f"{temp_path}.fai"
            if os.path.exists(index_path):
                os.unlink(index_path)

    def test_init(self, example_fasta):
        """Test initialization with a FASTA file."""
        reader = FastaReader(example_fasta)
        # Initialization should not raise any exceptions
        assert reader is not None

    def test_get_sequence(self, example_fasta):
        """Test getting sequences from a FASTA file."""
        reader = FastaReader(example_fasta)

        # Get full sequence
        assert reader.get_sequence("chr1", 0, 16) == "ACGTACGTACGTACGT"

        # Get subsequence
        assert reader.get_sequence("chr1", 4, 8) == "ACGT"

        # Get from second chromosome
        assert reader.get_sequence("chr2", 4, 12) == "ACGTACGT"

    def test_get_invalid_sequence(self, example_fasta):
        """Test error handling for invalid sequence requests."""
        reader = FastaReader(example_fasta)

        # Invalid chromosome
        with pytest.raises(KeyError):
            reader.get_sequence("chr3", 0, 10)

        # Invalid coordinates
        with pytest.raises(IndexError):
            reader.get_sequence("chr1", -1, 10)

        with pytest.raises(IndexError):
            reader.get_sequence("chr1", 0, 100)

    def test_get_chromosomes(self, example_fasta):
        """Test getting the list of chromosomes."""
        reader = FastaReader(example_fasta)
        chroms = reader.get_chromosomes()

        assert len(chroms) == 2
        assert "chr1" in chroms
        assert "chr2" in chroms

    def test_get_chromosome_length(self, example_fasta):
        """Test getting chromosome lengths."""
        reader = FastaReader(example_fasta)

        assert reader.get_chromosome_length("chr1") == 16
        assert reader.get_chromosome_length("chr2") == 16

        # Invalid chromosome
        with pytest.raises(KeyError):
            reader.get_chromosome_length("chr3")

    def test_get_all_chromosome_lengths(self, example_fasta):
        """Test getting all chromosome lengths."""
        reader = FastaReader(example_fasta)
        lengths = reader.get_all_chromosome_lengths()

        assert len(lengths) == 2
        assert lengths["chr1"] == 16
        assert lengths["chr2"] == 16


if __name__ == "__main__":
    pytest.main(["-xvs", __file__])
