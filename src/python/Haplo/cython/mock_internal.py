#!/usr/bin/env python3
"""
Mock implementation of the Haplo.cython._internal module.

It can be used for testing Python code that depends on the C++ functionality
without requiring the actual C++ components to be built.

This is especially useful during the Python 3 migration process.
"""

import logging
from typing import Any, Dict, Union

# Configure module logger
logger = logging.getLogger(__name__)


def ensure_str(text: Union[str, bytes]) -> str:
    """
    Ensure the input is a string, converting from bytes if necessary.

    This function mimics the string handling at the C++/Python boundary
    for Python 3 compatibility.

    Args:
        text: Input string or bytes

    Returns:
        Unicode string
    """
    if isinstance(text, bytes):
        return text.decode("utf-8")
    return str(text)


def ensure_bytes(text: Union[str, bytes]) -> bytes:
    """
    Ensure the input is bytes, converting from string if necessary.

    Args:
        text: Input string or bytes

    Returns:
        Bytes
    """
    if isinstance(text, str):
        return text.encode("utf-8")
    return bytes(text)


class MockVariantRecord:
    """Mock implementation of a variant record"""

    def __init__(self, chrom: str, pos: int, ref: str, alt: str, qual: float = 0):
        self.chrom = chrom
        self.pos = pos
        self.ref = ref
        self.alt = alt
        self.qual = qual
        self.info = {}
        self.format = {}

    def __str__(self) -> str:
        return f"{self.chrom}:{self.pos} {self.ref}>{self.alt}"


class MockHaploCompare:
    """Mock implementation of the haplotype comparison functionality"""

    def __init__(self, reference_file: str = None):
        self.reference_file = reference_file
        self.truth_variants = []
        self.query_variants = []
        self.results = None

    def add_truth_variant(self, variant: Union[MockVariantRecord, Dict]) -> None:
        """Add a truth variant"""
        if isinstance(variant, dict):
            variant = MockVariantRecord(**variant)
        self.truth_variants.append(variant)

    def add_query_variant(self, variant: Union[MockVariantRecord, Dict]) -> None:
        """Add a query variant"""
        if isinstance(variant, dict):
            variant = MockVariantRecord(**variant)
        self.query_variants.append(variant)

    def compare(self) -> Dict[str, Any]:
        """Perform the comparison"""
        # Simple mock implementation
        tp = min(len(self.truth_variants), len(self.query_variants))
        fp = max(0, len(self.query_variants) - len(self.truth_variants))
        fn = max(0, len(self.truth_variants) - len(self.query_variants))

        self.results = {
            "total_truth": len(self.truth_variants),
            "total_query": len(self.query_variants),
            "true_positives": tp,
            "false_positives": fp,
            "false_negatives": fn,
        }
        return self.results


def complement_sequence(seq: str) -> str:
    """
    Return the complementary DNA sequence.

    Args:
        seq: A DNA sequence string

    Returns:
        Complementary DNA sequence
    """
    # Handle potential bytes input for Python 3 compatibility
    if isinstance(seq, bytes):
        seq = seq.decode("utf-8")

    # DNA complementation
    trans = str.maketrans("ACGTRYMKWSBDHVN", "TGCAYRKMWSVHDBN")
    return seq.upper().translate(trans)


def reverse_complement(seq: str) -> str:
    """
    Return the reverse complement of a DNA sequence.

    Args:
        seq: A DNA sequence string

    Returns:
        Reverse complemented DNA sequence
    """
    return complement_sequence(seq)[::-1]


def read_fasta_index(fasta_file: str) -> Dict[str, Dict[str, int]]:
    """
    Mock implementation of reading a FASTA index file.

    Args:
        fasta_file: Path to a FASTA file

    Returns:
        Dictionary with sequence information
    """
    logger.warning("Using mock FASTA index for %s", fasta_file)
    # Return a minimal mock index with some standard chromosomes
    return {
        "chr1": {"length": 248956422, "offset": 0},
        "chr2": {"length": 242193529, "offset": 0},
        "chr3": {"length": 198295559, "offset": 0},
        "chrX": {"length": 156040895, "offset": 0},
        "chrY": {"length": 57227415, "offset": 0},
    }


def get_reference_sequence(fasta_file: str, chrom: str, start: int, end: int) -> str:
    """
    Mock implementation of fetching reference sequence from a FASTA file.

    Args:
        fasta_file: Path to a FASTA file
        chrom: Chromosome name
        start: Start position (0-based)
        end: End position (exclusive)

    Returns:
        Reference sequence string
    """
    logger.warning("Using mock reference sequence for %s:%d-%d", chrom, start, end)
    # Return a mock sequence of appropriate length
    length = end - start
    # Generate a mock sequence with balanced nucleotide content
    bases = ["A", "C", "G", "T"]
    return "".join(bases[i % 4] for i in range(length))


def parse_vcf_file(vcf_file: str) -> Dict[str, Any]:
    """
    Mock implementation of parsing a VCF file.

    Args:
        vcf_file: Path to a VCF file

    Returns:
        Dictionary with VCF contents
    """
    logger.warning("Using mock VCF parser for %s", vcf_file)
    # Return a minimal mock VCF structure
    return {
        "records": [],
        "samples": ["SAMPLE"],
        "contigs": ["chr1", "chr2", "chr3", "chrX", "chrY"],
    }


def write_vcf_file(vcf_data: Dict[str, Any], output_file: str) -> None:
    """
    Mock implementation of writing a VCF file.

    Args:
        vcf_data: VCF data structure
        output_file: Path to output file
    """
    logger.warning("Using mock VCF writer for %s", output_file)
    # In a real implementation, this would write the VCF to a file
    # For mock purposes, just log the action
    record_count = len(vcf_data.get("records", []))
    logger.info("Would write %d records to %s", record_count, output_file)


from typing import Any, Dict, Union


class MockVariantRecord:
    """Mock implementation of a variant record"""

    def __init__(self, chrom: str, pos: int, ref: str, alt: str, qual: float = 0):
        self.chrom = chrom
        self.pos = pos
        self.ref = ref
        self.alt = alt
        self.qual = qual
        self.info = {}
        self.format = {}

    def __str__(self) -> str:
        return f"{self.chrom}:{self.pos} {self.ref}>{self.alt}"


class MockHaploCompare:
    """Mock implementation of the haplotype comparison functionality"""

    def __init__(self, reference_file: str = None):
        self.reference_file = reference_file
        self.truth_variants = []
        self.query_variants = []
        self.results = None

    def add_truth_variant(self, variant: Union[MockVariantRecord, Dict]) -> None:
        """Add a truth variant"""
        if isinstance(variant, dict):
            variant = MockVariantRecord(**variant)
        self.truth_variants.append(variant)

    def add_query_variant(self, variant: Union[MockVariantRecord, Dict]) -> None:
        """Add a query variant"""
        if isinstance(variant, dict):
            variant = MockVariantRecord(**variant)
        self.query_variants.append(variant)

    def compare(self) -> Dict[str, Any]:
        """Perform the comparison"""
        # Simple mock implementation
        self.results = {
            "total_truth": len(self.truth_variants),
            "total_query": len(self.query_variants),
            "tp": 0,
            "fp": 0,
            "fn": 0,
        }

        # Count exact matches as TP
        for qvar in self.query_variants:
            matched = False
            for tvar in self.truth_variants:
                if (
                    qvar.chrom == tvar.chrom
                    and qvar.pos == tvar.pos
                    and qvar.ref == tvar.ref
                    and qvar.alt == tvar.alt
                ):
                    matched = True
                    self.results["tp"] += 1
                    break
            if not matched:
                self.results["fp"] += 1

        # Count unmatched truth variants as FN
        for tvar in self.truth_variants:
            matched = False
            for qvar in self.query_variants:
                if (
                    qvar.chrom == tvar.chrom
                    and qvar.pos == tvar.pos
                    and qvar.ref == tvar.ref
                    and qvar.alt == tvar.alt
                ):
                    matched = True
                    break
            if not matched:
                self.results["fn"] += 1

        # Calculate metrics
        if self.results["tp"] + self.results["fn"] > 0:
            self.results["recall"] = self.results["tp"] / (
                self.results["tp"] + self.results["fn"]
            )
        else:
            self.results["recall"] = 0

        if self.results["tp"] + self.results["fp"] > 0:
            self.results["precision"] = self.results["tp"] / (
                self.results["tp"] + self.results["fp"]
            )
        else:
            self.results["precision"] = 0

        if self.results["precision"] + self.results["recall"] > 0:
            self.results["f1"] = (
                2
                * (self.results["precision"] * self.results["recall"])
                / (self.results["precision"] + self.results["recall"])
            )
        else:
            self.results["f1"] = 0

        return self.results


# Mock version functions
def get_version() -> str:
    """Get the hap.py version string (mock implementation)"""
    return "0.3.15-mock"


def get_build_time() -> str:
    """Get the hap.py build timestamp (mock implementation)"""
    return "2023-01-01 00:00:00"


def complement_sequence(seq: str) -> str:
    """Returns the complement of a DNA sequence (mock implementation)"""
    # Handle potential bytes input for Python 3 compatibility
    if isinstance(seq, bytes):
        seq = seq.decode("utf-8")

    # DNA complementation
    trans = str.maketrans("ACGTRYMKWSBDHVN", "TGCAYRKMWSVHDBN")
    return seq.upper().translate(trans)


def reverse_complement(seq: str) -> str:
    """Returns the reverse complement of a DNA sequence (mock implementation)"""
    # Handle potential bytes input for Python 3 compatibility
    if isinstance(seq, bytes):
        seq = seq.decode("utf-8")

    # Complement and reverse
    return complement_sequence(seq)[::-1]


def test_module() -> Dict[str, str]:
    """Test if the module is working properly (mock implementation)"""
    return {"version": get_version(), "build_time": get_build_time()}
