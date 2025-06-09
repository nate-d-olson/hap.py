"""
Mock implementation of cpp_internal for Python 3 testing without C++ components

This module provides fallback functionality when the Cython module
cannot be imported or compiled.
"""

import logging
import warnings

import numpy as np

# Log warning about using mock implementation
warnings.warn("Using mock implementation of cpp_internal", stacklevel=2)
logging.warning("Using mock implementation of cpp_internal")


def test_string_handling():
    """Mock test for string handling"""
    return "Mock implementation working"


def test_basic_functionality():
    """Mock test for basic functionality"""
    return "Mock basic functionality working"


def complement_sequence(seq):
    """
    Return the complementary DNA sequence.

    Args:
        seq: A DNA sequence string or bytes

    Returns:
        Complementary DNA sequence
    """
    # Handle potential bytes input for Python 3 compatibility
    if isinstance(seq, bytes):
        seq = seq.decode("utf-8")

    # DNA complementation
    trans = str.maketrans("ACGTRYMKWSBDHVN", "TGCAYRKMWSVHDBN")
    return seq.upper().translate(trans)


def reverse_complement(seq):
    """
    Return the reverse complement of a DNA sequence.

    Args:
        seq: A DNA sequence string or bytes

    Returns:
        Reverse complemented DNA sequence
    """
    # Handle potential bytes input for Python 3 compatibility
    if isinstance(seq, bytes):
        seq = seq.decode("utf-8")

    return complement_sequence(seq)[::-1]


class VariantProcessor:
    """Mock implementation of the VariantProcessor Cython class."""

    def __init__(self):
        """Initialize an empty variant processor."""
        self.variants = []

    def add_variant(self, variant):
        """
        Add a variant to the processor.

        Args:
            variant: Variant object with chrom, pos, ref, alt attributes
        """
        self.variants.append(variant)

    def get_variant_chrom(self, idx):
        """
        Get chromosome name for a variant.

        Args:
            idx: Index of the variant

        Returns:
            str: Chromosome name
        """
        if idx >= len(self.variants):
            raise IndexError(f"Index {idx} out of range")
        return self.variants[idx].chrom

    def process_variants(self, threads=1):
        """
        Process variants with mock implementation.

        Args:
            threads: Number of threads to use (ignored in mock)

        Returns:
            list: Processed variant data
        """
        results = []
        for variant in self.variants:
            # Create a simple result dictionary
            result = {
                "chrom": variant.chrom,
                "position": variant.pos,
                "ref": variant.ref,
                "alt": variant.alt,
                "processed": True,
            }
            results.append(result)
        return results


def compute_roc_points(tp, fp, fn, tn=None, resolution=100):
    """
    Compute ROC curve points.

    Args:
        tp: True positives array or list
        fp: False positives array or list
        fn: False negatives array or list
        tn: True negatives array or list (optional)
        resolution: Number of points in the ROC curve

    Returns:
        tuple: (recalls, precisions, f1_scores, specificities)
    """
    # Convert inputs to numpy arrays if they aren't already
    tp = np.asarray(tp)
    fp = np.asarray(fp)
    fn = np.asarray(fn)
    if tn is not None:
        tn = np.asarray(tn)

    # Calculate metrics
    with np.errstate(divide="ignore", invalid="ignore"):
        recalls = tp / (tp + fn)
        precisions = tp / (tp + fp)

    # Handle divide by zero
    recalls = np.nan_to_num(recalls)
    precisions = np.nan_to_num(precisions)

    # Calculate F1 scores
    f1_scores = 2 * precisions * recalls / (precisions + recalls)
    f1_scores = np.nan_to_num(f1_scores)

    # Calculate specificities if tn is provided
    if tn is not None:
        specificities = tn / (tn + fp)
        specificities = np.nan_to_num(specificities)
    else:
        specificities = np.zeros_like(recalls)

    return recalls, precisions, f1_scores, specificities


# Add mock implementations of the functions in the Cython module
# TODO: Analyze the original module and add appropriate mock functions
