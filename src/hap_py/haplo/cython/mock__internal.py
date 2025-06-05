"""
Mock implementation of _internal for Python 3 testing without C++ components

This module provides fallback functionality when the Cython module
cannot be imported or compiled.
"""

import logging
import warnings

# Log warning about using mock implementation
warnings.warn("Using mock implementation of _internal", stacklevel=2)
logging.warning("Using mock implementation of _internal")


def test_string_handling():
    """Mock test for string handling"""
    return "Mock implementation working"


def test_basic_functionality():
    """Mock test for basic functionality"""
    return "Mock basic functionality working"


def complement_sequence(sequence: str) -> str:
    """Return the DNA complement for a sequence."""
    trans = str.maketrans("ACGTacgt", "TGCAtgca")
    return sequence.translate(trans)


def reverse_complement(sequence: str) -> str:
    """Return the reverse complement of a sequence."""
    return complement_sequence(sequence)[::-1]

def get_version() -> str:
    """Return a mock hap.py version string."""
    return "0.0.mock"


def get_build_time() -> str:
    """Return a mock build timestamp."""
    return "1970-01-01T00:00:00"


def is_available() -> bool:
    """Indicate that the mock implementation is available."""
    return True


def test_module() -> dict:
    """Test if the module is working properly."""
    return {"version": get_version(), "build_time": get_build_time()}

