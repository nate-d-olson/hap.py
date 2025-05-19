# cython: language_level=3
# distutils: language=c++

# Removed __future__ imports; targeting Python 3 only

from libc.stdlib cimport free, malloc
from libc.string cimport memcpy

import numpy as np

cimport numpy as np


# Fast sequence handling with proper Python 3 bytes/string conversion
cpdef bytes complement_sequence(seq_input):
    """
    Return the complement of a DNA sequence.

    Args:
        seq_input: String or bytes sequence to complement

    Returns:
        bytes: Complemented sequence as bytes
    """
    # Convert to bytes if input is string
    cdef bytes seq
    if isinstance(seq_input, str):
        seq = seq_input.encode('ascii')
    else:
        seq = seq_input  # Assume bytes

    cdef:
        size_t i, n = len(seq)
        char* c_seq = seq  # This works in Python 3 because bytes maps to char*
        char* result = <char*>malloc(n + 1)

    if result == NULL:
        raise MemoryError("Failed to allocate memory for sequence complement")

    try:
        for i in range(n):
            if c_seq[i] == b'A'[0]:
                result[i] = b'T'[0]
            elif c_seq[i] == b'T'[0]:
                result[i] = b'A'[0]
            elif c_seq[i] == b'G'[0]:
                result[i] = b'C'[0]
            elif c_seq[i] == b'C'[0]:
                result[i] = b'G'[0]
            else:
                result[i] = c_seq[i]

        result[n] = 0  # Null terminator
        return result[:n]  # Convert back to Python bytes
    finally:
        free(result)

cpdef bytes reverse_complement(seq_input):
    """
    Return the reverse complement of a DNA sequence.

    Args:
        seq_input: String or bytes sequence to reverse complement

    Returns:
        bytes: Reverse complemented sequence as bytes
    """
    # First get complement
    cdef bytes comp_seq = complement_sequence(seq_input)

    # Then reverse it
    return comp_seq[::-1]

# Convert to proper Python string if needed
def process_sequence(seq_input):
    """
    Process a DNA sequence, returning a string.

    Args:
        seq_input: Input sequence (string or bytes)

    Returns:
        str: Processed sequence as a Python string
    """
    cdef bytes result = reverse_complement(seq_input)

    # Convert back to string if input was string
    if isinstance(seq_input, str):
        return result.decode('ascii')
    return result
