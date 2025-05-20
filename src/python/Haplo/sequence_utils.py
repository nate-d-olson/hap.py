"""
Pure-Python sequence handling utilities for complement and reverse complement.
"""


def complement_sequence(seq_input):
    """
    Return the complement of a DNA sequence.

    Args:
        seq_input: String or bytes sequence to complement

    Returns:
        Same type as input: Complemented sequence
    """
    is_bytes = isinstance(seq_input, bytes)
    if is_bytes:
        seq = seq_input.decode("ascii")
    else:
        seq = seq_input

    comp_dict = {
        "A": "T",
        "T": "A",
        "G": "C",
        "C": "G",
        "N": "N",
        "a": "t",
        "t": "a",
        "g": "c",
        "c": "g",
        "n": "n",
    }
    result = "".join(comp_dict.get(base, base) for base in seq)

    if is_bytes:
        return result.encode("ascii")
    return result


def reverse_complement(seq_input):
    """
    Return the reverse complement of a DNA sequence.

    Args:
        seq_input: String or bytes sequence to reverse complement

    Returns:
        Same type as input: Reverse complemented sequence
    """
    comp_seq = complement_sequence(seq_input)
    if isinstance(comp_seq, bytes):
        return comp_seq[::-1]
    return comp_seq[::-1]


def process_sequence(seq_input):
    """
    Process a DNA sequence, returning a string or bytes.

    Args:
        seq_input: Input sequence (string or bytes)

    Returns:
        str or bytes: Processed sequence
    """
    result = reverse_complement(seq_input)
    if isinstance(seq_input, str):
        return result.decode("ascii") if isinstance(result, bytes) else result
    return result
