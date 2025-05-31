"""
Pure-Python implementations of sequence utility functions.
"""
# Complement mapping
_COMP_MAP = {
    "A": "T",
    "C": "G",
    "G": "C",
    "T": "A",
    "a": "t",
    "c": "g",
    "g": "c",
    "t": "a",
    "N": "N",
    "n": "n",
}


def complement_sequence(seq):
    """Return the complement of a DNA sequence (string or bytes)."""
    is_bytes = isinstance(seq, (bytes, bytearray))
    s = seq.decode("ascii") if is_bytes else seq
    comp = "".join(_COMP_MAP.get(ch, ch) for ch in s)
    return comp.encode("ascii") if is_bytes else comp


def reverse_complement(seq):
    """Return the reverse complement of a DNA sequence (string or bytes)."""
    comp = complement_sequence(seq)
    is_bytes = isinstance(comp, (bytes, bytearray))
    s = comp if is_bytes else comp
    rev = s[::-1]
    return rev


def process_sequence(seq):
    """Alias for reverse_complement -- normalize and reverse sequence."""
    return reverse_complement(seq)
