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


from typing import Union


def complement_sequence(seq: Union[str, bytes, bytearray]) -> Union[str, bytes]:
    """Return the complement of a DNA sequence (string or bytes)."""
    is_bytes = isinstance(seq, (bytes, bytearray))
    s: str = seq.decode("ascii") if is_bytes else seq  # type: ignore
    comp = "".join(_COMP_MAP.get(ch, ch) for ch in s)
    return comp.encode("ascii") if is_bytes else comp


def reverse_complement(seq: Union[str, bytes, bytearray]) -> Union[str, bytes]:
    """Return the reverse complement of a DNA sequence (string or bytes)."""
    comp = complement_sequence(seq)
    is_bytes = isinstance(comp, (bytes, bytearray))
    # comp is str or bytes
    if is_bytes:
        rev = comp[::-1]  # type: ignore
    else:
        rev = comp[::-1]  # type: ignore
    return rev


def process_sequence(seq: Union[str, bytes, bytearray]) -> Union[str, bytes]:
    """Alias for reverse_complement -- normalize and reverse sequence."""
    return reverse_complement(seq)
