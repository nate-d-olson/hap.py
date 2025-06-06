"""
Tests for string handling and mock Cython implementations in Python 3.
"""


def test_string_handling_module():
    """Test the string handling utilities."""
    from hap_py.haplo.string_handling import (
        ensure_bytes,
        ensure_str,
        ensure_text_io,
    )

    # Test ensure_str
    assert ensure_str(b"test") == "test"
    assert ensure_str("test") == "test"
    assert ensure_str(None) is None

    # Test ensure_bytes
    assert ensure_bytes("test") == b"test"
    assert ensure_bytes(b"test") == b"test"
    assert ensure_bytes(None) is None

    # Test ensure_text_io
    assert isinstance(ensure_text_io("test", "r"), str)
    assert isinstance(ensure_text_io("test", "rb"), bytes)
    assert isinstance(ensure_text_io(b"test", "r"), str)
    assert isinstance(ensure_text_io(b"test", "rb"), bytes)


def test_cython_mock_import():
    """Test using mock Cython implementations."""
    from hap_py.haplo import cython_mock as cython_module

    seq = "ACGTACGT"
    comp_seq = cython_module.complement_sequence(seq)
    assert comp_seq == "TGCATGCA"

    bytes_seq = b"ACGT"
    str_result = cython_module.complement_sequence(bytes_seq)
    assert isinstance(str_result, bytes)
    assert str_result == b"TGCA"
