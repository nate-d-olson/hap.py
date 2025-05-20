import os


def test_string_handling_module():
    """Test the string handling utilities."""
    from Haplo.string_handling import ensure_bytes, ensure_str, ensure_text_io

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
    # Store original environment value
    original_value = os.environ.get("HAPLO_USE_MOCK", None)

    try:
        # Set environment variable to use mocks
        os.environ["HAPLO_USE_MOCK"] = "1"

        # Import the Cython module package
        import Haplo.cython

        # Verify we're using mocks
        assert Haplo.cython.USING_MOCK is True

        # Test complement_sequence function
        seq = "ACGTACGT"
        comp_seq = Haplo.cython.complement_sequence(seq)
        assert comp_seq == "TGCATGCA"

        # Test with bytes input (Python 3 compatibility test)
        bytes_seq = b"ACGT"
        str_result = Haplo.cython.complement_sequence(bytes_seq)
        assert isinstance(str_result, str)
        assert str_result == "TGCA"

    finally:
        # Restore original environment
        if original_value is None:
            del os.environ["HAPLO_USE_MOCK"]
        else:
            os.environ["HAPLO_USE_MOCK"] = original_value


def test_mock_variant_classes():
    """Test the mock variant record and comparison classes."""
    # Store original environment value
    original_value = os.environ.get("HAPLO_USE_MOCK", None)

    try:
        # Set environment variable to use mocks
        os.environ["HAPLO_USE_MOCK"] = "1"

        # Import the Cython module package
        import Haplo.cython

        # Test record class
        record = Haplo.cython.MockVariantRecord(
            chrom="chr1", pos=100, ref="A", alt="T", qual=30
        )
        assert record.chrom == "chr1"
        assert record.pos == 100
        assert record.ref == "A"
        assert record.alt == "T"
        assert str(record) == "chr1:100 A>T"

        # Test mock comparison
        comparator = Haplo.cython.MockHaploCompare()
        comparator.add_truth_variant(record)
        comparator.add_query_variant(record)
        results = comparator.compare()

        # Check the results - ensure they have expected keys
        assert "total_truth" in results
        assert "total_query" in results
        assert results["total_truth"] == 1
        assert results["total_query"] == 1

        # Check either true_positives or true_positives are present
        # (allowing either naming convention)
        tp_key = next(
            (k for k in results if k.lower().replace("_", "") == "truepositives"), None
        )
        assert tp_key is not None, "No true positives key found in results"
        assert results[tp_key] == 1

    finally:
        # Restore original environment
        if original_value is None:
            del os.environ["HAPLO_USE_MOCK"]
        else:
            os.environ["HAPLO_USE_MOCK"] = original_value
