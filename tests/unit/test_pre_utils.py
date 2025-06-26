import pytest

from src.hap_py.pre import hasChrPrefix


@pytest.mark.unit
def test_hasChrPrefix_empty_list():
    # No chromosomes provided: undecided
    assert hasChrPrefix([]) is None


@pytest.mark.unit
def test_hasChrPrefix_prefixed_only():
    # Only prefixed names should return True
    assert hasChrPrefix(["chr1", "chr2", "chrX", "chrY"]) is True


@pytest.mark.unit
def test_hasChrPrefix_non_prefixed_only():
    # Only non-prefixed names should return False
    assert hasChrPrefix(["1", "2", "X", "Y"]) is False


@pytest.mark.unit
def test_hasChrPrefix_mixed_equal_counts():
    # Equal counts of prefixed and non-prefixed: undecided
    assert hasChrPrefix(["1", "chr1"]) is None


@pytest.mark.unit
def test_hasChrPrefix_mixed_unequal_counts():
    # More prefixed than non-prefixed: True
    assert hasChrPrefix(["chr1", "chr2", "1"]) is True
    # More non-prefixed than prefixed: False
    assert hasChrPrefix(["2", "3", "chr2"]) is False


@pytest.mark.unit
def test_hasChrPrefix_accepts_set_input():
    # Iterable types like set should be accepted
    assert hasChrPrefix(set(["1", "2", "3"])) is False
    # Iterable types like set should be accepted
    assert hasChrPrefix(set(["1", "2", "3"])) is False
