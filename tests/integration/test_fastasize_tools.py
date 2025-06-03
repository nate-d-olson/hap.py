import pytest

from happy.Tools.fastasize import calculateLength


@pytest.mark.unit
def test_calculate_length_simple():
    locations = "chrMT chrY:1-10"
    fastacontiglengths = {"chrY": 59373566, "chrM": 16571}
    length = calculateLength(fastacontiglengths, locations)
    assert length == 10, f"Expected length 10, got {length}"


@pytest.mark.unit
def test_calculate_length_full_contig():
    locations = "chrM"
    fastacontiglengths = {"chrM": 500}
    length = calculateLength(fastacontiglengths, locations)
    assert length == 500, f"Expected length 500, got {length}"
