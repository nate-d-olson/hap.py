import pytest

from hap_py.tools.fastasize import calculateLength


@pytest.mark.unit
def test_fastasize_calculation():
    """Ensure calculateLength returns expected result."""
    locations = "chrMT chrY:1-10"
    contigs = "{'chrY': 59373566, 'chrM': 16571}"
    assert calculateLength(contigs, locations) == 10
