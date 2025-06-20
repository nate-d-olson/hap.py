import pytest

from hap_py.tools.fastasize import (
    calculateLength,
    fastaNonNContigLengths,
    fastaSampleRegions,
)


@pytest.mark.unit
def test_fastasize_calculation():
    """Ensure calculateLength returns expected result."""
    locations = "chrMT chrY:1-10"
    contigs = "{'chrY': 59373566, 'chrM': 16571}"
    assert calculateLength(contigs, locations) == 10


@pytest.mark.unit
def test_fasta_non_n_lengths():
    lengths = fastaNonNContigLengths("tests/data/common/test.fa")
    assert lengths["chrQ"] == 8
    assert lengths["all"] == 8


def test_fasta_sample_regions(monkeypatch):
    # Force deterministic random selection
    monkeypatch.setattr("random.randint", lambda a, b: 0)
    regions = fastaSampleRegions(
        "tests/data/common/test.fa", n_regions=1, region_length=2
    )
    assert regions == ["chrQ:0-2"]
