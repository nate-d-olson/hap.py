import pytest

from hap_py.haplo.quantify_models import StratificationRegion


def test_stratification_region_with_bed_file(sample_bed_file):
    """Construct StratificationRegion with a bed file path."""
    region = StratificationRegion(name="test", bed_file=sample_bed_file)

    assert region.name == "test"
    assert region.bed_file == sample_bed_file
    assert region.filter_expression is None


def test_stratification_region_with_filter_expression():
    """Construct StratificationRegion with a filter expression."""
    expr = "QUAL > 30"
    region = StratificationRegion(name="filter_region", filter_expression=expr)

    assert region.name == "filter_region"
    assert region.bed_file is None
    assert region.filter_expression == expr


def test_stratification_region_requires_definition():
    """ValueError raised when neither bed_file nor filter_expression provided."""
    with pytest.raises(ValueError):
        StratificationRegion(name="invalid")
