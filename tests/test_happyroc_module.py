import tempfile

import pytest

# skip if pandas not installed
pytest.importorskip("pandas", reason="pandas is required for happyroc tests")
import pandas as pd
from Haplo.happyroc import RESULT_ALLCOLUMNS, roc


def test_happyroc_roc_all_key(tmp_path):
    # Create a simple ROC-like tab-delimited file
    header_cols = [
        "Type",
        "Subtype",
        "Subset",
        "Filter",
        "Genotype",
        "QQ",
        "TRUTH.TOTAL",
    ]
    header = "\t".join(header_cols)
    lines = [
        header,
        "SNP\t*\t*\tALL\t*\t1\t10",
        "INDEL\t*\t*\tPASS\t*\t2\t20",
    ]
    roc_file = tmp_path / "roc.txt"
    roc_file.write_text("\n".join(lines))
    # Run ROC parser
    result = roc(str(roc_file), output_path=None)
    # 'all' key should always be present
    assert "all" in result
    df_all = result["all"]
    # DataFrame must contain the expected columns
    for col in RESULT_ALLCOLUMNS:
        assert col in df_all.columns, f"Missing column {col} in ROC DataFrame"


@pytest.mark.parametrize(
    "filter_val, expected_key",
    [
        (None, "all"),
        ("PASS", "Locations.SNP.PASS"),
    ],
)
def test_happyroc_filter_handling(tmp_path, filter_val, expected_key):
    # Verify that specifying a filter restricts output keys
    header = "Type\tSubtype\tSubset\tFilter\tGenotype\tQQ\tTRUTH.TOTAL"
    records = [
        "SNP\t*\t*\tPASS\t*\t1\t5",
        "SNP\t*\t*\tALL\t*\t2\t10",
    ]
    f = tmp_path / "r2.txt"
    f.write_text("\n".join([header] + records))
    result = roc(str(f), output_path=None, filter_handling=filter_val)
    # Check that expected key is in result when filter matches
    assert expected_key in result
