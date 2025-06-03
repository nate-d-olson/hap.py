"""
Unit tests for the Tools.vcfextract module: field conversion, INFO parsing, and header extraction.
"""

import gzip
import os

import pytest

from happy.Tools.vcfextract import extract_header, field, getInfo


def test_field_simple():
    assert field("123") == 123
    assert field("3.14") == 3.14
    assert field("abc") == "abc"


def test_field_list():
    assert field("1,2,3") == [1, 2, 3]
    assert field("a,b,3") == ["a", "b", 3]


def test_getInfo():
    info_str = "DP=100;AF=0.5;FLAG"
    info = getInfo(info_str)
    assert info["DP"] == 100
    assert pytest.approx(info["AF"], rel=1e-6) == 0.5
    assert info["FLAG"] is True


def test_extract_header(tmp_path):
    vcf = tmp_path / "test.vcf"
    lines = [
        "##fileformat=VCFv4.2",
        '##INFO=<ID=DP,Number=1,Type=Integer,Description="Read depth">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        "#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT S1",
        "1\t100\t.\tA\tT\t.\t.\tDP=10\tGT\t0/1",
    ]
    vcf.write_text("\n".join(lines) + "\n")
    hdr = extract_header(
        str(vcf), extract_columns=True, extract_info=True, extract_formats=True
    )
    assert "columns" in hdr
    assert hdr["columns"][0] == "CHROM"
    assert "DP" in hdr["info"]
    assert hdr["info"]["DP"]["type"] == "Integer"
    assert "GT" in hdr["formats"]
    assert hdr["formats"]["GT"]["number"] == "1"
