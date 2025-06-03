"""
Unit tests for the pure-Python VariantProcessor in Haplo.variant_processor.
"""

import pytest
from Haplo.variant_processor import VariantProcessor, create_standard_processor


class DummyVariant:
    def __init__(self, chrom, pos, ref, alt):
        self.chrom = chrom
        self.pos = pos
        self.ref = ref
        self.alt = alt


def test_add_and_get_chromosome():
    vp = VariantProcessor()
    var = DummyVariant("chr1", 100, "A", "T")
    vp.add_variant(var)
    assert vp.get_variant_chrom(0) == "chr1"
    with pytest.raises(IndexError):
        vp.get_variant_chrom(1)


def test_process_variants_default():
    vp = create_standard_processor()
    # add two dummy variants
    v1 = DummyVariant("chrX", 5, "G", "C")
    v2 = DummyVariant("chrY", 10, "T", "TA")
    vp.add_variant(v1)
    vp.add_variant(v2)
    results = vp.process_variants()
    assert isinstance(results, list)
    assert len(results) == 2
    for res, orig in zip(results, [v1, v2]):
        assert res["chrom"] == orig.chrom
        assert res["position"] == orig.pos
        assert res["ref"] == orig.ref
        assert res["alt"] == orig.alt
        assert res["processed"] is True
