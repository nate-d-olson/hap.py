"""
Microbenchmarks for the pure-Python VariantProcessor implementation.
"""
import pytest
from Haplo.variant_processor import VariantProcessor


class DummyVariant:
    """Lightweight variant object for benchmarking."""

    __slots__ = ("chrom", "pos", "ref", "alt")

    def __init__(self, chrom, pos, ref, alt):
        self.chrom = chrom
        self.pos = pos
        self.ref = ref
        self.alt = alt


def make_processor(n):
    """Initialize a VariantProcessor with n dummy variants."""
    vp = VariantProcessor()
    for i in range(n):
        # cycle through a simple substitution variant
        vp.add_variant(DummyVariant("chr1", i + 1, "A", "T"))
    return vp


@pytest.mark.parametrize("n", [100, 1000, 10000])
def test_process_variants_benchmark(benchmark, n):
    """Benchmark processing of n variants."""
    vp = make_processor(n)
    # benchmark the process_variants method
    benchmark(vp.process_variants)
