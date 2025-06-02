"""
Microbenchmarks for Tools.vcfextract module.
"""
import tempfile
import textwrap

import pytest

from happy.Tools.vcfextract import extract_header, getInfo


@pytest.fixture(scope="module")
def small_vcf(tmp_path_factory):
    # Create a tiny VCF file with headers and one record
    d = tmp_path_factory.mktemp("vcf")
    v = d / "small.vcf"
    content = textwrap.dedent(
        """
        ##fileformat=VCFv4.2
        ##INFO=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
        ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
        #CHROM POS ID REF ALT QUAL FILTER INFO FORMAT S1
        1 100 . A T . PASS DP=100 GT 0/1
        """
    )
    # normalize spaces and write
    lines = [
        "#" + line if line.startswith("#CHROM") else line
        for line in content.splitlines()
    ]
    v.write_text("\n".join(lines) + "\n")
    return str(v)


@pytest.mark.parametrize("n_calls", [100, 1000, 5000])
def test_getInfo_benchmark(benchmark, n_calls):
    sample = "DP=100;AF=0.123;FLAG"

    # Benchmark parsing INFO fields repeatedly
    def fn():
        for _ in range(n_calls):
            _ = getInfo(sample)

    benchmark(fn)


def test_extract_header_benchmark(benchmark, small_vcf):
    # Benchmark header extraction once
    benchmark(extract_header, small_vcf, True, True, True)
