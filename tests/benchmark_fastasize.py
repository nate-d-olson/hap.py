"""
Microbenchmarks for Tools.fastasize module.
"""
import os
import tempfile

import pytest

from happy.Tools.fastasize import calculateLength, fastaContigLengths


@pytest.fixture(scope="module")
def fasta_and_index(tmp_path_factory):
    # Create a temporary FASTA and .fai index
    d = tmp_path_factory.mktemp("fasta")
    fa = d / "test.fa"
    content = ">chrA\nAAAAA\n>chrB\nCCCCCCCCC\n>chrC\nGGGG"  # lengths 5,9,4
    fa.write_text(content + "\n")
    # Create corresponding .fai file
    fai = d / "test.fa.fai"
    # fai columns: name, length, offset, linebases, linewidth
    lines = [
        f"chrA\t5\t0\t6\t6",
        f"chrB\t9\t0\t9\t9",
        f"chrC\t4\t0\t4\t4",
    ]
    fai.write_text("\n".join(lines) + "\n")
    return str(fa)


@pytest.mark.parametrize("n_regions", [1, 10, 100])
def test_fasta_operations_benchmark(benchmark, fasta_and_index, n_regions):
    # Benchmark contig length extraction and region length calculation
    fasta = fasta_and_index
    regions = ["chrA:1-3"] * n_regions

    def bench_all():
        lengths = fastaContigLengths(fasta)
        calculateLength(lengths, regions)

    benchmark(bench_all)
