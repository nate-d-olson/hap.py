"""
Unit tests for the fastaContigLengths and fastaSampleRegions utilities.
"""
import os
import re

import pytest

from happy.Tools.fastasize import fastaContigLengths, fastaSampleRegions


def test_fasta_contig_lengths(tmp_path):
    # Create a dummy FASTA file and its .fai index
    fa = tmp_path / "test.fa"
    fa.write_text(">chr1\n" "ATGC\n" ">chr2\n" "AAA\n")
    fai = tmp_path / "test.fa.fai"
    # Minimal .fai with contig names and lengths
    fai.write_text("chr1\t4\n" "chr2\t3\n")
    # Verify contig lengths are read correctly
    result = fastaContigLengths(str(fa))
    assert result == {"chr1": 4, "chr2": 3}


def test_fasta_contig_lengths_missing_index(tmp_path):
    # FASTA without .fai should raise
    fa = tmp_path / "noindex.fa"
    fa.write_text(">chrX\nATGC\n")
    with pytest.raises(Exception):
        fastaContigLengths(str(fa))


def test_fasta_sample_regions_empty(tmp_path):
    # n_regions=0 should return empty list
    fa = tmp_path / "empty.fa"
    fa.write_text(">chr1\n" "ATGCATGC\n")
    fai = tmp_path / "empty.fa.fai"
    fai.write_text("chr1\t8\n")
    regions = fastaSampleRegions(str(fa), n_regions=0)
    assert regions == []


def test_fasta_sample_regions_one(tmp_path):
    # n_regions=1 should return exactly one region of specified length
    fa = tmp_path / "one.fa"
    seq = "A" * 12
    fa.write_text(f">chrA\n{seq}\n")
    fai = tmp_path / "one.fa.fai"
    fai.write_text("chrA\t12\n")
    regions = fastaSampleRegions(str(fa), n_regions=1, region_length=5)
    assert isinstance(regions, list) and len(regions) == 1
    region = regions[0]
    # Format should be 'chrom:start-end'
    assert re.match(r"^chrA:\d+-\d+$", region)
    # Verify region length == requested length
    coords = region.split(":")[1]
    start, end = map(int, coords.split("-"))
    assert (end - start) == 5
