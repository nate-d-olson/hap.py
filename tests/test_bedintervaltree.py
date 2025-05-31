import gzip
import os
import tempfile

import pytest
from Tools.bedintervaltree import BedIntervalTree


@pytest.fixture
def sample_bed(tmp_path):
    content = """
chr1	0	10	label1
chr1	5	15	label2
"""
    bed_file = tmp_path / "test.bed"
    bed_file.write_text(content.strip() + "\n")
    return str(bed_file)


def test_add_and_counts(sample_bed):
    tree = BedIntervalTree()
    tree.addFromBed(sample_bed)
    # Two intervals added
    assert tree.count() == 2
    assert tree.count("label1") == 1
    assert tree.count("label2") == 1
    # Base counts: each interval covers 10 bases
    total_bases = tree.countbases()
    assert total_bases == 20
    # Chromosome-specific countbases
    bases_chr = tree.countbases(chrom="chr1", start=1, end=20)
    assert bases_chr == 20


def test_intersect(sample_bed):
    tree = BedIntervalTree()
    tree.addFromBed(sample_bed)
    # Query a region overlapping both intervals
    ivals = tree.intersect("chr1", start=7, end=8)
    # Expect two intervals overlapping position 7
    assert len(ivals) == 2
    labels = {tuple(iv.value)[0] for iv in ivals}
    assert labels == {"label1", "label2"}


def test_gzipped_bed(tmp_path, sample_bed):
    # Create gzipped bed file
    gz_path = tmp_path / "test2.bed.gz"
    with gzip.open(gz_path, "wt") as gz:
        gz.write(open(sample_bed).read())
    tree = BedIntervalTree()
    tree.addFromBed(str(gz_path))
    assert tree.count() == 2
