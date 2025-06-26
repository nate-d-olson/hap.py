from pathlib import Path

import pytest

pytest.importorskip("bx")

from hap_py.tools.bedintervaltree import BedIntervalTree


def test_intersect_and_counts(tmp_path: Path) -> None:
    bed = tmp_path / "ints.bed"
    bed.write_text(
        "chr1\t0\t10\nchr1\t20\t30\nchr2\t5\t15\n",
        encoding="utf-8",
    )
    tree = BedIntervalTree()
    tree.addFromBed(str(bed), label="test")

    overlaps = tree.intersect("chr1", 5, 25)
    assert len(overlaps) == 2
    assert [iv.value for iv in overlaps] == [["test"], ["test"]]

    assert tree.countbases("chr1", 1, 31) == 20
    assert tree.count(label="test") == 3
