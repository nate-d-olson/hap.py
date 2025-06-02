import pytest

from happy.Haplo.sequence_utils import (
    complement_sequence,
    process_sequence,
    reverse_complement,
)


@pytest.mark.parametrize(
    "seq, expected_comp",
    [
        ("ACGTNacgtn", "TGCANtgcan"),
        ("", ""),
        ("AAAaaa", "TTTttt"),
        ("XYZ", "XYZ"),
    ],
)
def test_complement_sequence_str(seq, expected_comp):
    assert complement_sequence(seq) == expected_comp


@pytest.mark.parametrize(
    "seq_bytes, expected_comp_bytes",
    [
        (b"ACGTN", b"TGCAN"),
        (b"", b""),
        (b"AaNn", b"TtNn"),
    ],
)
def test_complement_sequence_bytes(seq_bytes, expected_comp_bytes):
    assert complement_sequence(seq_bytes) == expected_comp_bytes


@pytest.mark.parametrize(
    "seq, expected_revcomp",
    [
        ("ACGT", "ACGT"[::-1].translate(str.maketrans("ACGT", "TGCA"))),
        ("AACCgg", "ccGGTT"),
        ("", ""),
    ],
)
def test_reverse_complement_str(seq, expected_revcomp):
    assert reverse_complement(seq) == expected_revcomp


@pytest.mark.parametrize(
    "seq_bytes, expected_revcomp_bytes",
    [
        (b"ACGT", b"ACGT"[::-1].translate(bytes.maketrans(b"ACGT", b"TGCA"))),
        (b"AaNn", b"nNtT"),
    ],
)
def test_reverse_complement_bytes(seq_bytes, expected_revcomp_bytes):
    assert reverse_complement(seq_bytes) == expected_revcomp_bytes


def test_process_sequence_alias():
    # process_sequence should equal reverse_complement
    seq = "AGTC"
    assert process_sequence(seq) == reverse_complement(seq)
    seq_b = b"AGTC"
    assert process_sequence(seq_b) == reverse_complement(seq_b)
