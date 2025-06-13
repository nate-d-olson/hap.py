import gzip
from pathlib import Path

import pytest

from hap_py.tools.bcftools import countVCFRows, parseStats, runShellCommand


def test_run_shell_command_success():
    stdout, stderr, rc = runShellCommand("echo", "hello")
    assert stdout.strip() == "hello"
    assert rc == 0


def test_run_shell_command_failure():
    with pytest.raises(Exception):
        runShellCommand("nonexistent_command_xyz")


def test_parse_stats():
    output = "SN\t0\tnumber of records:\t10\nSN\t0\tnumber of indels:\t3\n"
    df = parseStats(output)
    assert df.loc[df["type"] == "records", "count"].iloc[0] == 10
    assert df.loc[df["type"] == "indels", "count"].iloc[0] == 3


def test_count_vcf_rows(tmp_path: Path):
    text = (
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr1\t1\t.\tA\tC\t.\tPASS\t.\n"
        "chr1\t2\t.\tG\tT\t.\tPASS\t.\n"
    )
    vcf = tmp_path / "a.vcf"
    vcf.write_text(text)
    assert countVCFRows(str(vcf)) == 2

    gz_path = tmp_path / "a.vcf.gz"
    with gzip.open(gz_path, "wt") as f:
        f.write(text)
    assert countVCFRows(str(gz_path)) == 2
