import subprocess
import sys

import pytest


@pytest.mark.parametrize(
    "mod_name, expected",
    [
        ("happy.hap", b"Haplotype Comparison"),
        ("happy.qfy", b"Usage"),
        ("happy.pre", b"Preprocessing for a VCF file"),
        ("happy.ftx", b"Usage"),
        ("happy.cnx", b"Usage"),
        ("happy.ovc", b"Usage"),
    ],
)
def test_module_help(mod_name, expected):
    # Run python -m modulename --help
    cmd = [sys.executable, "-m", mod_name, "--help"]
    try:
        output = subprocess.check_output(cmd, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        output = e.output
    assert expected in output
