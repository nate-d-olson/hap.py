import os
import subprocess
import sys

import pytest


@pytest.mark.parametrize(
    "mod_name, expected",
    [
        ("happy.hap", b"Haplotype Comparison"),
        ("happy.qfy", b"Usage"),
        ("happy.pre", b"Preprocessing for a VCF file"),
    ],
)
def test_module_help(mod_name, expected):
    # Run python -m modulename --help with src/python on PYTHONPATH
    cmd = [sys.executable, "-m", mod_name, "--help"]
    # Prepare environment
    env = os.environ.copy()
    # Add src/python directory for module resolution
    root_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
    src_dir = os.path.join(root_dir, "src", "python")
    if os.path.isdir(src_dir):
        prev = env.get("PYTHONPATH", "")
        env["PYTHONPATH"] = src_dir + (os.pathsep + prev if prev else "")
    try:
        output = subprocess.check_output(cmd, stderr=subprocess.STDOUT, env=env)
    except subprocess.CalledProcessError as e:
        output = e.output
    assert expected in output
