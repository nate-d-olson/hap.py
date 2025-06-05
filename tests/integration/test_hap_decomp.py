import filecmp
import gzip
import os
import subprocess
import sys

import pytest


@pytest.mark.integration
def test_hap_decomposition(tmp_path):
    root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
    example = os.path.join(root, "example", "decomp")
    # Input files
    truth = os.path.join(example, "decomp_test.truth.vcf.gz")
    query = os.path.join(example, "decomp_test.query.vcf.gz")
    conf_bed = os.path.join(example, "decomp_test.conf.bed.gz")
    # Output prefix
    prefix = tmp_path / "decomp_out"
    # Run hap.py with decomposition preprocessing
    cmd = [
        sys.executable,
        "-m",
        "happy.hap",
        "-T",
        "20",
        truth,
        query,
        "-f",
        conf_bed,
        "-o",
        str(prefix),
        "--preprocess-truth",
        "-X",
        "-V",
        "--force-interactive",
    ]
    subprocess.check_call(cmd)
    # Compare summary
    out_summary = f"{prefix}.summary.csv"
    expected_summary = os.path.join(example, "expected.summary.csv")
    assert filecmp.cmp(out_summary, expected_summary), "Summary CSV mismatch"
    # Compare VCF (ignore header lines)
    out_vcf = tmp_path / "decomp.vcf"
    # Extract body lines from output VCF
    with gzip.open(f"{prefix}.vcf.gz", "rt") as fi:
        out_lines = [l for l in fi if not l.startswith("#")]
    # Extract body lines from expected VCF
    expected_vcf = os.path.join(example, "expected.vcf")
    with open(expected_vcf) as fe:
        exp_lines = [l for l in fe if not l.startswith("#")]
    assert out_lines == exp_lines, "Decomposition VCF mismatch"
