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
    with gzip.open(f"{prefix}.vcf.gz", "rt") as fi, open(out_vcf, "w") as fo:
        for line in fi:
            if not line.startswith("#"):
                fo.write(line)
    expected_vcf = os.path.join(example, "expected.vcf")
    assert filecmp.cmp(str(out_vcf), expected_vcf), "Decomposition VCF mismatch"
