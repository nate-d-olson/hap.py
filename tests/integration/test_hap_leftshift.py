import filecmp
import gzip
import os
import subprocess
import sys

import pytest


@pytest.mark.integration
def test_hap_leftshift(tmp_path):
    root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
    data_dir = os.path.join(root, "src", "data", "leftshifting_example")
    # Input files
    truth = os.path.join(data_dir, "truth.vcf")
    query = os.path.join(data_dir, "query.vcf")
    ref = os.path.join(data_dir, "ref.fa")
    # Output prefix
    prefix = tmp_path / "ls_out"
    # Run hap.py with left-shift preprocessing
    cmd = [
        sys.executable,
        "-m",
        "happy.hap",
        truth,
        query,
        "-o",
        str(prefix),
        "-X",
        "--reference",
        ref,
        "-l",
        "chrT",
        "--preprocess-truth",
        "--leftshift",
        "--force-interactive",
    ]
    subprocess.check_call(cmd)
    # Compare extended counts CSV
    out_ext = f"{prefix}.extended.csv"
    expected_ext = os.path.join(data_dir, "expected.extended.csv")
    # Invoke comparison script
    cmp_script = os.path.join(root, "src", "sh", "compare_extended.py")
    subprocess.check_call([sys.executable, cmp_script, out_ext, expected_ext])
    # Extract VCF bodies
    out_vcf = tmp_path / "ls.vcf"
    with gzip.open(f"{prefix}.vcf.gz", "rt") as fi, open(out_vcf, "w") as fo:
        for line in fi:
            if not line.startswith("#"):
                fo.write(line)
    expected_vcf = os.path.join(data_dir, "expected.vcf")
    assert filecmp.cmp(str(out_vcf), expected_vcf), "Leftshift VCF mismatch"
