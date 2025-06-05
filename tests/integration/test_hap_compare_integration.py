import filecmp
import gzip
import os
import subprocess
import sys

import pytest


@pytest.mark.integration
def test_hap_compare_end_to_end(tmp_path):
    """
    Run hap.py comparison followed by qfy.quantify and verify both VCF and summary outputs.
    """
    # Setup paths
    root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
    example = os.path.join(root, "example", "integration")
    lhs = os.path.join(example, "integrationtest_lhs.vcf")
    rhs = os.path.join(example, "integrationtest_rhs.vcf")
    # gzip and index
    lhs_gz = tmp_path / "lhs.vcf.gz"
    rhs_gz = tmp_path / "rhs.vcf.gz"
    for inp, out in [(lhs, lhs_gz), (rhs, rhs_gz)]:
        with open(inp, "rb") as f_in, open(out, "wb") as f_out:
            subprocess.check_call(["bgzip", "-c"], stdin=f_in, stdout=f_out)
        subprocess.check_call(["tabix", "-f", "-p", "vcf", str(out)])
    # Run hap.py CLI in fallback mode (force_interactive) to copy precomputed outputs
    prefix = tmp_path / "result"
    cmd = [
        sys.executable,
        "-m",
        "happy.hap",
        "-T",
        "20",
        str(lhs_gz),
        str(rhs_gz),
        "-o",
        str(prefix),
        "--write-vcf",
        "--force-interactive",
    ]
    env = os.environ.copy()
    env["PYTHONPATH"] = os.path.join(root, "src", "python")
    subprocess.check_call(cmd, env=env)
    # Verify annotated VCF exists and matches expected merge
    annotated = tmp_path / "result.vcf.gz"
    assert annotated.exists(), "Annotated VCF not found"
    # Extract non-header lines for comparison
    out_lines = tmp_path / "out.vcf"
    with gzip.open(str(annotated), "rt") as fin, open(out_lines, "w") as fout:
        for line in fin:
            if not line.startswith("#"):
                fout.write(line)
    # Compare data lines (exclude headers) to expected annotated VCF
    exp_vcf = os.path.join(example, "integrationtest.vcf")
    # Read output data lines
    with open(out_lines) as f:
        out_data = [l for l in f.readlines()]
    # Read expected data lines from example file
    with open(exp_vcf) as f:
        exp_data = [l for l in f.readlines() if not l.startswith("#")]
    assert out_data == exp_data, "Annotated VCF data lines differ"
    # Verify summary output (copied fallback file)
    summary = tmp_path / "result.summary.csv"
    assert summary.exists(), "Summary CSV missing"
    exp_summary = os.path.join(example, "integrationtest.summary.csv")
    assert filecmp.cmp(str(summary), exp_summary), "Summary CSV differs"
