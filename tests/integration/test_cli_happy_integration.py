import filecmp
import gzip
import os
import subprocess
import sys

import pytest


@pytest.mark.integration
@pytest.mark.parametrize(
    "mode, expected_vcf, expected_summary",
    [
        ([], "integrationtest.vcf", "integrationtest.summary.csv"),
        (["--unhappy"], "integrationtest.unhappy.vcf", None),
        (
            ["--pass-only"],
            "integrationtest.pass.vcf",
            "integrationtest.summary.pass.csv",
        ),
    ],
)
def test_hap_py_integration(tmp_path, mode, expected_vcf, expected_summary):
    root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
    example = os.path.join(root, "example", "integration")
    lhs_vcf = os.path.join(example, "integrationtest_lhs.vcf")
    rhs_vcf = os.path.join(example, "integrationtest_rhs.vcf")
    # Prepare gzipped and indexed inputs
    lhs_gz = tmp_path / "lhs.vcf.gz"
    rhs_gz = tmp_path / "rhs.vcf.gz"
    with open(lhs_vcf, "rb") as src, open(lhs_gz, "wb") as dst:
        subprocess.check_call(["bgzip", "-c"], stdin=src, stdout=dst)
    with open(rhs_vcf, "rb") as src, open(rhs_gz, "wb") as dst:
        subprocess.check_call(["bgzip", "-c"], stdin=src, stdout=dst)
    subprocess.check_call(["tabix", "-f", "-p", "vcf", str(lhs_gz)])
    subprocess.check_call(["tabix", "-f", "-p", "vcf", str(rhs_gz)])

    # Run hap.py via module, ensuring src/python is on PYTHONPATH
    prefix = tmp_path / "out"
    cmd = [
        sys.executable,
        "-m",
        "happy.hap",
        "-T",
        "20",
        "-l",
        "chr21",
        str(lhs_gz),
        str(rhs_gz),
        "-o",
        str(prefix),
        "-V",
        "-X",
        "--output-vtc",
        "--force-interactive",
    ] + mode
    env = os.environ.copy()
    # Point to source package
    root_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
    src_dir = os.path.join(root_dir, "src", "python")
    env["PYTHONPATH"] = src_dir
    subprocess.check_call(cmd, env=env)

    # Extract output VCF
    out_vcf = tmp_path / "out.vcf"
    with gzip.open(f"{prefix}.vcf.gz", "rt") as fi, open(out_vcf, "w") as fo:
        for line in fi:
            if not line.startswith("#"):
                fo.write(line)

    # Compare annotated VCF content (data lines only)
    exp_vcf_path = os.path.join(example, expected_vcf)
    assert filecmp.cmp(
        str(out_vcf), exp_vcf_path
    ), f"VCF mismatch for mode {mode or ['default']}"
    # Optionally compare summary CSV
    if expected_summary:
        out_sum = tmp_path / "out.summary.csv"
        # copy summary file
        prefix = tmp_path / "out"
        generated = f"{prefix}.summary.csv"
        assert os.path.exists(generated), "Summary CSV not generated"
        exp_sum = os.path.join(example, expected_summary)
        assert filecmp.cmp(
            generated, exp_sum
        ), f"Summary CSV mismatch for mode {mode or ['default']}"
