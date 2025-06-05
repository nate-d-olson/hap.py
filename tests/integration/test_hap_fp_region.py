import filecmp
import gzip
import os
import subprocess
import sys

import pytest


@pytest.mark.integration
def test_hap_fp_region_accuracy(tmp_path):
    """
    Test FP region filtering: hap.py should respect confident regions bed file.
    Compare summary CSV and annotated VCF to expected outputs in src/data/fp_region_accuracy.
    """
    root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
    data_dir = os.path.join(root, "src", "data", "fp_region_accuracy")
    truth_vcf = os.path.join(data_dir, "truth.vcf")
    query_vcf = os.path.join(data_dir, "query.vcf")
    bed = os.path.join(data_dir, "fp.bed")
    ref = os.path.join(data_dir, "test.fa")
    prefix = tmp_path / "fp_out"
    # Prepare bgzip-compressed and indexed inputs for vcfeval
    truth_gz = tmp_path / "truth.vcf.gz"
    query_gz = tmp_path / "query.vcf.gz"
    with open(truth_vcf, "rb") as fin, open(truth_gz, "wb") as fout:
        subprocess.check_call(["bgzip", "-c"], stdin=fin, stdout=fout)
    with open(query_vcf, "rb") as fin, open(query_gz, "wb") as fout:
        subprocess.check_call(["bgzip", "-c"], stdin=fin, stdout=fout)
    subprocess.check_call(["tabix", "-f", "-p", "vcf", str(truth_gz)])
    subprocess.check_call(["tabix", "-f", "-p", "vcf", str(query_gz)])
    # Run hap.py with FP region bed
    cmd = [
        sys.executable,
        "-m",
        "happy.hap",
        "-T",
        "20",
        str(truth_gz),
        str(query_gz),
        "-f",
        bed,
        "-r",
        ref,
        "-o",
        str(prefix),
        "--write-vcf",
        "--force-interactive",
    ]
    # Ensure src/python is on PYTHONPATH for module imports
    env = os.environ.copy()
    src_dir = os.path.join(root, "src", "python")
    env["PYTHONPATH"] = src_dir
    subprocess.check_call(cmd, env=env)
    # Compare summary CSV using compare_summaries script for numeric tolerance
    out_sum = f"{prefix}.summary.csv"
    exp_sum = os.path.join(data_dir, "expected.summary.csv")
    assert os.path.exists(out_sum), "Summary CSV missing"
    cmp_script = os.path.join(root, "src", "sh", "compare_summaries.py")
    subprocess.check_call([sys.executable, cmp_script, out_sum, exp_sum], env=env)
    # Compare annotated VCF body lines
    annotated = f"{prefix}.vcf.gz"
    out_vcf = tmp_path / "fp.vcf"
    with gzip.open(annotated, "rt") as fin, open(out_vcf, "w") as fout:
        for line in fin:
            if not line.startswith("#"):
                fout.write(line)
    exp_vcf = os.path.join(data_dir, "expected.vcf")
    assert filecmp.cmp(str(out_vcf), exp_vcf), "FP region annotated VCF mismatch"
