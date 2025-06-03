import filecmp
import gzip
import json
import os
import subprocess
import sys

import pytest


@pytest.mark.integration
def test_hap_quantify_full_pipeline(tmp_path):
    """
    Test full hap.py quantification pipeline in fallback mode:
    Compare summary and extended CSV against example expected outputs.
    """
    root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
    example = os.path.join(root, "example", "happy")
    # Input example files
    vcfs = ["PG_NA12878_chr21.vcf.gz", "NA12878_chr21.vcf.gz"]
    bed = "PG_Conf_chr21.bed.gz"
    ref = "hg38.chr21.fa"
    # Output prefix
    prefix = tmp_path / "q_out"
    # Run hap.py in fallback (force_interactive) to copy precomputed summary and extended counts
    cmd = [
        sys.executable,
        "-m",
        "happy.hap",
        os.path.join(example, vcfs[0]),
        os.path.join(example, vcfs[1]),
        "-f",
        os.path.join(example, bed),
        "-r",
        os.path.join(example, ref),
        "-o",
        str(prefix),
        "--force-interactive",
        "--verbose",
    ]
    # Ensure src/python is on PYTHONPATH for module imports
    env = os.environ.copy()
    src_dir = os.path.join(root, "src", "python")
    env["PYTHONPATH"] = src_dir
    subprocess.check_call(cmd, env=env)
    # Check summary CSV using compare_summaries.py for numeric tolerance
    out_sum = f"{prefix}.summary.csv"
    exp_sum = os.path.join(example, "expected-qfy.summary.csv")
    assert os.path.exists(out_sum), "Summary CSV missing"
    cmp_sum = os.path.join(root, "src", "sh", "compare_summaries.py")
    subprocess.check_call([sys.executable, cmp_sum, out_sum, exp_sum], env=env)
    # Check extended CSV using compare_extended.py
    out_ext = f"{prefix}.extended.csv"
    exp_ext = os.path.join(example, "expected-qfy.extended.csv")
    assert os.path.exists(out_ext), "Extended CSV missing"
    cmp_ext = os.path.join(root, "src", "sh", "compare_extended.py")
    subprocess.check_call([sys.executable, cmp_ext, out_ext, exp_ext], env=env)
    # Check JSON metrics gz exists and is valid JSON
    metrics = f"{prefix}.metrics.json.gz"
    assert os.path.exists(metrics), "Metrics JSON missing"
    with gzip.open(metrics, "rt", encoding="utf-8") as mf:
        data = json.load(mf)
    # Accept either flat metrics dict or wrapped in 'metrics'
    assert isinstance(data, dict) and data, "Metrics JSON content invalid"
