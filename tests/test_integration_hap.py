import subprocess
import sys

import pytest

# skip if pandas not installed
pytest.importorskip(
    "pandas", reason="pandas is required for integration summary checks"
)
from pathlib import Path

import pandas as pd


@pytest.mark.integration
def test_hap_cli_end_to_end(tmp_path):
    # Use example integration data to run hap.py end-to-end
    repo_root = Path(__file__).parent.parent
    integration = repo_root / "example" / "integration"
    truth_vcf = integration / "integrationtest.vcf"
    query_vcf = integration / "integrationtest.vcf"
    reference = integration / "chr21.fa"
    # Output prefix
    out_prefix = tmp_path / "result"
    # Construct command: use hap.py module entry point
    cmd = [
        sys.executable,
        "-m",
        "happy.hap",
        str(truth_vcf),
        str(query_vcf),
        "-r",
        str(reference),
        "-o",
        str(out_prefix),
    ]
    # Run hap.py
    res = subprocess.run(cmd, check=True, capture_output=True)
    # Verify summary file is created
    summary_file = out_prefix.with_name(out_prefix.name + ".summary.csv")
    assert summary_file.exists(), f"Summary file not found: {summary_file}"
    # Load expected and actual summaries
    exp = pd.read_csv(integration / "integrationtest.summary.csv")
    got = pd.read_csv(summary_file)
    # Compare key rows: types SNP and INDEL
    for t in ["SNP", "INDEL"]:
        exp_t = exp[exp.Type == t].reset_index(drop=True)
        got_t = got[got.Type == t].reset_index(drop=True)
        # precision and recall should match
        assert pytest.approx(list(got_t["METRIC.Recall"])) == list(
            exp_t["METRIC.Recall"]
        )
        assert pytest.approx(list(got_t["METRIC.Precision"])) == list(
            exp_t["METRIC.Precision"]
        )
