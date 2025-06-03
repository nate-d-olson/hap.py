"""Heavy whole-genome smoke test (skipped by default).

This test executes the full *hap.py* pipeline on the Q100 example dataset to
validate that the modernised codebase still runs correctly on large inputs.
Execution is opt-in because the dataset is several GB and the runtime is
significant (> 10 min).
"""

from __future__ import annotations

import os
import subprocess
from pathlib import Path

import pytest

Q100DIR = Path("example/whole-genome-test").resolve()


@pytest.mark.heavy
def test_q100_smoke(tmp_path: Path) -> None:  # pragma: no cover – heavy path
    if os.environ.get("RUN_HEAVY") != "1":
        pytest.skip("RUN_HEAVY env not set – skipping heavy dataset test")

    cmd = [
        "hap.py",
        "--threads",
        "2",  # limit during CI to keep job short
        "--engine",
        "vcfeval",
        "--gender",
        "male",
        "-r",
        str(Q100DIR / "GRCh38.fa"),
        "-f",
        str(Q100DIR / "GRCh38_HG2-T2TQ100-V1.1_smvar_dipcall-z2k.benchmark.bed"),
        "--stratification",
        str(Q100DIR / "GRCh38@all" / "GRCh38-all-stratifications.tsv"),
        "-o",
        str(tmp_path / "q100-test"),
        "--quiet",
        str(Q100DIR / "GRCh38_HG2-T2TQ100-V1.1_smvar_dipcall-z2k.vcf.gz"),
        str(Q100DIR / "GRCh38_HG2-DRAGENv4.3.4-smvar.vcf.gz"),
    ]

    # Ensure hap.py is importable via the editable install used in tests.
    env = os.environ.copy()
    env["PYTHONUNBUFFERED"] = "1"

    result = subprocess.run(cmd, env=env, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr[-1000:]

    # Quick sanity: annotated VCF and summary CSV exist
    assert (tmp_path / "q100-test.vcf.gz").exists()
    assert (tmp_path / "q100-test.vcf.gz.tbi").exists()
    assert (tmp_path / "q100-test.summary.csv").exists()
