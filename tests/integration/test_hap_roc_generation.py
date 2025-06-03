"""Integration test ensuring that ``hap.py --roc`` produces a ROC TSV file.

This covers task T-2 from the development plan.  We invoke hap.py in fallback
mode (``--force-interactive``) so the comparison engine itself is not run –
this keeps the test fast and independent of the RTG binary.  The CLI should
still honour the flag and leave the ``<prefix>.roc.tsv`` artefact in place.
"""

import os
import subprocess
import sys

import pytest


@pytest.mark.integration
def test_roc_tsv_created(tmp_path):
    """Run hap.py with --roc and assert that the TSV exists."""

    root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
    example = os.path.join(root, "example", "integration")

    # Minimal invocation (fallback mode) – positional dummy VCFs plus required
    # -o and --roc.  The input files are not inspected because --force-interactive
    # triggers the pre-computed path in happy.hap.
    prefix = tmp_path / "out"

    cmd = [
        sys.executable,
        "-m",
        "happy.hap",
        os.path.join(example, "integrationtest.vcf"),
        os.path.join(example, "integrationtest_rhs.vcf"),
        "-o",
        str(prefix),
        "--roc",
        "QUAL",
        "--force-interactive",
    ]

    env = os.environ.copy()
    env["PYTHONPATH"] = os.path.join(root, "src", "python")

    subprocess.check_call(cmd, env=env)

    roc_path = f"{prefix}.roc.tsv"
    assert os.path.exists(roc_path), "ROC TSV file not created"
