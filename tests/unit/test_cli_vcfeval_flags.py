"""Unit tests for newly exposed CLI flags (tasks T-1).

The tests run the ``happy.hap`` module with ``--help`` or error-eliciting
arguments to ensure that:

* The new options ``--threads``, ``--roc / --no-roc`` and
  ``--Xloose-match-distance`` are surfaced in the help text.
* Invalid values cause a non-zero exit status and an informative message.
"""

import subprocess
import sys


def _run_hap_py(*extra_args):  # type: ignore[typing-arg-types]
    """Helper: run ``python -m happy.hap`` with the provided arguments."""

    cmd = [sys.executable, "-m", "happy.hap", "--help"] + list(extra_args)
    return subprocess.run(cmd, capture_output=True, text=True, check=False)


def test_help_includes_new_flags():
    """The ``-h/--help`` output must mention the new flags."""

    res = _run_hap_py()
    assert res.returncode == 0, res.stderr

    help_text = res.stdout
    # Presence of the three new CLI options in help text.
    for token in ("--threads", "--roc", "--no-roc", "--Xloose-match-distance"):
        assert token in help_text, f"Help output missing {token}"


def test_invalid_threads_value_errors():
    """Non-integer values for ``--threads`` should exit non-zero."""

    cmd = [
        sys.executable,
        "-m",
        "happy.hap",
        "--threads",
        "not_an_int",
    ]
    res = subprocess.run(cmd, capture_output=True, text=True, check=False)

    assert res.returncode != 0, "Expect non-zero exit code for bad threads value"
    assert "invalid int value" in res.stderr.lower() or "error" in res.stderr.lower()
