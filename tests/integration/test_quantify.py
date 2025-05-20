"""
Integration tests for quantify_test functionality.
Migrated from src/sh/run_quantify_test.sh
"""

import gzip
import json

import pytest
from tests.utils import (
    compare_summary_files,
    get_bin_dir,
    get_example_dir,
    get_project_root,
    get_python_executable,
    run_shell_command,
)


@pytest.mark.integration
@pytest.mark.slow
def test_quantify_test(tmp_path):
    """Test quantification and GA4GH intermediate file format compliance."""
    # Get paths to required files and tools
    get_project_root()
    bin_dir = get_bin_dir()
    example_dir = get_example_dir()
    python_exe = get_python_executable()

    # Prepare file paths
    reference_fa = example_dir / "chr21.fa"
    truth_vcf = example_dir / "happy" / "PG_NA12878_chr21.vcf.gz"
    query_vcf = example_dir / "happy" / "NA12878_chr21.vcf.gz"
    conf_bed = example_dir / "happy" / "PG_Conf_chr21.bed.gz"
    expected_summary = example_dir / "happy" / "expected-qfy.summary.csv"
    expected_extended = example_dir / "happy" / "expected-qfy.extended.csv"

    # Set up temporary output prefix
    output_prefix = tmp_path / "happy_output"

    # Run hap.py
    hap_py_cmd = [
        python_exe,
        str(bin_dir / "hap.py"),
        "-l",
        "chr21",
        str(truth_vcf),
        str(query_vcf),
        "-f",
        str(conf_bed),
        "-r",
        str(reference_fa),
        "-o",
        str(output_prefix),
        "--no-adjust-conf-regions",
        "--force-interactive",
        "--verbose",
    ]

    returncode, _, stderr = run_shell_command(" ".join(hap_py_cmd))
    assert returncode == 0, f"hap.py failed with error: {stderr}"

    # Compare summary precision/recall
    assert compare_summary_files(
        output_prefix.with_suffix(".summary.csv"), expected_summary
    ), "Summary files differ"

    # Compare extended precision/recall
    assert compare_summary_files(
        output_prefix.with_suffix(".extended.csv"), expected_extended
    ), "Extended files differ"

    # Run qfy.py for re-quantification using GA4GH spec
    qfy_output_prefix = output_prefix.with_suffix(".qfy")
    qfy_cmd = [
        python_exe,
        str(bin_dir / "qfy.py"),
        str(output_prefix.with_suffix(".vcf.gz")),
        "-o",
        str(qfy_output_prefix),
        "-r",
        str(reference_fa),
        "-f",
        str(conf_bed),
        "-t",
        "ga4gh",
        "-X",
        "--verbose",
    ]

    returncode, _, stderr = run_shell_command(" ".join(qfy_cmd))
    assert returncode == 0, f"qfy.py failed with error: {stderr}"

    # Extract metrics and filter out variable elements
    # Original hap.py run
    hap_metrics_file = output_prefix.with_suffix(".metrics.json.gz")
    hap_metrics_filtered = tmp_path / "hap_metrics_filtered.json"

    with gzip.open(hap_metrics_file, "rt") as f:
        metrics_data = json.load(f)
        # Remove variable elements that would cause comparison issues
        if "metrics" in metrics_data and "header" in metrics_data["metrics"]:
            if "timestamp" in metrics_data["metrics"]["header"]:
                del metrics_data["metrics"]["header"]["timestamp"]
            if "hap.py" in metrics_data["metrics"]["header"]:
                del metrics_data["metrics"]["header"]["hap.py"]

    with open(hap_metrics_filtered, "w", encoding="utf-8") as f:
        json.dump(metrics_data, f, sort_keys=True, indent=2)

    # Re-quantified run
    qfy_metrics_file = qfy_output_prefix.with_suffix(".metrics.json.gz")
    qfy_metrics_filtered = tmp_path / "qfy_metrics_filtered.json"

    with gzip.open(qfy_metrics_file, "rt") as f:
        metrics_data = json.load(f)
        # Remove variable elements that would cause comparison issues
        if "metrics" in metrics_data and "header" in metrics_data["metrics"]:
            if "timestamp" in metrics_data["metrics"]["header"]:
                del metrics_data["metrics"]["header"]["timestamp"]
            if "qfy.py" in metrics_data["metrics"]["header"]:
                del metrics_data["metrics"]["header"]["qfy.py"]

    with open(qfy_metrics_filtered, "w", encoding="utf-8") as f:
        json.dump(metrics_data, f, sort_keys=True, indent=2)

    # Compare re-quantified summary to expected
    assert compare_summary_files(
        qfy_output_prefix.with_suffix(".summary.csv"), expected_summary
    ), "Re-quantified summary files differ from expected"

    # Compare filtered metrics from both runs
    with open(hap_metrics_filtered, encoding="utf-8") as f1:
        hap_data = json.load(f1)
    with open(qfy_metrics_filtered, encoding="utf-8") as f2:
        qfy_data = json.load(f2)

    assert hap_data == qfy_data, "Re-quantified metrics are different from original run"


@pytest.mark.integration
@pytest.mark.slow
def test_quantify_test(tmp_path):
    """Test quantification and GA4GH intermediate file format compliance."""
    # Get paths to required files and tools
    project_root = get_project_root()
    bin_dir = get_bin_dir()
    example_dir = get_example_dir()
    python_exe = get_python_executable()

    # Get paths to the shell scripts directory
    project_root / "src" / "sh"

    # Prepare file paths
    reference_fa = example_dir / "chr21.fa"
    truth_vcf = example_dir / "happy" / "PG_NA12878_chr21.vcf.gz"
    query_vcf = example_dir / "happy" / "NA12878_chr21.vcf.gz"
    conf_bed = example_dir / "happy" / "PG_Conf_chr21.bed.gz"
    expected_summary = example_dir / "happy" / "expected-qfy.summary.csv"
    expected_extended = example_dir / "happy" / "expected-qfy.extended.csv"

    # Set up temporary output prefix
    output_prefix = tmp_path / "happy_output"

    # Run hap.py
    hap_py_cmd = [
        python_exe,
        str(bin_dir / "hap.py"),
        "-l",
        "chr21",
        str(truth_vcf),
        str(query_vcf),
        "-f",
        str(conf_bed),
        "-r",
        str(reference_fa),
        "-o",
        str(output_prefix),
        "--no-adjust-conf-regions",
        "--force-interactive",
        "--verbose",
    ]

    returncode, stdout, stderr = run_shell_command(" ".join(hap_py_cmd))
    assert returncode == 0, f"hap.py failed with error: {stderr}"

    # Compare summary precision/recall
    assert compare_summary_files(
        output_prefix.with_suffix(".summary.csv"), expected_summary
    ), "Summary files differ"

    # Compare extended precision/recall
    assert compare_summary_files(
        output_prefix.with_suffix(".extended.csv"), expected_extended
    ), "Extended files differ"

    # Run qfy.py for re-quantification using GA4GH spec
    qfy_output_prefix = output_prefix.with_suffix(".qfy")
    qfy_cmd = [
        python_exe,
        str(bin_dir / "qfy.py"),
        str(output_prefix.with_suffix(".vcf.gz")),
        "-o",
        str(qfy_output_prefix),
        "-r",
        str(reference_fa),
        "-f",
        str(conf_bed),
        "-t",
        "ga4gh",
        "-X",
        "--verbose",
    ]

    returncode, stdout, stderr = run_shell_command(" ".join(qfy_cmd))
    assert returncode == 0, f"qfy.py failed with error: {stderr}"

    # Extract metrics and filter out variable elements
    # Original hap.py run
    hap_metrics_file = output_prefix.with_suffix(".metrics.json.gz")
    hap_metrics_filtered = tmp_path / "hap_metrics_filtered.json"

    with gzip.open(hap_metrics_file, "rt") as f:
        metrics_data = json.load(f)
        # Remove variable elements that would cause comparison issues
        if "metrics" in metrics_data and "header" in metrics_data["metrics"]:
            if "timestamp" in metrics_data["metrics"]["header"]:
                del metrics_data["metrics"]["header"]["timestamp"]
            if "hap.py" in metrics_data["metrics"]["header"]:
                del metrics_data["metrics"]["header"]["hap.py"]

    with open(hap_metrics_filtered, "w") as f:
        json.dump(metrics_data, f, sort_keys=True, indent=2)

    # Re-quantified run
    qfy_metrics_file = qfy_output_prefix.with_suffix(".metrics.json.gz")
    qfy_metrics_filtered = tmp_path / "qfy_metrics_filtered.json"

    with gzip.open(qfy_metrics_file, "rt") as f:
        metrics_data = json.load(f)
        # Remove variable elements that would cause comparison issues
        if "metrics" in metrics_data and "header" in metrics_data["metrics"]:
            if "timestamp" in metrics_data["metrics"]["header"]:
                del metrics_data["metrics"]["header"]["timestamp"]
            if "qfy.py" in metrics_data["metrics"]["header"]:
                del metrics_data["metrics"]["header"]["qfy.py"]

    with open(qfy_metrics_filtered, "w") as f:
        json.dump(metrics_data, f, sort_keys=True, indent=2)

    # Compare re-quantified summary to expected
    assert compare_summary_files(
        qfy_output_prefix.with_suffix(".summary.csv"), expected_summary
    ), "Re-quantified summary files differ from expected"

    # Compare filtered metrics from both runs
    with open(hap_metrics_filtered) as f1, open(qfy_metrics_filtered) as f2:
        hap_data = json.load(f1)
        qfy_data = json.load(f2)
        assert hap_data == qfy_data, (
            "Re-quantified metrics are different from original run"
        )
