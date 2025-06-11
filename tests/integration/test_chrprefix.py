"""
Integration tests for chromosome prefix detection.
Migrated from src/sh/run_chrprefix_test.sh
"""

import filecmp
import gzip
import subprocess
from pathlib import Path

import pytest

from tests.utils import (
    compare_files_content,
    get_bin_dir,
    get_python_executable,
    get_src_data_dir,
    validate_output_files,
)


@pytest.mark.integration
def test_numeric_chrs(tmp_path, rtg_executable):
    """Test chr prefix detection with numeric chromosomes."""
    # Get paths to required files
    src_data_dir = get_src_data_dir() / "numeric_chrs"
    # Using CLI commands instead of script paths

    # Input files
    truth_vcf = src_data_dir / "truth.vcf"
    query_vcf = src_data_dir / "query.vcf"
    fp_bed = src_data_dir / "fp.bed"
    reference = src_data_dir / "test.fa"
    expected_vcf = src_data_dir / "expected.vcf"
    expected_summary = src_data_dir / "expected.summary.csv"

    # Check that required files exist
    assert truth_vcf.exists(), f"Truth VCF {truth_vcf} not found"
    assert query_vcf.exists(), f"Query VCF {query_vcf} not found"
    assert fp_bed.exists(), f"FP BED {fp_bed} not found"
    assert reference.exists(), f"Reference {reference} not found"
    assert expected_vcf.exists(), f"Expected VCF {expected_vcf} not found"
    assert expected_summary.exists(), f"Expected summary {expected_summary} not found"

    # Output file paths
    output_prefix = tmp_path / "numeric_test"
    output_vcf_gz = output_prefix.with_suffix(".vcf.gz")
    output_vcf = output_prefix.with_suffix(".vcf")
    output_summary = Path(str(output_prefix) + ".summary.csv")

    # Run hap.py on numeric chromosome files using CLI command
    python_exe = get_python_executable()
    hap_py_script = get_bin_dir() / "hap.py"
    cmd = [
        python_exe,
        str(hap_py_script),
        str(truth_vcf),
        str(query_vcf),
        "-f",
        str(fp_bed),
        "-o",
        str(output_prefix),
        "-X",
        "--reference",
        str(reference),
        "-V",
        "--force-interactive",
        "--engine-vcfeval-path",
        rtg_executable,  # Use fixture instead of hardcoded path
    ]

    result = subprocess.run(cmd, capture_output=True, check=True)
    assert (
        result.returncode == 0
    ), f"hap.py failed with numeric chromosomes: {result.stderr.decode()}"

    # Validate output files with robust checking
    required_files = [".summary.csv", ".vcf.gz"]
    optional_files = [".roc.tsv", ".extended.csv", ".metrics.json.gz"]

    missing_required, failed_comparisons = validate_output_files(
        str(output_prefix), required_files, optional_files
    )

    assert not missing_required, f"Missing required output files: {missing_required}"

    # Compare summary files
    assert compare_files_content(
        str(output_summary), str(expected_summary)
    ), "Summary output differs from expected"

    # Compare VCF files - extract content from gzipped output
    with gzip.open(output_vcf_gz, "rt") as f_gz:
        vcf_content = [line for line in f_gz if not line.startswith("#")]

    with open(output_vcf, "w", encoding="utf-8") as f_out:
        f_out.writelines(vcf_content)

    assert compare_files_content(
        str(output_vcf), str(expected_vcf)
    ), "VCF output differs from expected"


@pytest.mark.integration
def test_chr_prefixed(tmp_path, rtg_executable):
    """Test chr prefix detection with 'chr' prefixed chromosomes."""
    # Get paths to required files
    src_data_dir = get_src_data_dir() / "numeric_chrs"
    # Using CLI commands instead of script paths

    # Input files
    truth_vcf = src_data_dir / "chrtruth.vcf"
    query_vcf = src_data_dir / "chrquery.vcf"
    fp_bed = src_data_dir / "chrfp.bed"
    reference = src_data_dir / "chrtest.fa"
    expected_vcf = src_data_dir / "chrexpected.vcf"
    expected_summary = src_data_dir / "expected.summary.csv"

    # Check that required files exist
    assert truth_vcf.exists(), f"Truth VCF {truth_vcf} not found"
    assert query_vcf.exists(), f"Query VCF {query_vcf} not found"
    assert fp_bed.exists(), f"FP BED {fp_bed} not found"
    assert reference.exists(), f"Reference {reference} not found"
    assert expected_vcf.exists(), f"Expected VCF {expected_vcf} not found"
    assert expected_summary.exists(), f"Expected summary {expected_summary} not found"

    # Output file paths
    output_prefix = tmp_path / "chr_test"
    output_vcf_gz = output_prefix.with_suffix(".vcf.gz")
    output_vcf = output_prefix.with_suffix(".vcf")
    output_summary = Path(str(output_prefix) + ".summary.csv")

    # Run hap.py on chr-prefixed files using CLI command
    python_exe = get_python_executable()
    hap_py_script = get_bin_dir() / "hap.py"
    cmd = [
        python_exe,
        str(hap_py_script),
        str(truth_vcf),
        str(query_vcf),
        "-f",
        str(fp_bed),
        "-o",
        str(output_prefix),
        "-X",
        "--reference",
        str(reference),
        "-V",
        "--force-interactive",
        "--engine-vcfeval-path",
        rtg_executable,  # Use fixture instead of hardcoded path
    ]

    result = subprocess.run(cmd, capture_output=True, check=True)
    assert (
        result.returncode == 0
    ), f"hap.py failed with chr-prefixed chromosomes: {result.stderr.decode()}"

    # Compare summary files using a simple diff approach
    # Instead of calling a separate compare script, we'll compare directly
    with open(output_summary, encoding="utf-8") as f_out:
        output_lines = f_out.readlines()
    with open(expected_summary, encoding="utf-8") as f_exp:
        expected_lines = f_exp.readlines()

    assert output_lines == expected_lines, "Summary output differs from expected"

    # Compare VCF files
    with gzip.open(output_vcf_gz, "rt") as f_gz:
        vcf_content = [line for line in f_gz if not line.startswith("#")]

    with open(output_vcf, "w", encoding="utf-8") as f_out:
        f_out.writelines(vcf_content)

    assert filecmp.cmp(output_vcf, expected_vcf), "VCF output differs from expected"


@pytest.mark.integration
def test_mixed_chr_prefix(tmp_path, rtg_executable):
    """Test chr prefix detection with mixed chromosome naming \\
    (chr in truth, numeric in query)."""
    # Get paths to required files
    src_data_dir = get_src_data_dir() / "numeric_chrs"
    # Using CLI commands instead of script paths

    # Input files - mixing chr in truth with numeric in query
    truth_vcf = src_data_dir / "chrtruth.vcf"
    query_vcf = src_data_dir / "query.vcf"  # numeric chromosomes
    fp_bed = src_data_dir / "chrfp.bed"
    reference = src_data_dir / "chrtest.fa"
    expected_vcf = src_data_dir / "chrexpected.vcf"
    expected_summary = src_data_dir / "expected.summary.csv"

    # Check that required files exist
    assert truth_vcf.exists(), f"Truth VCF {truth_vcf} not found"
    assert query_vcf.exists(), f"Query VCF {query_vcf} not found"
    assert fp_bed.exists(), f"FP BED {fp_bed} not found"
    assert reference.exists(), f"Reference {reference} not found"
    assert expected_vcf.exists(), f"Expected VCF {expected_vcf} not found"
    assert expected_summary.exists(), f"Expected summary {expected_summary} not found"

    # Output file paths
    output_prefix = tmp_path / "mixed_test"
    output_vcf_gz = output_prefix.with_suffix(".vcf.gz")
    output_vcf = output_prefix.with_suffix(".vcf")
    output_summary = Path(str(output_prefix) + ".summary.csv")

    # Run hap.py with mixed chromosome naming using CLI command
    python_exe = get_python_executable()
    hap_py_script = get_bin_dir() / "hap.py"
    cmd = [
        python_exe,
        str(hap_py_script),
        str(truth_vcf),
        str(query_vcf),
        "-f",
        str(fp_bed),
        "-o",
        str(output_prefix),
        "-X",
        "--reference",
        str(reference),
        "-V",
        "--force-interactive",
        "--engine-vcfeval-path",
        rtg_executable,  # Use fixture instead of hardcoded path
    ]

    result = subprocess.run(cmd, capture_output=True, check=True)
    assert (
        result.returncode == 0
    ), f"hap.py failed with mixed chromosome naming: {result.stderr.decode()}"

    # Compare summary files using a simple diff approach
    # Instead of calling a separate compare script, we'll compare directly
    with open(output_summary, encoding="utf-8") as f_out:
        output_lines = f_out.readlines()
    with open(expected_summary, encoding="utf-8") as f_exp:
        expected_lines = f_exp.readlines()

    assert output_lines == expected_lines, "Summary output differs from expected"

    # Compare VCF files
    with gzip.open(output_vcf_gz, "rt") as f_gz:
        vcf_content = [line for line in f_gz if not line.startswith("#")]

    with open(output_vcf, "w", encoding="utf-8") as f_out:
        f_out.writelines(vcf_content)

    assert filecmp.cmp(output_vcf, expected_vcf), "VCF output differs from expected"
