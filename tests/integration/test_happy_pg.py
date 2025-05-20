"""
Integration tests for happy_pg_test functionality.
Migrated from src/sh/run_happy_pg_test.sh
"""

import pytest
from tests.utils import (
    check_vcfeval_availability,
    compare_summary_files,
    get_bin_dir,
    get_example_dir,
    get_python_executable,
    run_shell_command,
)


@pytest.mark.integration
@pytest.mark.slow
def test_happy_pg_test(tmp_path):
    """Test PG evaluation for hap.py with different engine options."""
    # Get paths to required files and tools
    bin_dir = get_bin_dir()
    example_dir = get_example_dir()
    python_exe = get_python_executable()

    # Prepare file paths
    reference_fa = example_dir / "chr21.fa"
    truth_vcf = example_dir / "happy" / "PG_NA12878_chr21.vcf.gz"
    query_vcf = example_dir / "happy" / "NA12878_chr21.vcf.gz"
    conf_bed = example_dir / "happy" / "PG_Conf_chr21.bed.gz"
    expected_summary = example_dir / "happy" / "expected.summary.csv"
    expected_pass_summary = example_dir / "happy" / "expected.pass.summary.csv"
    expected_unhappy_summary = example_dir / "happy" / "expected.unhappy.summary.csv"
    expected_vcfeval_summary = example_dir / "happy" / "expected.vcfeval.summary.csv"

    # Set up temporary output prefixes
    output_prefix = tmp_path / "happy_output"
    output_vcfeval_prefix = tmp_path / "happy_output.vcfeval"
    output_pass_prefix = tmp_path / "happy_output.pass"
    output_unhappy_prefix = tmp_path / "happy_output.unhappy"

    # Check if vcfeval is available
    has_vcfeval = check_vcfeval_availability()

    # If vcfeval is available, run hap.py with vcfeval engine
    if has_vcfeval:
        vcfeval_cmd = [
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
            str(output_vcfeval_prefix),
            "--engine=vcfeval",
            "-X",
            "--force-interactive",
        ]

        cmd_str = " ".join(vcfeval_cmd)
        returncode, _, stderr = run_shell_command(cmd_str)
        assert returncode == 0, f"hap.py with vcfeval failed with error: {stderr}"

        # Compare vcfeval summary
        assert compare_summary_files(
            output_vcfeval_prefix.with_suffix(".summary.csv"),
            expected_vcfeval_summary,
        ), "vcfeval files differ"

    # Run standard hap.py
    standard_cmd = [
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
        "-X",
        "--force-interactive",
    ]

    cmd_str = " ".join(standard_cmd)
    returncode, _, stderr = run_shell_command(cmd_str)
    assert returncode == 0, f"standard hap.py failed with error: {stderr}"

    # Compare standard summary
    assert compare_summary_files(
        output_prefix.with_suffix(".summary.csv"),
        expected_summary,
    ), "Standard summary files differ"

    # Run hap.py with pass-only option
    pass_cmd = [
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
        str(output_pass_prefix),
        "-X",
        "--pass-only",
        "--force-interactive",
    ]

    cmd_str = " ".join(pass_cmd)
    returncode, _, stderr = run_shell_command(cmd_str)
    assert returncode == 0, f"hap.py with pass-only failed with error: {stderr}"

    # Compare pass-only summary
    assert compare_summary_files(
        output_pass_prefix.with_suffix(".summary.csv"),
        expected_pass_summary,
    ), "PASS summary files differ"

    # Run hap.py with unhappy and ROC option
    unhappy_cmd = [
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
        str(output_unhappy_prefix),
        "-X",
        "--unhappy",
        "--roc",
        "INFO.VQSLOD",
        "--force-interactive",
    ]

    cmd_str = " ".join(unhappy_cmd)
    returncode, _, stderr = run_shell_command(cmd_str)
    assert returncode == 0, f"hap.py with unhappy failed with error: {stderr}"

    # Compare unhappy summary
    assert compare_summary_files(
        output_unhappy_prefix.with_suffix(".summary.csv"),
        expected_unhappy_summary,
    ), "UNHAPPY summary files differ"
