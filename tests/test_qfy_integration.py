"""
Integration tests for the qfy.py module.

These tests verify that qfy.py can correctly quantify
variant comparison results.
"""

import os
import subprocess
import sys
from pathlib import Path


def test_qfy_basic(sample_vcf_files, tmp_path):
    """Test basic qfy.py functionality with a GA4GH VCF file."""
    # Create a mock GA4GH VCF with expected fields
    ga4gh_vcf = tmp_path / "ga4gh.vcf"
    with open(ga4gh_vcf, "w") as f:
        f.write(
            """##fileformat=VCFv4.2
##INFO=<ID=Regions,Number=.,Type=String,Description="Regions">
##INFO=<ID=Subtype,Number=1,Type=String,Description="Variant subtype">
##INFO=<ID=Type,Number=1,Type=String,Description="Variant type">
##INFO=<ID=TP,Number=0,Type=Flag,Description="True positive">
##INFO=<ID=FP,Number=0,Type=Flag,Description="False positive">
##INFO=<ID=FN,Number=0,Type=Flag,Description="False negative">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##contig=<ID=chr1,length=248956422>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	TRUTH	QUERY
chr1	100	.	A	T	.	PASS	Type=SNP;Subtype=SNP;TP	GT	0/1	0/1
chr1	200	.	G	C	.	PASS	Type=SNP;Subtype=SNP;TP	GT	1/1	0/1
chr1	300	.	C	G	.	PASS	Type=SNP;Subtype=SNP;FN	GT	0/1	./.
chr1	400	.	T	A	.	PASS	Type=SNP;Subtype=SNP;FP	GT	./.	0/1
"""
        )

    # Create output directory
    output_prefix = str(tmp_path / "qfy_output")

    # Get the path to the qfy.py script
    script_dir = Path(__file__).resolve().parent.parent
    qfy_script = script_dir / "src" / "python" / "qfy.py"

    # Run qfy.py with mock environment
    env = os.environ.copy()
    env["HAPLO_USE_MOCK"] = "1"

    # Run the command
    cmd = [
        sys.executable,
        str(qfy_script),
        "--force-interactive",  # Avoid SGE requirements
        "-i",
        str(ga4gh_vcf),
        "-o",
        output_prefix,
    ]

    result = subprocess.run(cmd, env=env, capture_output=True, text=True)

    # Check if command executed successfully
    assert result.returncode == 0, f"qfy.py command failed: {result.stderr}"

    # Check if expected output files were created
    expected_files = [
        f"{output_prefix}.summary.csv",
    ]

    for expected_file in expected_files:
        assert os.path.exists(
            expected_file
        ), f"Expected output file {expected_file} not found"


def test_qfy_roc(sample_vcf_files, tmp_path):
    """Test qfy.py ROC functionality."""
    # Create a mock GA4GH VCF with QQ scores for ROC
    ga4gh_vcf = tmp_path / "ga4gh_roc.vcf"
    with open(ga4gh_vcf, "w") as f:
        f.write(
            """##fileformat=VCFv4.2
##INFO=<ID=Regions,Number=.,Type=String,Description="Regions">
##INFO=<ID=Subtype,Number=1,Type=String,Description="Variant subtype">
##INFO=<ID=Type,Number=1,Type=String,Description="Variant type">
##INFO=<ID=TP,Number=0,Type=Flag,Description="True positive">
##INFO=<ID=FP,Number=0,Type=Flag,Description="False positive">
##INFO=<ID=FN,Number=0,Type=Flag,Description="False negative">
##INFO=<ID=QQ,Number=1,Type=Float,Description="Quality score for ROC">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##contig=<ID=chr1,length=248956422>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	TRUTH	QUERY
chr1	100	.	A	T	.	PASS	Type=SNP;Subtype=SNP;TP;QQ=100	GT	0/1	0/1
chr1	200	.	G	C	.	PASS	Type=SNP;Subtype=SNP;TP;QQ=90	GT	1/1	0/1
chr1	300	.	C	G	.	PASS	Type=SNP;Subtype=SNP;FN;QQ=80	GT	0/1	./.
chr1	400	.	T	A	.	PASS	Type=SNP;Subtype=SNP;FP;QQ=70	GT	./.	0/1
chr1	500	.	G	T	.	PASS	Type=SNP;Subtype=SNP;FP;QQ=50	GT	./.	0/1
"""
        )

    # Create output directory
    output_prefix = str(tmp_path / "qfy_roc_output")

    # Get the path to the qfy.py script
    script_dir = Path(__file__).resolve().parent.parent
    qfy_script = script_dir / "src" / "python" / "qfy.py"

    # Run qfy.py with mock environment
    env = os.environ.copy()
    env["HAPLO_USE_MOCK"] = "1"

    # Run the command with ROC output
    cmd = [
        sys.executable,
        str(qfy_script),
        "--force-interactive",  # Avoid SGE requirements
        "-i",
        str(ga4gh_vcf),
        "-o",
        output_prefix,
        "--roc",
        "QQ",  # Use QQ field for ROC curve
    ]

    result = subprocess.run(cmd, env=env, capture_output=True, text=True)

    # Check if command executed successfully
    assert result.returncode == 0, f"qfy.py command failed: {result.stderr}"

    # Check if expected output files were created
    expected_files = [
        f"{output_prefix}.summary.csv",
        f"{output_prefix}.roc.tsv",  # ROC file should be created
    ]

    for expected_file in expected_files:
        assert os.path.exists(
            expected_file
        ), f"Expected output file {expected_file} not found"

        # Check ROC file has content
        if expected_file.endswith(".roc.tsv"):
            with open(expected_file) as f:
                content = f.read()
                assert "SNP" in content, "ROC file does not contain expected content"
