"""
Integration tests for the qfy.py module.

These tests verify that qfy.py can correctly quantify
variant comparison results.
"""

import os
import subprocess
import sys
from pathlib import Path


def test_qfy_basic(tmp_path):
    """Test basic qfy.py functionality with a GA4GH VCF file."""
    # Create a mock GA4GH VCF with expected fields
    ga4gh_vcf = tmp_path / "ga4gh.vcf"
    with open(ga4gh_vcf, "w", encoding="utf-8") as f:
        f.write(
            """##fileformat=VCFv4.2
##INFO=<ID=BS,Number=1,Type=Integer,Description="Benchmarking superlocus ID">
##INFO=<ID=Regions,Number=.,Type=String,Description="Regions">
##INFO=<ID=Subtype,Number=1,Type=String,Description="Variant subtype">
##INFO=<ID=Type,Number=1,Type=String,Description="Variant type">
##INFO=<ID=TP,Number=0,Type=Flag,Description="True positive">
##INFO=<ID=FP,Number=0,Type=Flag,Description="False positive">
##INFO=<ID=FN,Number=0,Type=Flag,Description="False negative">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=BD,Number=1,Type=String,Description="Decision">
##FORMAT=<ID=BK,Number=1,Type=String,Description="Decision subtype">
##FORMAT=<ID=BI,Number=1,Type=String,Description="Additional info">
##FORMAT=<ID=QQ,Number=1,Type=Float,Description="Quality">
##FORMAT=<ID=BVT,Number=1,Type=String,Description="Variant type">
##FORMAT=<ID=BLT,Number=1,Type=String,Description="Location type">
##contig=<ID=chr1,length=248956422>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTRUTH\tQUERY
chr1\t100\t.\tA\tT\t50\tPASS\tBS=1;Type=SNP;Subtype=SNP;TP\tGT:BD:BK:BI:QQ:BVT:BLT\t0/1:TP:gm:tv:100:SNP:het\t0/1:TP:gm:tv:100:SNP:het
chr1\t200\t.\tG\tC\t50\tPASS\tBS=2;Type=SNP;Subtype=SNP;TP\tGT:BD:BK:BI:QQ:BVT:BLT\t1/1:TP:gm:tv:100:SNP:homalt\t0/1:TP:gm:tv:90:SNP:het
chr1\t300\t.\tC\tG\t50\tPASS\tBS=3;Type=SNP;Subtype=SNP;FN\tGT:BD:BK:BI:QQ:BVT:BLT\t0/1:FN:gm:tv:80:SNP:het\t./.:.:.:.:.:.:.
chr1\t400\t.\tT\tA\t50\tPASS\tBS=4;Type=SNP;Subtype=SNP;FP\tGT:BD:BK:BI:QQ:BVT:BLT\t./.:.:.:.:.:.:.\t0/1:FP:gm:tv:70:SNP:het
"""
        )

    # Create output directory
    output_prefix = str(tmp_path / "qfy_output")

    # Get the path to the qfy.py script
    script_dir = Path(__file__).resolve().parent.parent
    qfy_script = script_dir / "src" / "hap_py" / "qfy.py"

    # Run qfy.py with mock environment
    env = os.environ.copy()

    # Run the command
    cmd = [
        sys.executable,
        str(qfy_script),
        "--force-interactive",  # Avoid SGE requirements
        "-i",
        str(ga4gh_vcf),
        "-o",
        output_prefix,
        "-t",
        "ga4gh",
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


def test_qfy_roc(tmp_path):
    """Test qfy.py ROC functionality."""
    # Create a mock GA4GH VCF with QQ scores for ROC
    ga4gh_vcf = tmp_path / "ga4gh_roc.vcf"
    with open(ga4gh_vcf, "w", encoding="utf-8") as f:
        f.write(
            """##fileformat=VCFv4.2
##INFO=<ID=BS,Number=1,Type=Integer,Description="Benchmarking superlocus ID">
##INFO=<ID=Regions,Number=.,Type=String,Description="Regions">
##INFO=<ID=Subtype,Number=1,Type=String,Description="Variant subtype">
##INFO=<ID=Type,Number=1,Type=String,Description="Variant type">
##INFO=<ID=TP,Number=0,Type=Flag,Description="True positive">
##INFO=<ID=FP,Number=0,Type=Flag,Description="False positive">
##INFO=<ID=FN,Number=0,Type=Flag,Description="False negative">
##INFO=<ID=QQ,Number=1,Type=Float,Description="Quality score for ROC">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=BD,Number=1,Type=String,Description="Decision">
##FORMAT=<ID=BK,Number=1,Type=String,Description="Decision subtype">
##FORMAT=<ID=BI,Number=1,Type=String,Description="Additional info">
##FORMAT=<ID=QQ,Number=1,Type=Float,Description="Quality">
##FORMAT=<ID=BVT,Number=1,Type=String,Description="Variant type">
##FORMAT=<ID=BLT,Number=1,Type=String,Description="Location type">
##contig=<ID=chr1,length=248956422>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTRUTH\tQUERY
chr1\t100\t.\tA\tT\t50\tPASS\tBS=1;Type=SNP;Subtype=SNP;TP;QQ=100\tGT:BD:BK:BI:QQ:BVT:BLT\t0/1:TP:gm:tv:100:SNP:het\t0/1:TP:gm:tv:100:SNP:het
chr1\t200\t.\tG\tC\t50\tPASS\tBS=2;Type=SNP;Subtype=SNP;TP;QQ=90\tGT:BD:BK:BI:QQ:BVT:BLT\t1/1:TP:gm:tv:100:SNP:homalt\t0/1:TP:gm:tv:90:SNP:het
chr1\t300\t.\tC\tG\t50\tPASS\tBS=3;Type=SNP;Subtype=SNP;FN;QQ=80\tGT:BD:BK:BI:QQ:BVT:BLT\t0/1:FN:gm:tv:80:SNP:het\t./.:.:.:.:.:.:.
chr1\t400\t.\tT\tA\t50\tPASS\tBS=4;Type=SNP;Subtype=SNP;FP;QQ=70\tGT:BD:BK:BI:QQ:BVT:BLT\t./.:.:.:.:.:.:.\t0/1:FP:gm:tv:70:SNP:het
chr1\t500\t.\tG\tT\t50\tPASS\tBS=5;Type=SNP;Subtype=SNP;FP;QQ=50\tGT:BD:BK:BI:QQ:BVT:BLT\t./.:.:.:.:.:.:.\t0/1:FP:gm:tv:50:SNP:het
"""
        )

    # Create output directory
    output_prefix = str(tmp_path / "qfy_roc_output")

    # Get the path to the qfy.py script
    script_dir = Path(__file__).resolve().parent.parent
    qfy_script = script_dir / "src" / "hap_py" / "qfy.py"

    # Run qfy.py with mock environment
    env = os.environ.copy()

    # Run the command with ROC output
    cmd = [
        sys.executable,
        str(qfy_script),
        "--force-interactive",  # Avoid SGE requirements
        "-i",
        str(ga4gh_vcf),
        "-o",
        output_prefix,
        "-t",
        "ga4gh",
        "--roc",
        "QQ",  # Use QQ field for ROC curve
    ]

    result = subprocess.run(cmd, env=env, capture_output=True, text=True, check=False)

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
