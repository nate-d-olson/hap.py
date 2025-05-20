"""
Integration tests for blocksplit functionality.
Migrated from src/sh/run_blocksplit_test.sh
"""

import subprocess
import sys
import tempfile
from pathlib import Path

import pytest


@pytest.mark.integration
@pytest.mark.cpp
def test_blocksplit():
    """Test blocksplit functionality on VCF files."""
    # Get paths to required files
    project_root = Path(__file__).parent.parent.parent
    example_dir = project_root / "example" / "happy"
    vcf1_path = example_dir / "PG_NA12878_hg38-chr21.vcf.gz"
    vcf2_path = example_dir / "NA12878-GATK3-chr21.vcf.gz"

    # Get path to tools
    bin_dir = project_root / "build" / "bin"
    blocksplit_bin = bin_dir / "blocksplit"
    bcftools_bin = bin_dir / "bcftools"
    ovc_script = project_root / "src" / "python" / "ovc.py"

    # Check that required files exist
    assert vcf1_path.exists(), f"Test VCF1 {vcf1_path} not found"
    assert vcf2_path.exists(), f"Test VCF2 {vcf2_path} not found"
    assert blocksplit_bin.exists(), f"Blocksplit binary {blocksplit_bin} not found"
    assert bcftools_bin.exists(), f"Bcftools binary {bcftools_bin} not found"
    assert ovc_script.exists(), f"OVC script {ovc_script} not found"

    # Create temporary files
    with tempfile.NamedTemporaryFile(suffix=".bed") as temp_result:
        # Run blocksplit
        blocksplit_cmd = [
            str(blocksplit_bin),
            str(vcf1_path),
            str(vcf2_path),
            "-o",
            temp_result.name,
            "-l",
            "chr21",
            "-w",
            "10000",
        ]
        subprocess.run(blocksplit_cmd, check=True)

        # Check for overlaps using ovc.py
        ovc_cmd = [sys.executable, str(ovc_script), temp_result.name]
        subprocess.run(ovc_cmd, check=True)

        # Check that vcf1 variants are all captured in the blocks
        with tempfile.NamedTemporaryFile(
            suffix=".vcf"
        ) as tf_x1, tempfile.NamedTemporaryFile(suffix=".vcf") as tf_x2:
            # Extract all variants from vcf1
            subprocess.run(
                [str(bcftools_bin), "view", "-H", str(vcf1_path), "chr21"],
                stdout=tf_x1,
                check=True,
            )

            # Extract variants in the bed regions from vcf1
            subprocess.run(
                [
                    str(bcftools_bin),
                    "view",
                    "-H",
                    str(vcf1_path),
                    "chr21",
                    "-R",
                    temp_result.name,
                ],
                stdout=tf_x2,
                check=True,
            )

            # Compare to make sure all variants are included
            with open(tf_x1.name) as f1, open(tf_x2.name) as f2:
                assert f1.read() == f2.read(), (
                    "Not all VCF1 variants are covered by the blocks"
                )

            # Same check for vcf2
            subprocess.run(
                [str(bcftools_bin), "view", "-H", str(vcf2_path), "chr21"],
                stdout=tf_x1,
                check=True,
            )

            subprocess.run(
                [
                    str(bcftools_bin),
                    "view",
                    "-H",
                    str(vcf2_path),
                    "chr21",
                    "-R",
                    temp_result.name,
                ],
                stdout=tf_x2,
                check=True,
            )

            with open(tf_x1.name) as f1, open(tf_x2.name) as f2:
                assert f1.read() == f2.read(), (
                    "Not all VCF2 variants are covered by the blocks"
                )
