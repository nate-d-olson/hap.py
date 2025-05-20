#!/usr/bin/env python3
"""
Integration test for the ftx.py CLI script.

This test verifies that the ftx.py script can process VCF files correctly
and produce the expected output files.
"""

import os
import subprocess
import sys
import tempfile
import pytest
from pathlib import Path


@pytest.mark.integration
def test_ftx_basic_cli():
    """Test basic ftx.py CLI functionality."""
    print("Starting ftx.py integration test...")

    # Set up temporary output directory
    with tempfile.TemporaryDirectory() as temp_dir:
        # Use example VCF file if available
        example_dir = Path(__file__).resolve().parent.parent / "example"
        if (example_dir / "hc.vcf.gz").exists():
            input_vcf = str(example_dir / "hc.vcf.gz")
            print(f"Using example VCF: {input_vcf}")
        else:
            print("Warning: Example VCF file not found, using a dummy input")
            # Create a dummy input file
            input_vcf = os.path.join(temp_dir, "dummy.vcf")
            with open(input_vcf, "w") as f:
                f.write("##fileformat=VCFv4.2\n")
                f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
                f.write("1\t1000\t.\tA\tG\t100\tPASS\t.\n")

        # Create output filename
        output_csv = os.path.join(temp_dir, "features.csv")
        print(f"Output will be written to: {output_csv}")

        # Get the path to the ftx.py script
        script_dir = Path(__file__).resolve().parent.parent / "src" / "python"
        ftx_script = script_dir / "ftx.py"
        print(f"Using ftx.py script: {ftx_script}")

        # Run ftx.py with mock environment
        env = os.environ.copy()
        env["HAPLO_USE_MOCK"] = "1"

        # Run the command with minimal required arguments
        cmd = [
            sys.executable,
            str(ftx_script),
            input_vcf,
            "-o",
            output_csv,
            "--feature-table",
            "generic",
        ]

        print(f"Running command: {' '.join(cmd)}")
        try:
            result = subprocess.run(
                cmd,
                env=env,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                check=False,
            )

            print(f"Exit code: {result.returncode}")
            print(f"STDOUT: {result.stdout}")
            print(f"STDERR: {result.stderr}")

            # Check if command executed successfully
            if result.returncode != 0:
                print(
                    f"ERROR: ftx.py command failed with exit code {result.returncode}"
                )
                return False

            # Check if output file was created
            if not os.path.exists(output_csv):
                print(f"ERROR: Expected output file {output_csv} not found")
                return False

            print(f"SUCCESS: Output file created: {output_csv}")
            return True
        except Exception as e:
            print(f"ERROR: Failed to run ftx.py: {str(e)}")
            return False


if __name__ == "__main__":
    print("Testing ftx.py script...")
    success = test_ftx_basic_cli()

    if success:
        print("\nAll ftx.py tests passed!")
        sys.exit(0)
    else:
        print("\nSome ftx.py tests failed!")
        sys.exit(1)
