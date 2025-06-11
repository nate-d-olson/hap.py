"""
Integration test fixtures for hap.py.

These fixtures provide test resources for integration tests of the hap.py toolkit.
"""

import os
import shutil
from pathlib import Path

import pytest
from pytest import mark

# Markers for test categorization
mark.core = mark.CoreTest = mark.Core
mark.integration = mark.IntegrationTest = mark.Integration
mark.unit = mark.UnitTest = mark.Unit

# Tests to skip
SKIP_TESTS = [
    # Integration tests
    "test_blocksplit.py",
    "test_multimerge.py",
    "test_performance.py",
    "test_quantify_stratification.py",
    "test_pathtraversal.py",
    # Unit tests
    "test_cython_integration.py",
    "test_string_handling.py",
]


def pytest_collection_modifyitems(items):
    """Skip non-core tests."""
    for item in items:
        # Skip tests that are not essential
        if any(skip_test in str(item.fspath) for skip_test in SKIP_TESTS):
            item.add_marker(pytest.mark.skip(reason="Skipping non-essential test"))

        # Add appropriate markers
        if "integration" in str(item.fspath):
            item.add_marker(mark.integration)
        elif "unit" in str(item.fspath):
            item.add_marker(mark.unit)


def get_rtg_path():
    """Get the RTG tools path for integration tests."""
    project_root = Path(__file__).parent.parent

    # Check common environment variables first
    for var in ["RTG", "RTGTOOLS_PATH", "RTG_PATH"]:
        if var in os.environ:
            candidate = Path(os.environ[var])
            if candidate.is_dir():
                candidate = candidate / "rtg"
            if candidate.exists():
                return str(candidate)

    # Try external directory (modernized location)
    rtg_path = project_root / "external" / "rtg-tools-3.12.1" / "rtg"
    if rtg_path.exists():
        print(f"Found RTG tools at {rtg_path}")
        return str(rtg_path)

    # Try project root (where the rtg symlink might be)
    rtg_path = project_root / "rtg"
    if rtg_path.exists():
        return str(rtg_path)

    # Try libexec directory (legacy location)
    rtg_path = project_root / "libexec" / "rtg"
    if rtg_path.exists():
        return str(rtg_path)

    # Check if rtg is available in PATH
    rtg_in_path = shutil.which("rtg")
    if rtg_in_path:
        return rtg_in_path

    # Fallback to system RTG, but let the user know it's not found
    print("WARNING: RTG tools not found. Integration tests requiring RTG may fail.")
    return "rtg"


@pytest.fixture(scope="session")
def rtg_executable():
    """Provide RTG executable path for tests."""
    rtg = get_rtg_path()
    # If we resolved to plain "rtg" ensure it exists in PATH
    if rtg == "rtg" and shutil.which("rtg") is None:
        pytest.skip(
            "RTG tools not found. Set RTG or RTGTOOLS_PATH to run integration tests"
        )
    return rtg


@pytest.fixture(scope="session")
def project_root():
    """Provide project root directory."""
    return Path(__file__).parent.parent


@pytest.fixture(scope="session")
def example_dir(project_root):
    """Provide example data directory."""
    candidate = project_root / "tests" / "data" / "example"
    if candidate.exists():
        return candidate
    return project_root / "example"


@pytest.fixture(scope="session")
def reference_file(example_dir):
    """Provide reference file for tests."""
    # Check environment variable first
    if "HGREF" in os.environ:
        return os.environ["HGREF"]

    # Use chr21.fa from example directory
    ref_path = example_dir / "chr21.fa"
    if ref_path.exists():
        return str(ref_path)

    # Fallback to example.fa
    ref_path = example_dir / "example.fa"
    if ref_path.exists():
        return str(ref_path)

    return None


@pytest.fixture
def sample_vcf_files(tmp_path):
    """
    Create sample VCF files for testing.

    Returns:
        tuple: (truth_vcf_path, query_vcf_path)
    """
    # Create a simple truth VCF
    truth_vcf = tmp_path / "truth.vcf"
    with open(truth_vcf, "w") as f:
        f.write(
            """##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##contig=<ID=chr1,length=248956422>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE
chr1	100	.	A	T	.	PASS	.	GT	0/1
chr1	200	.	G	C	.	PASS	.	GT	1/1
chr1	300	.	C	G	.	PASS	.	GT	0/1
"""
        )

    # Create a simple query VCF with one error
    query_vcf = tmp_path / "query.vcf"
    with open(query_vcf, "w") as f:
        f.write(
            """##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##contig=<ID=chr1,length=248956422>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE
chr1	100	.	A	T	.	PASS	.	GT	0/1
chr1	200	.	G	C	.	PASS	.	GT	0/1
chr1	400	.	T	A	.	PASS	.	GT	0/1
"""
        )

    return str(truth_vcf), str(query_vcf)


@pytest.fixture
def sample_reference(tmp_path):
    """
    Create a sample reference FASTA file for testing.

    Returns:
        str: Path to the reference FASTA file
    """
    # Create a simple reference FASTA
    ref_fasta = tmp_path / "reference.fa"
    with open(ref_fasta, "w") as f:
        f.write(
            """>chr1
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
"""
        )

    # Create a simple FASTA index (fai)
    with open(str(ref_fasta) + ".fai", "w") as f:
        f.write("chr1\t120\t6\t60\t61\n")

    return str(ref_fasta)


@pytest.fixture
def sample_bed_file(tmp_path):
    """
    Create a sample BED file for testing.

    Returns:
        str: Path to the BED file
    """
    # Create a simple BED file
    bed_file = tmp_path / "regions.bed"
    with open(bed_file, "w") as f:
        f.write(
            """chr1	50	150
chr1	175	225
chr1	275	325
"""
        )

    return str(bed_file)


@pytest.fixture
def mock_rtg_tools(tmp_path):
    """
    Create a mock RTG Tools directory for testing.

    Returns:
        str: Path to the mock RTG Tools directory
    """
    # Create a mock RTG Tools directory
    rtg_dir = tmp_path / "rtg"
    rtg_dir.mkdir()

    # Create a mock RTG executable
    rtg_exe = rtg_dir / "rtg"
    with open(rtg_exe, "w") as f:
        f.write(
            """#!/bin/sh
echo "RTG Tools mock implementation"
echo "Command: $@"
exit 0
"""
        )

    # Make the mock executable executable
    os.chmod(rtg_exe, 0o755)

    return str(rtg_dir)
