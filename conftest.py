"""
Pytest configuration for hap.py test suite.
"""

import os
import shutil
import sys
import tempfile
from pathlib import Path

import pytest

# Add the repository's ``src`` directory to ``sys.path`` so tests import the
# ``hap_py`` package from the working tree rather than any installed version.
repo_root = Path(__file__).resolve().parent
src_path = repo_root / "src"
sys.path.insert(0, str(src_path))


@pytest.fixture(scope="session")
def rtg_tools_path():
    """
    Robust RTG tools finder fixture.
    Returns the path to RTG tools executable for tests.
    """
    # List of potential RTG locations in order of preference
    repo_root = Path(__file__).parent
    potential_paths = [
        # Local build directory paths
        repo_root / "external" / "rtg-tools-3.12.1" / "rtg",
        repo_root / "build" / "external" / "rtg-tools" / "rtg",
        repo_root / "external" / "rtg-tools" / "rtg",
        # Check if rtg is in PATH
        shutil.which("rtg"),
        # Other common locations
        "/usr/local/bin/rtg",
        "/opt/rtg/rtg",
    ]

    for path in potential_paths:
        if path is not None:
            path_str = str(path)
            if Path(path_str).exists() and os.access(path_str, os.X_OK):
                return path_str

    # If no RTG found, skip tests that require it
    pytest.skip("RTG tools not found. Please install RTG tools or check the path.")


@pytest.fixture(scope="session")
def temp_dir_with_rtg(rtg_tools_path):
    """Create a temporary directory and ensure RTG tools are available."""
    temp_dir = tempfile.mkdtemp(prefix="happy_rtg_test_")

    # Create a symbolic link to RTG in the temp directory for tests
    rtg_link = os.path.join(temp_dir, "rtg")
    try:
        os.symlink(rtg_tools_path, rtg_link)
        yield temp_dir, rtg_tools_path
    finally:
        shutil.rmtree(temp_dir)


@pytest.fixture(scope="session")
def example_data_dir():
    """Return the path to the example data directory."""
    return os.path.abspath(os.path.join(os.path.dirname(__file__), "example"))


@pytest.fixture(scope="function")
def temp_dir():
    """Create a temporary directory for test data that is cleaned up after the test."""
    temp_dir = tempfile.mkdtemp(prefix="happy_test_")
    yield temp_dir
    shutil.rmtree(temp_dir)


@pytest.fixture(scope="session")
def reference_file(example_data_dir):
    """Return the path to the reference file used in tests."""
    return os.path.join(example_data_dir, "chr21.fa")


@pytest.fixture(scope="session")
def truth_vcf(example_data_dir):
    """Return the path to the truth VCF file used in tests."""
    return os.path.join(example_data_dir, "performance.vcf.gz")


@pytest.fixture(scope="session")
def query_vcf(example_data_dir):
    """Return the path to the query VCF file used in tests."""
    return os.path.join(example_data_dir, "hc.vcf.gz")


@pytest.fixture(scope="session")
def confident_regions(example_data_dir):
    """Return the path to the confident regions BED file used in tests."""
    return os.path.join(example_data_dir, "performance.confident.bed.gz")
