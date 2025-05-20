"""
Pytest configuration for hap.py test suite.
"""

import os
import shutil
import sys
import tempfile

import pytest

# Add src/python to path for imports during tests
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "src", "python"))


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
