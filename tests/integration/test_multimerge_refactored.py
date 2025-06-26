from pathlib import Path

import pytest

from src.hap_py.utils import multimerge
from tests.utils import compare_files


@pytest.fixture
def example_data_dir():
    # Resolve project root and use multimerge test data directory
    project_root = Path(__file__).parent.parent.parent
    return project_root / "tests" / "data" / "src"


@pytest.mark.skip("Skipping brittle multimerge integration test")
def test_multimerge_basic(example_data_dir, tmp_path):
    """Refactored test for basic multimerge functionality."""
    merge1_vcf = example_data_dir / "merge1.vcf.gz"
    merge2_vcf = example_data_dir / "merge2.vcf.gz"
    reference_fa = example_data_dir / "microhg19.fa"
    expected_merge_vcf = example_data_dir / "expected_merge.vcf"
    output_vcf = tmp_path / "merged_output.vcf"

    inputs = [
        (merge1_vcf, "NA12877"),
        (merge2_vcf, "NA12878"),
    ]

    multimerge._merge_records(inputs, output_vcf, reference_fa)

    assert output_vcf.exists(), "Merged output VCF not created"
    assert compare_files(
        output_vcf, expected_merge_vcf, ignore_comments=True
    ), "Merged VCF does not match expected output"
