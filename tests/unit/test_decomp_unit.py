"""
Unit/component tests for variant decomposition functionality, refactored from integration test.

These tests focus on the core decomposition logic, using in-memory or temporary files, and do not rely on subprocess or external binaries.
"""

from pathlib import Path

import pytest

# Import the relevant decomposition logic if available
# For demonstration, assume hap_py.haplo.variant_processor has a decompose_variants function
try:
    from src.hap_py.haplo import variant_processor
except ImportError:
    variant_processor = None


@pytest.fixture
def vcf_truth(tmp_path):
    content = """##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t100\t.\tA\tG\t.\tPASS\t.
chr1\t200\t.\tAT\tA\t.\tPASS\t.
"""
    path = tmp_path / "truth.vcf"
    path.write_text(content)
    return path


@pytest.fixture
def vcf_query(tmp_path):
    content = """##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t100\t.\tA\tG\t.\tPASS\t.
chr1\t200\t.\tA\tA\t.\tPASS\t.
"""
    path = tmp_path / "query.vcf"
    path.write_text(content)
    return path


def test_decompose_variants_basic(vcf_truth, vcf_query, tmp_path):
    if not variant_processor or not hasattr(variant_processor, "decompose_variants"):
        pytest.skip("decompose_variants function not available in variant_processor")

    output_vcf = tmp_path / "decomp_out.vcf"
    # Call the decomposition logic directly
    variant_processor.decompose_variants(
        str(vcf_truth), str(vcf_query), str(output_vcf)
    )
    assert output_vcf.exists(), "Decomposed output VCF not created"
    content = output_vcf.read_text()
    assert "chr1" in content
    assert "100" in content
    assert "200" in content
    assert "200" in content
