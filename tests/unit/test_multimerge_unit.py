from pathlib import Path

import pytest

from src.hap_py.utils import multimerge


@pytest.fixture
def vcf_content_1():
    return """##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1
chr1\t100\t.\tA\tG\t.\tPASS\t.\tGT\t0/1
"""


@pytest.fixture
def vcf_content_2():
    return """##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample2
chr1\t100\t.\tA\tT\t.\tPASS\t.\tGT\t0/1
"""


@pytest.fixture
def reference_fasta(tmp_path):
    # Create a minimal reference fasta file
    fasta_path = tmp_path / "ref.fa"
    fasta_path.write_text(">chr1\n" + "A" * 200 + "\n")
    return fasta_path


def write_vcf(content: str, tmp_path: Path, filename: str) -> Path:
    path = tmp_path / filename
    path.write_text(content)
    return path


def test_multimerge_unit(vcf_content_1, vcf_content_2, reference_fasta, tmp_path):
    vcf1_path = write_vcf(vcf_content_1, tmp_path, "vcf1.vcf")
    vcf2_path = write_vcf(vcf_content_2, tmp_path, "vcf2.vcf")
    output_vcf = tmp_path / "merged.vcf"

    inputs = [
        (vcf1_path, "sample1"),
        (vcf2_path, "sample2"),
    ]

    multimerge._merge_records(inputs, output_vcf, reference_fasta)

    assert output_vcf.exists(), "Merged output VCF not created"

    # Read output and check for expected merged content
    merged_content = output_vcf.read_text()
    assert "sample1" in merged_content
    assert "sample2" in merged_content
    assert "chr1" in merged_content
    assert "100" in merged_content


def test_multimerge_import_error(tmp_path, reference_fasta):
    """_merge_records should raise an error for malformed VCF input."""
    # Create malformed VCF content
    vcf_content = """##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1
chr1\tBAD\t.\tA\tG\t.\tPASS\t.\tGT\t0/1
"""
    vcf_path = write_vcf(vcf_content, tmp_path, "import_error.vcf")
    output_vcf = tmp_path / "import_result.vcf"
    inputs = [(vcf_path, "sample1")]

    import pytest

    with pytest.raises(Exception) as excinfo:
        multimerge._merge_records(inputs, output_vcf, reference_fasta)
    err = str(excinfo.value)
    assert any(sub in err for sub in ("BAD", "invalid", "error"))
