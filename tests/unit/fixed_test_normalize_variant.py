import tempfile

import pytest

pytest.importorskip("pysam")
import pysam

from hap_py.haplo.python_preprocess import PreprocessEngine


def test_normalize_variant():
    """Test variant normalization."""
    # Create dummy VCF and FASTA files for initialization
    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".vcf", delete=False
    ) as dummy_vcf_file, tempfile.NamedTemporaryFile(
        mode="w", suffix=".fa", delete=False
    ) as dummy_fa_file:
        dummy_vcf_path = dummy_vcf_file.name
        dummy_fa_path = dummy_fa_file.name

        # Write valid VCF content
        dummy_vcf_file.write(
            """##fileformat=VCFv4.2
##reference=file://test.fa
##contig=<ID=chr1,length=200>
##FILTER=<ID=PASS,Description="All filters passed">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample1
chr1\t100\t.\tA\tG\t.\tPASS\t.\tGT\t0/1
"""
        )
        dummy_vcf_file.flush()

        # Write minimal valid content to FASTA to avoid pysam errors
        # and ensure it's indexable by pysam.faidx
        dummy_fa_file.write(">chr1\n")
        dummy_fa_file.write("A" * 60 + "\n")  # Standard FASTA line length
        dummy_fa_file.write("C" * 60 + "\n")
        dummy_fa_file.flush()
        pysam.faidx(dummy_fa_path)  # Create .fai index

    engine = PreprocessEngine(input_vcf=dummy_vcf_path, reference_fasta=dummy_fa_path)

    # Test normalization of SNV - should not change
    pos, ref, alt = engine.normalize_variant("chr1", 100, "A", "G")
    assert pos == 100
    assert ref == "A"
    assert alt == "G"

    # Test normalization with common prefix and suffix
    pos, ref, alt = engine.normalize_variant("chr1", 100, "ATCG", "ATTG")
    assert pos == 101
    assert ref == "TC"
    assert alt == "TT"

    # Test normalization with only common prefix
    pos, ref, alt = engine.normalize_variant("chr1", 100, "ATCG", "ATTT")
    assert pos == 102
    assert ref == "CG"
    assert alt == "TT"

    # Test normalization with only common suffix
    pos, ref, alt = engine.normalize_variant("chr1", 100, "ATCG", "TTCG")
    assert pos == 100
    assert ref == "AT"
    assert alt == "TT"
    assert pos == 100
    assert ref == "AT"
    assert alt == "TT"
    assert pos == 100
    assert ref == "AT"
    assert alt == "TT"
    assert pos == 100
    assert ref == "AT"
    assert alt == "TT"
