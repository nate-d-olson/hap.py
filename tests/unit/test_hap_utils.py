import pytest

from hap_py.hap import find_common_chromosomes, generate_summary


@pytest.mark.unit
def test_generate_summary_correct_format():
    """Summary string should match exact required format."""
    total = 10
    meeting = 5
    failing = 2
    coverage = 75.5
    summary = generate_summary(total, meeting, failing, coverage)
    expected = (
        "Numeric Chrs Summary Report:\n"
        "• Total variants processed: 10\n"
        "• Variants meeting criteria: 5\n"
        "• Variants failing criteria: 2\n"
        "• Overall coverage: 75.5%"
    )
    assert summary == expected


@pytest.mark.unit
def test_find_common_chromosomes_with_chr_prefix_and_mt():
    """Should normalize 'chr' prefixes and handle mitochondrial names."""
    reference_contigs = {"1", "2", "chr3", "MT", "chrX"}
    vcf_chromosomes = {"chr1", "3", "m", "chrY"}
    common = find_common_chromosomes(reference_contigs, vcf_chromosomes)
    # Intersection normalized keys: "1","3","mt"
    # VCF naming for those keys: "chr1","3","m"
    assert set(common) == {"chr1", "3", "m"}


@pytest.mark.unit
def test_find_common_chromosomes_empty_when_no_match():
    """Should return empty list when no common normalized names."""
    reference_contigs = {"chrA", "chrB"}
    vcf_chromosomes = {"1", "2"}
    common = find_common_chromosomes(reference_contigs, vcf_chromosomes)
    assert common == []
    common = find_common_chromosomes(reference_contigs, vcf_chromosomes)
    assert common == []
