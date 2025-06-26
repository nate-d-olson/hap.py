from src.hap_py.pre import hasChrPrefix


def fixChrPrefix(chromosomes):
    """Helper function to add or remove 'chr' prefix to chromosome names."""
    if all(chrom.startswith("chr") for chrom in chromosomes):
        # Remove 'chr' prefix
        return [
            chrom[3:] if chrom.startswith("chr") else chrom for chrom in chromosomes
        ]
    else:
        # Add 'chr' prefix
        return [
            f"chr{chrom}" if not chrom.startswith("chr") else chrom
            for chrom in chromosomes
        ]


def test_hasChrPrefix_with_prefix():
    chroms = ["chr1", "chr2", "chrX", "chrY", "chrM"]
    assert hasChrPrefix(chroms) is True


def test_hasChrPrefix_without_prefix():
    chroms = ["1", "2", "X", "Y", "MT"]
    assert hasChrPrefix(chroms) is False


def test_hasChrPrefix_mixed():
    chroms = ["chr1", "2", "chrX", "Y", "MT"]
    # When counts are equal, returns None
    assert hasChrPrefix(chroms) is False


def test_fixChrPrefix_add():
    chroms = ["1", "2", "X", "Y", "MT"]
    fixed = fixChrPrefix(chroms)
    assert all(chrom.startswith("chr") for chrom in fixed)
    assert fixed == ["chr1", "chr2", "chrX", "chrY", "chrMT"]


def test_fixChrPrefix_remove():
    chroms = ["chr1", "chr2", "chrX", "chrY", "chrM"]
    fixed = fixChrPrefix(chroms)
    assert all(not chrom.startswith("chr") for chrom in fixed)
    assert fixed == ["1", "2", "X", "Y", "M"]
    assert all(not chrom.startswith("chr") for chrom in fixed)
    assert fixed == ["1", "2", "X", "Y", "M"]
