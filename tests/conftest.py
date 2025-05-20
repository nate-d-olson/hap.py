"""
Integration test fixtures for hap.py.

These fixtures provide test resources for integration tests of the hap.py toolkit.
"""

import os

import pytest


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
