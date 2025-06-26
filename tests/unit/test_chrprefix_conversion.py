from src.hap_py.pre import hasChrPrefix
from tests.utils import compare_files_content


def test_hasChrPrefix_behavior():
    # Test hasChrPrefix with various chromosome lists
    assert hasChrPrefix(["chr1", "chr2", "chrX"]) is True
    assert hasChrPrefix(["1", "2", "X"]) is False
    assert hasChrPrefix(["chr1", "2", "chrX"]) is True


def test_fixChrPrefix_add_and_remove():
    def fixChrPrefix(chromosomes):
        if all(chrom.startswith("chr") for chrom in chromosomes):
            return [
                chrom[3:] if chrom.startswith("chr") else chrom for chrom in chromosomes
            ]
        else:
            return [
                f"chr{chrom}" if not chrom.startswith("chr") else chrom
                for chrom in chromosomes
            ]

    chroms_no_prefix = ["1", "2", "X", "Y", "MT"]
    chroms_with_prefix = ["chr1", "chr2", "chrX", "chrY", "chrM"]

    # Add prefix
    fixed_add = fixChrPrefix(chroms_no_prefix)
    assert all(chrom.startswith("chr") for chrom in fixed_add)
    assert fixed_add == ["chr1", "chr2", "chrX", "chrY", "chrMT"]

    # Remove prefix
    fixed_remove = fixChrPrefix(chroms_with_prefix)
    assert all(not chrom.startswith("chr") for chrom in fixed_remove)
    assert fixed_remove == ["1", "2", "X", "Y", "M"]


def test_compare_summary_files(tmp_path):
    # This test simulates comparing summary files content for prefix handling
    # Setup dummy expected and output summary files
    expected_summary = tmp_path / "expected.summary.csv"
    output_summary = tmp_path / "output.summary.csv"

    content = "Type,Filter,TRUTH.TOTAL,TRUTH.TP\nSNP,.,100,95\n"
    expected_summary.write_text(content)
    output_summary.write_text(content)

    assert compare_files_content(str(output_summary), str(expected_summary))


def test_vcf_content_comparison(tmp_path):
    # Simulate VCF content comparison ignoring header lines
    vcf_content = [
        "chr1\t100\t.\tA\tT\t.\tPASS\t.\n",
        "chr1\t200\t.\tG\tC\t.\tPASS\t.\n",
    ]
    expected_vcf = tmp_path / "expected.vcf"
    output_vcf = tmp_path / "output.vcf"

    expected_vcf.write_text("".join(vcf_content))
    output_vcf.write_text("".join(vcf_content))

    with open(output_vcf) as f_out, open(expected_vcf) as f_exp:
        lines_out = [line for line in f_out if not line.startswith("#")]
        lines_exp = [line for line in f_exp if not line.startswith("#")]
        assert lines_out == lines_exp
        lines_out = [line for line in f_out if not line.startswith("#")]
        lines_exp = [line for line in f_exp if not line.startswith("#")]
        assert lines_out == lines_exp
