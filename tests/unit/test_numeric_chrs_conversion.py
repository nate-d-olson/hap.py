from unittest.mock import MagicMock, patch

# Assuming the core logic for chromosome prefix detection and processing is in src.hap_py.pre
from src.hap_py.pre import hasChrPrefix


def test_numeric_chrs_logic(tmp_path):
    # Setup dummy file paths
    src_data_dir = tmp_path / "numeric_chrs"
    src_data_dir.mkdir()

    truth_vcf = src_data_dir / "truth.vcf.gz"
    query_vcf = src_data_dir / "query.vcf.gz"
    fp_bed = src_data_dir / "fp.bed"
    reference = src_data_dir / "test.fa"
    expected_vcf = src_data_dir / "expected.vcf"
    expected_summary = src_data_dir / "expected.summary.csv"

    # Create dummy files
    for f in [truth_vcf, query_vcf, fp_bed, reference, expected_vcf, expected_summary]:
        f.write_text("dummy content")

    # Mock subprocess.run to simulate hap.py CLI call success
    with patch("subprocess.run") as mock_run:
        mock_result = MagicMock()
        mock_result.returncode = 0
        mock_run.return_value = mock_result

        # Simulate running hap.py CLI command
        cmd = [
            "hap.py",
            str(truth_vcf),
            str(query_vcf),
            "-f",
            str(fp_bed),
            "-o",
            str(tmp_path / "numeric_test"),
            "-X",
            "--reference",
            str(reference),
            "-V",
            "--force-interactive",
            "--engine-vcfeval-path",
            "/path/to/rtg",
        ]

        result = mock_run(cmd, capture_output=True, check=True)
        assert result.returncode == 0

    # Test hasChrPrefix function on sample chromosomes
    chroms = ["1", "2", "3"]
    assert hasChrPrefix(chroms) is False

    chroms_prefixed = ["chr1", "chr2", "chr3"]
    assert hasChrPrefix(chroms_prefixed) is True
    chroms_prefixed = ["chr1", "chr2", "chr3"]
    assert hasChrPrefix(chroms_prefixed) is True
