from pathlib import Path

import pysam
import pytest

from hap_py.utils import multimerge


def create_simple_vcf(path: Path, sample: str) -> None:
    header = pysam.VariantHeader()
    header.add_meta("fileformat", "VCFv4.2")
    header.contigs.add("chr1")
    header.add_sample(sample)
    header.formats.add("GT", number=1, type="String", description="Genotype")

    with pysam.VariantFile(str(path), "w", header=header) as out:
        rec = header.new_record(contig="chr1", start=0, stop=1, alleles=("A", "C"))
        rec.samples[sample]["GT"] = (0, 1)
        out.write(rec)


def create_alt_vcf(path: Path, sample: str, alt: str) -> None:
    """Create VCF with a configurable ALT allele."""
    header = pysam.VariantHeader()
    header.add_meta("fileformat", "VCFv4.2")
    header.contigs.add("chr1")
    header.add_sample(sample)
    header.formats.add("GT", number=1, type="String", description="Genotype")

    with pysam.VariantFile(str(path), "w", header=header) as out:
        rec = header.new_record(contig="chr1", start=0, stop=1, alleles=("A", alt))
        rec.samples[sample]["GT"] = (0, 1)
        out.write(rec)


def test_parse_args() -> None:
    ns = multimerge._parse_args(
        [
            "a.vcf:S1",
            "b.vcf:S2",
            "-o",
            "out.vcf",
            "-r",
            "ref.fa",
        ]
    )
    assert ns.inputs == ["a.vcf:S1", "b.vcf:S2"]
    assert ns.output == "out.vcf"
    assert ns.reference == "ref.fa"


def test_parse_args_missing_output() -> None:
    with pytest.raises(SystemExit):
        multimerge._parse_args(["a.vcf:S1", "-r", "ref.fa"])


def test_merge_records(tmp_path: Path) -> None:
    v1 = tmp_path / "a.vcf"
    v2 = tmp_path / "b.vcf"
    create_simple_vcf(v1, "S1")
    create_simple_vcf(v2, "S2")
    out_vcf = tmp_path / "merged.vcf"

    multimerge._merge_records([(v1, "S1"), (v2, "S2")], out_vcf, tmp_path / "ref.fa")

    with pysam.VariantFile(str(out_vcf)) as vf:
        records = list(vf.fetch())
        assert len(records) == 1
        assert list(vf.header.samples) == ["S1", "S2"]
        gt1 = records[0].samples["S1"]["GT"]
        gt2 = records[0].samples["S2"]["GT"]
        assert gt1 == (0, 1) and gt2 == (0, 1)


def test_merge_records_union_alt(tmp_path: Path) -> None:
    v1 = tmp_path / "a.vcf"
    v2 = tmp_path / "b.vcf"
    create_alt_vcf(v1, "S1", "T")
    create_alt_vcf(v2, "S2", "G")
    out_vcf = tmp_path / "merged.vcf"

    multimerge._merge_records([(v1, "S1"), (v2, "S2")], out_vcf, tmp_path / "ref.fa")

    with pysam.VariantFile(str(out_vcf)) as vf:
        record = next(vf.fetch())
        assert set(record.alts) == {"T", "G"}


def test_main(tmp_path: Path) -> None:
    v1 = tmp_path / "a.vcf"
    v2 = tmp_path / "b.vcf"
    create_simple_vcf(v1, "S1")
    create_simple_vcf(v2, "S2")
    out_vcf = tmp_path / "merged.vcf"
    rc = multimerge.main(
        [
            f"{v1}:S1",
            f"{v2}:S2",
            "-o",
            str(out_vcf),
            "-r",
            str(tmp_path / "ref.fa"),
        ]
    )
    assert rc == 0
    assert out_vcf.exists()


def test_main_default_samples(tmp_path: Path) -> None:
    v1 = tmp_path / "sample1.vcf"
    v2 = tmp_path / "sample2.vcf"
    create_simple_vcf(v1, "S1")
    create_simple_vcf(v2, "S2")
    out_vcf = tmp_path / "merged.vcf"
    rc = multimerge.main(
        [
            str(v1),
            str(v2),
            "-o",
            str(out_vcf),
            "-r",
            str(tmp_path / "ref.fa"),
        ]
    )
    assert rc == 0
    with pysam.VariantFile(str(out_vcf)) as vf:
        assert list(vf.header.samples) == ["S1", "sample1", "sample2"]
