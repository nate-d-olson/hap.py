from pathlib import Path

import pysam

from hap_py.haplo import gvcf2bed


def create_vcf(path: Path) -> None:
    header = pysam.VariantHeader()
    header.add_meta("fileformat", "VCFv4.2")
    header.contigs.add("chr1")
    header.add_sample("SAMPLE")
    header.formats.add("GT", number=1, type="String", description="Genotype")
    header.filters.add("FAIL", None, None, "Fail filter")

    with pysam.VariantFile(str(path), "w", header=header) as out:
        rec1 = header.new_record(contig="chr1", start=0, stop=1, alleles=("A", "C"))
        rec1.filter.add("PASS")
        rec1.samples["SAMPLE"]["GT"] = (0, 1)
        out.write(rec1)

        rec2 = header.new_record(contig="chr1", start=10, stop=11, alleles=("G", "T"))
        rec2.filter.add("FAIL")
        rec2.samples["SAMPLE"]["GT"] = (0, 1)
        out.write(rec2)

def create_overlap_vcf(path: Path) -> None:
    """Create a VCF with adjacent confident calls."""
    header = pysam.VariantHeader()
    header.add_meta("fileformat", "VCFv4.2")
    header.contigs.add("chr1")
    header.add_sample("SAMPLE")
    header.formats.add("GT", number=1, type="String", description="Genotype")

    with pysam.VariantFile(str(path), "w", header=header) as out:
        r1 = header.new_record(contig="chr1", start=0, stop=1, alleles=("A", "C"))
        r1.filter.add("PASS")
        r1.samples["SAMPLE"]["GT"] = (0, 1)
        out.write(r1)

        r2 = header.new_record(contig="chr1", start=1, stop=2, alleles=("G", "T"))
        r2.filter.add("PASS")
        r2.samples["SAMPLE"]["GT"] = (0, 1)
        out.write(r2)


def create_uncertain_vcf(path: Path) -> None:
    """VCF with non-PASS and missing genotype records."""
    header = pysam.VariantHeader()
    header.add_meta("fileformat", "VCFv4.2")
    header.contigs.add("chr1")
    header.add_sample("SAMPLE")
    header.formats.add("GT", number=1, type="String", description="Genotype")
    header.filters.add("LowQual", None, None, "Low quality")

    with pysam.VariantFile(str(path), "w", header=header) as out:
        fail_rec = header.new_record(contig="chr1", start=0, stop=1, alleles=("A", "C"))
        fail_rec.filter.add("LowQual")
        fail_rec.samples["SAMPLE"]["GT"] = (0, 1)
        out.write(fail_rec)

        nocall = header.new_record(contig="chr1", start=5, stop=6, alleles=("G", "T"))
        nocall.filter.add("PASS")
        nocall.samples["SAMPLE"]["GT"] = (None, None)
        out.write(nocall)


def test_gvcf2bed_basic(tmp_path: Path) -> None:
    vcf = tmp_path / "test.vcf"
    create_vcf(vcf)
    ref = Path("tests/data/common/test.fa")
    bed_path = gvcf2bed.gvcf2bed(str(vcf), str(ref), scratch_prefix=str(tmp_path))
    with open(bed_path, encoding="utf-8") as fh:
        lines = [line.strip() for line in fh.readlines()]
    assert lines == ["chr1\t0\t1"]


def test_gvcf2bed_regions_filter(tmp_path: Path) -> None:
    vcf = tmp_path / "test.vcf"
    create_vcf(vcf)
    regions = tmp_path / "regions.bed"
    regions.write_text("chr1\t5\t6\n")
    ref = Path("tests/data/common/test.fa")
    bed_path = gvcf2bed.gvcf2bed(
        str(vcf), str(ref), regions=str(regions), scratch_prefix=str(tmp_path)
    )
    with open(bed_path, encoding="utf-8") as fh:
        content = fh.read().strip()
    assert content == ""


def test_gvcf2bed_merge_adjacent(tmp_path: Path) -> None:
    vcf = tmp_path / "merge.vcf"
    create_overlap_vcf(vcf)
    ref = Path("tests/data/common/test.fa")
    bed_path = gvcf2bed.gvcf2bed(str(vcf), str(ref), scratch_prefix=str(tmp_path))
    with open(bed_path, encoding="utf-8") as fh:
        lines = [line.strip() for line in fh.readlines()]
    assert lines == ["chr1\t0\t2"]


def test_gvcf2bed_skip_uncertain(tmp_path: Path) -> None:
    vcf = tmp_path / "uncertain.vcf"
    create_uncertain_vcf(vcf)
    ref = Path("tests/data/common/test.fa")
    bed_path = gvcf2bed.gvcf2bed(str(vcf), str(ref), scratch_prefix=str(tmp_path))
    with open(bed_path, encoding="utf-8") as fh:
        content = fh.read().strip()
    assert content == ""
