from pathlib import Path

import pytest

import hap_py.pre as premod
from hap_py.pre import preprocess


@pytest.fixture
def basic_vcf(tmp_path):
    vcf = tmp_path / "in.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
        "chr1\t1\t.\tA\tT\t.\tPASS\t.\tGT\t0/1\n"
    )
    (tmp_path / "in.vcf.fai").write_text("chr1\t1\t0\t1\t2\n")
    return str(vcf)


@pytest.fixture
def basic_ref(tmp_path):
    ref = tmp_path / "ref.fa"
    ref.write_text(">chr1\nA\n")
    (tmp_path / "ref.fa.fai").write_text("chr1\t1\t0\t1\t2\n")
    return str(ref)


def test_preprocess_decompose_only(tmp_path, monkeypatch, basic_vcf, basic_ref):
    # Stub heavy operations to avoid external calls
    monkeypatch.setattr(premod, "runBcftools", lambda *args, **kwargs: None)
    monkeypatch.setattr(
        premod, "preprocessVCF", lambda *args, **kwargs: Path(args[1]).write_text("")
    )
    monkeypatch.setattr(
        premod.partialcredit,
        "partialCredit",
        lambda vtf, out, *args, **kwargs: Path(out).write_text(Path(vtf).read_text()),
    )
    monkeypatch.setattr(
        premod.vcfextract,
        "extractHeadersJSON",
        lambda x: {"tabix": {"chromosomes": ["chr1"]}, "fields": []},
    )
    monkeypatch.setattr(premod, "fastaContigLengths", lambda x: {"chr1": 1})

    output_vcf = tmp_path / "out.vcf"
    gender = preprocess(
        basic_vcf,
        str(output_vcf),
        basic_ref,
        locations=None,
        filters=None,
        fixchr=False,
        regions=None,
        targets=None,
        leftshift=False,
        decompose=True,
        bcftools_norm=False,
        threads=1,
        gender=None,
        somatic_allele_conversion=False,
    )

    assert output_vcf.exists(), "Output VCF was not created"
    assert gender is None or isinstance(gender, str)
