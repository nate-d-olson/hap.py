import os
import tempfile
from pathlib import Path

import pytest

# Ensure dummy executables for bgzip and tabix so tools.init() succeeds
tmp_tool_dir = tempfile.mkdtemp()
for tool in ("bgzip", "tabix"):
    tool_path = Path(tmp_tool_dir) / tool
    with open(tool_path, "w", encoding="utf-8") as t:
        t.write("#!/bin/sh\nexit 0\n")
    os.chmod(tool_path, 0o755)
os.environ["PATH"] = f"{tmp_tool_dir}:{os.environ.get('PATH','')}"

from hap_py.haplo.python_quantify import QuantifyEngine


def _write_vcf(
    path: Path, include_all_fields: bool = True, include_samples: bool = True
):
    header = [
        "##fileformat=VCFv4.2",
        '##INFO=<ID=BS,Number=1,Type=Integer,Description="Benchmarking superlocus ID">',
    ]
    format_fields = [
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '##FORMAT=<ID=BD,Number=1,Type=String,Description="Decision">',
        '##FORMAT=<ID=BK,Number=1,Type=String,Description="Decision subtype">',
        '##FORMAT=<ID=BI,Number=1,Type=String,Description="Additional info">',
        '##FORMAT=<ID=QQ,Number=1,Type=Float,Description="Quality">',
        '##FORMAT=<ID=BVT,Number=1,Type=String,Description="Variant type">',
        '##FORMAT=<ID=BLT,Number=1,Type=String,Description="Location type">',
    ]
    if include_all_fields:
        header.extend(format_fields)
    else:
        header.append(format_fields[0])  # only GT
    header_line = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"
    if include_samples:
        header_line += "\tTRUTH\tQUERY"
    else:
        header_line += "\tSAMPLE"
    header.append(header_line)

    body = (
        "chr1\t100\t.\tA\tT\t50\tPASS\tBS=1"
        "\tGT:BD:BK:BI:QQ:BVT:BLT"
        "\t0/1:TP:gm:info:100:SNP:het"
        "\t0/1:TP:gm:info:100:SNP:het\n"
    )
    with open(path, "w", encoding="utf-8") as f:
        for line in header:
            f.write(line + "\n")
        f.write(body)


def test_ga4gh_validation_success(tmp_path):
    vcf = tmp_path / "valid.vcf"
    _write_vcf(vcf)
    engine = QuantifyEngine(str(vcf), str(vcf), quantify_method="ga4gh")
    engine.quantify()


def test_ga4gh_validation_failure(tmp_path):
    vcf = tmp_path / "invalid.vcf"
    _write_vcf(vcf, include_all_fields=False, include_samples=False)
    with pytest.raises(ValueError):
        QuantifyEngine(str(vcf), str(vcf), quantify_method="ga4gh")
