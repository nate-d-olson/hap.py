"""VCF header generation utilities."""

from datetime import date
from typing import List, Union, TextIO, BinaryIO, Optional


def write_vcf_header(
    filelike: Union[TextIO, BinaryIO],
    extrainfo: Union[str, List[str]] = "",
    chrprefix: str = "chr",
    version: Optional[str] = None,
) -> None:
    """Write a VCF header with standard information.

    Args:
        filelike: File-like object to write header to
        extrainfo: Additional header lines (string or list of strings)
        chrprefix: Prefix to use for chromosome names
        version: Version string to include in header
    """
    header = [
        "CHROM",
        "POS",
        "ID",
        "REF",
        "ALT",
        "QUAL",
        "FILTER",
        "INFO",
        "FORMAT",
        "SIMPLE",
    ]

    infos = [
        "##fileformat=VCFv4.1",
        "##reference=hg19",
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    ]

    if extrainfo:
        if isinstance(extrainfo, list):
            infos.extend(extrainfo)
        else:
            infos.extend(extrainfo.split("\n"))

    # Standard human chromosomes with lengths
    contigs = [
        ["1", "249250621"],
        ["2", "243199373"],
        ["3", "198022430"],
        ["4", "191154276"],
        ["5", "180915260"],
        ["6", "171115067"],
        ["7", "159138663"],
        ["8", "146364022"],
        ["9", "141213431"],
        ["10", "135534747"],
        ["11", "135006516"],
        ["12", "133851895"],
        ["13", "115169878"],
        ["14", "107349540"],
        ["15", "102531392"],
        ["16", "90354753"],
        ["17", "81195210"],
        ["18", "78077248"],
        ["19", "59128983"],
        ["20", "63025520"],
        ["21", "48129895"],
        ["22", "51304566"],
        ["X", "155270560"],
    ]

    meta = [
        "##fileDate=" + date.today().isoformat(),
        "##source=HaploCompare",
    ]

    if version:
        meta.append(f"##source_version={version}")

    # Write all header components
    for line in infos:
        write_line(filelike, line + "\n")
    for chrom, length in contigs:
        write_line(filelike, f"##contig=<ID={chrprefix}{chrom},length={length}>\n")
    for line in meta:
        write_line(filelike, line + "\n")
    write_line(filelike, "#" + "\t".join(header) + "\n")


def write_line(filelike: Union[TextIO, BinaryIO], text: str) -> None:
    """Write a line to a file-like object, handling both text and binary modes.

    Args:
        filelike: File-like object to write to
        text: Text to write
    """
    if hasattr(filelike, "mode") and "b" in filelike.mode:
        filelike.write(text.encode("utf-8"))
    else:
        filelike.write(text)
