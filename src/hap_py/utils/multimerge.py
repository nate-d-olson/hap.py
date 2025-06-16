"""Simplified multimerge implementation.

This module provides a minimal Python implementation of the old C++ multimerge
utility. It currently performs a very basic merge of multiple VCF files using
pysam as a reference implementation. The goal is to provide enough
functionality for testing and compatibility with the previous interface.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pysam


def _parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser("multimerge")
    parser.add_argument("inputs", nargs="+", help="VCF files in the form file:sample")
    parser.add_argument("-o", "--output", required=True, help="Output VCF file")
    parser.add_argument("-r", "--reference", required=True, help="Reference FASTA")
    return parser.parse_args(argv)


def _read_vcf(path: Path) -> pysam.VariantFile:
    return pysam.VariantFile(str(path))


def _merge_records(
    inputs: list[tuple[Path, str]], out_vcf: Path, reference: Path
) -> None:
    header = None
    opened_vcfs: list[tuple[pysam.VariantFile, str]] = []
    records: dict[tuple[str, int, str], pysam.VariantRecord] = {}

    # Open all VCFs first to build the combined header.
    # TODO: The current implementation only preserves genotype (GT) fields and
    # ignores other annotations. A more complete merge should handle all INFO
    # and FORMAT fields consistently.
    for vcf_path, sample in inputs:
        vcf = _read_vcf(vcf_path)
        opened_vcfs.append((vcf, sample))
        if header is None:
            header = vcf.header.copy()
            if "GT" not in header.formats:
                header.formats.add(
                    "GT", number=1, type="String", description="Genotype"
                )

    for _, sample in opened_vcfs:
        if sample not in header.samples:
            header.add_sample(sample)

    for vcf, sample in opened_vcfs:
        # TODO: Streaming records would avoid constructing a potentially large
        # in-memory dictionary when merging many files.
        for rec in vcf.fetch():
            key = (rec.chrom, rec.pos, rec.ref)
            if key not in records:
                if rec.alts:
                    alleles = (rec.ref,) + tuple(rec.alts)
                else:
                    alleles = (rec.ref, ".")
                new_rec = header.new_record(
                    contig=rec.chrom,
                    start=rec.pos - 1,
                    stop=rec.stop,
                    id=rec.id,
                    alleles=alleles,
                )
                records[key] = new_rec
            else:
                new_rec = records[key]
                alts = {a for a in new_rec.alts or [] if a != "."}
                alts.update(a for a in (rec.alts or []) if a != ".")
                new_rec.alts = list(alts) if alts else ["."]
            gt = rec.samples[0].get("GT")
            records[key].samples[sample]["GT"] = gt
    with pysam.VariantFile(str(out_vcf), "w", header=header) as out:
        for rec in sorted(records.values(), key=lambda r: (r.contig, r.pos)):
            out.write(rec)


def main(argv: list[str] | None = None) -> int:
    args = _parse_args(argv or sys.argv[1:])
    inputs = []
    for spec in args.inputs:
        if ":" in spec:
            path, sample = spec.split(":", 1)
        else:
            path = spec
            sample = Path(path).stem
        inputs.append((Path(path), sample))
    _merge_records(inputs, Path(args.output), Path(args.reference))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
