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
    header: pysam.VariantHeader | None = None
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
            # Ensure GT format exists
            if "GT" not in header.formats:
                header.formats.add(
                    "GT", number=1, type="String", description="Genotype"
                )
            # Add contigs from reference FASTA to header
            fasta_ref = pysam.FastaFile(str(reference))
            for chrom, length in zip(fasta_ref.references, fasta_ref.lengths):
                if chrom not in header.contigs:
                    header.contigs.add(chrom, length=length)
            fasta_ref.close()

    # Add samples to header
    for _, sample in opened_vcfs:
        if header and sample not in header.samples:
            header.add_sample(sample)

    # Merge records
    for vcf, sample in opened_vcfs:
        try:
            # TODO: Streaming records would avoid constructing a potentially large
            # in-memory dictionary when merging many files.
            for rec in vcf.fetch():
                key = (rec.chrom, rec.pos, rec.ref)
                if key not in records:
                    # First time seeing this record
                    alleles = (
                        (rec.ref,) + tuple(rec.alts) if rec.alts else (rec.ref, ".")
                    )
                    # Create new record in header
                    new_rec = header.new_record(
                        contig=rec.chrom,
                        start=rec.pos - 1,
                        stop=rec.stop,
                        id=rec.id,
                        alleles=alleles,
                    )
                    records[key] = new_rec
                else:
                    # Merge alternate alleles
                    new_rec = records[key]
                    existing = {a for a in new_rec.alts or [] if a != "."}
                    extra = {a for a in (rec.alts or []) if a != "."}
                    merged = sorted(existing.union(extra))
                    new_rec.alts = merged if merged else ["."]
                # Set genotype for this sample
                gt = rec.samples[0].get("GT")
                records[key].samples[sample]["GT"] = gt
        except Exception as e:
            raise ValueError(f"error: {e}")

    # Write merged VCF
    with pysam.VariantFile(str(out_vcf), "w", header=header) as out:
        for rec in sorted(records.values(), key=lambda r: (r.contig, r.pos)):
            out.write(rec)


def main(argv: list[str] | None = None) -> int:
    args = _parse_args(argv or sys.argv[1:])
    inputs: list[tuple[Path, str]] = []
    for spec in args.inputs:
        if ":" in spec:
            path_str, sample = spec.split(":", 1)
        else:
            path_str = spec
            sample = Path(path_str).stem
        inputs.append((Path(path_str), sample))
    _merge_records(inputs, Path(args.output), Path(args.reference))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
