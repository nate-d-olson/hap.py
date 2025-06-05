"""Simplified hapenum tool placeholder.

This module sketches a Python replacement for the original C++ hapenum
utility. It parses a VCF and emits a rudimentary DOT representation of
possible haplotypes for a given region. The implementation is intentionally
minimal and only suitable for small test VCFs.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pysam


def _parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser("hapenum")
    parser.add_argument("vcf", help="Input VCF file")
    parser.add_argument("-r", "--reference", required=True, help="Reference FASTA")
    parser.add_argument("-l", "--location", required=True, help="Genomic location")
    parser.add_argument("--output-dot", required=True, help="Output DOT file")
    return parser.parse_args(argv)


def _enumerate(vcf: pysam.VariantFile, region: str, out_dot: Path) -> None:
    chrom, *rest = region.replace("chr", "").split(":")
    start = 1
    if rest:
        start = int(rest[0].split("-")[0])
    nodes = []
    edges = []
    idx = 0
    nodes.append(
        f'node_{idx}  [label="start: {chrom}:{start}", shape=diamond, penwidth=2, color=cadetblue4] ;'
    )
    prev = [idx]
    for rec in vcf.fetch(region=region):
        idx += 1
        nodes.append(
            f"node_{idx}  [label=\"alt: {','.join(rec.alts)} {rec.chrom}:{rec.pos}-{rec.stop-1}\", shape=parallelogram, penwidth=2, color=cadetblue4] ;"
        )
        for p in prev:
            edges.append(f"node_{p} -> node_{idx} [] ;")
        prev = [idx]
    idx += 1
    nodes.append(
        f"node_{idx}  [label=\"end: {chrom}:{rec.stop-1 if 'rec' in locals() else start}\", shape=triangle, penwidth=2, color=cadetblue4] ;"
    )
    for p in prev:
        edges.append(f"node_{p} -> node_{idx} [] ;")
    with open(out_dot, "w", encoding="utf-8") as out:
        out.write("digraph G {\n")
        for n in nodes:
            out.write(n + "\n")
        out.write("\n")
        for e in edges:
            out.write(e + "\n")
        out.write("}\n")


def main(argv: list[str] | None = None) -> int:
    args = _parse_args(argv or sys.argv[1:])
    vcf = pysam.VariantFile(args.vcf)
    _enumerate(vcf, args.location, Path(args.output_dot))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
