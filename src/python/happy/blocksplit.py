#!/usr/bin/env python3
"""CLI wrapper for blocksplit using Python implementation."""
import argparse
import shutil
import sys

from Haplo.blocksplit import blocksplitWrapper


def main():
    parser = argparse.ArgumentParser(
        prog="blocksplit", description="Split variant blocks into chunks"
    )
    parser.add_argument("vcf1", help="First VCF file")
    parser.add_argument("vcf2", help="Second VCF file")
    parser.add_argument("-o", dest="output", required=True, help="Output BED file path")
    parser.add_argument(
        "-l", dest="location", default=None, help="Chromosome or region to restrict"
    )
    parser.add_argument(
        "-w", dest="window", type=int, default=10000, help="Window size (bp)"
    )
    parser.add_argument(
        "--nblocks", dest="pieces", type=int, default=1, help="Number of blocks"
    )
    parser.add_argument(
        "--scratch",
        dest="scratch_prefix",
        default=None,
        help="Scratch directory for temp files",
    )
    args = parser.parse_args()

    class ArgsObj:
        pass

    obj = ArgsObj()
    obj.vcf1 = args.vcf1
    obj.vcf2 = args.vcf2
    obj.window = args.window
    obj.pieces = args.pieces
    obj.scratch_prefix = args.scratch_prefix

    try:
        bedfile, _ = blocksplitWrapper(args.location, obj)
        shutil.move(bedfile, args.output)
    except Exception as e:
        sys.exit(f"blocksplit error: {e}")


if __name__ == "__main__":  # pragma: no cover
    main()
