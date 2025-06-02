#!/usr/bin/env python3
"""
hap: CLI entry point for hap.py benchmarking tool.
"""
import argparse
import logging
import sys
import traceback

try:
    from happy import qfy
except ImportError:
    qfy = None


def main():
    parser = argparse.ArgumentParser(prog="hap.py", description="Haplotype Comparison")
    # Show version
    parser.add_argument(
        "-v",
        "--version",
        action="store_true",
        dest="version",
        help="Show version and exit",
    )
    parser.add_argument("-r", "--reference", dest="ref", help="Reference FASTA file")
    parser.add_argument(
        "-o", "--report-prefix", dest="reports_prefix", help="Output prefix"
    )
    parser.add_argument("--scratch-prefix", dest="scratch_prefix", help="Scratch dir")
    parser.add_argument(
        "--keep-scratch",
        dest="keep_scratch",
        action="store_true",
        help="Do not delete scratch files",
    )
    # Accept legacy CLI flags for chromosome limiting and modes
    parser.add_argument(
        "-l", "--chrom", dest="chrom", nargs="+", help=argparse.SUPPRESS
    )
    parser.add_argument(
        "--force-interactive",
        dest="force_interactive",
        action="store_true",
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "--unhappy", dest="unhappy", action="store_true", help=argparse.SUPPRESS
    )
    parser.add_argument(
        "--pass-only", dest="pass_only", action="store_true", help=argparse.SUPPRESS
    )
    # Ensure basic quantification flags are recognized
    parser.add_argument(
        "-V",
        "--write-vcf",
        dest="write_vcf",
        action="store_true",
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "-X",
        "--write-counts",
        dest="write_counts",
        action="store_true",
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "--output-vtc", dest="output_vtc", action="store_true", help=argparse.SUPPRESS
    )
    # Include full quantification args if available
    if qfy and hasattr(qfy, "updateArgs"):
        qfy.updateArgs(parser)
    args = parser.parse_args()
    try:
        # Run the quantification subcommand
        qfy.quantify(args)
    except Exception as e:
        logging.error(str(e))
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
