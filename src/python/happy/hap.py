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
    # Accept legacy CLI flag for chromosome limiting (single value)
    parser.add_argument("-l", "--chrom", dest="chrom", help=argparse.SUPPRESS)
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
    # Include full quantification args (write-vcf, write-counts, output-vtc, etc.)
    if qfy and hasattr(qfy, "updateArgs"):
        qfy.updateArgs(parser)
    # Parse options first; capture two positional VCF inputs: truth and query
    args, unknown = parser.parse_known_args()
    if len(unknown) < 2:
        parser.error("the following arguments are required: truth_vcf, query_vcf")
    args.truth_vcf, args.query_vcf = unknown[0], unknown[1]

    # Run quantification or comparison backend

    # Default: run quantification
    try:
        qfy.quantify(args)
    except Exception as e:
        logging.error(str(e))
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
