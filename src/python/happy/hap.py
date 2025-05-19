#!/usr/bin/env python3
"""
hap: CLI entry point for hap.py benchmarking tool.
"""
import argparse
import logging
import os
import sys
import tempfile
import traceback

from happy import pre, qfy
from happy.Haplo import blocksplit, gvcf2bed, partialcredit, quantify, vcfeval
from happy.Tools import bcftools, vcfextract
from happy.Tools.bcftools import bedOverlapCheck
from happy.Tools.fastasize import fastaContigLengths
from happy.Tools.parallel import getPool, runParallel
from happy.Tools.sessioninfo import sessionInfo


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
    # Include quantification args
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
