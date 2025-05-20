#!/usr/bin/env python3
#
# Copyright (c) 2010-2015 Illumina, Inc.
# All rights reserved.
#
# This file is distributed under the simplified BSD license.
# The full text can be found here (and in LICENSE.txt in the root folder of
# this distribution):
#
# https://github.com/Illumina/licenses/blob/master/Simplified-BSD-License.txt
#
# 9/9/2014
#
# Extract caller and aligner versions from VCF and BAM headers
#
# Usage:
#
# For usage instructions run with option --help
#
# Author:
#
# Peter Krusche <pkrusche@illumina.com>
#

import argparse
import json
import logging
import os
import sys
import traceback

scriptDir = os.path.abspath(os.path.dirname(os.path.realpath(__file__)))
# Update path for Python 3
lib_path = os.path.abspath(os.path.join(scriptDir, "..", "lib", "python3"))
if os.path.exists(lib_path):
    sys.path.append(lib_path)
else:
    fallback_path = os.path.abspath(os.path.join(scriptDir, "..", "lib"))
    sys.path.append(fallback_path)

import Tools
from Tools.vcfcallerinfo import CallerInfo


# noinspection PyBroadException
def main() -> int:
    """
    Extract caller and aligner information from VCF and BAM headers.

    Returns:
        int: 0 on success, 1 on failure
    """
    parser = argparse.ArgumentParser("Extract caller / aligner info from headers")

    parser.add_argument("input", help="Input VCF file")

    parser.add_argument(
        "-o",
        "--output",
        dest="output",
        required=True,
        help="Output file name (json format)",
    )

    parser.add_argument(
        "-b",
        "--bam",
        dest="bam",
        default=None,
        help="pass a BAM file for aligner name/version extraction",
    )

    args = parser.parse_args()

    ci = CallerInfo()
    ci.addVCF(args.input)
    if args.bam:
        ci.addBAM(args.bam)

    if not args.output.endswith(".json"):
        args.output += ".json"

    logging.info(f"Writing {args.output}")
    with open(args.output, "w", encoding="utf-8") as fp:
        json.dump(ci.asDict(), fp, sort_keys=True, indent=4, separators=(",", ": "))

    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as e:
        logging.error(e.decode("utf-8") if isinstance(e, bytes) else str(e))
        traceback.print_exc(file=Tools.LoggingWriter(logging.ERROR))
        sys.exit(1)
