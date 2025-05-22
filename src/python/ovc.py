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
# 10/1/2015
#
# Testing helper, will return exit code 1 if the passed bed file has
# overlapping intervals.
#
# Usage:
#
# ovc.py input.bed
#
# Author:
#
# Peter Krusche <pkrusche@illumina.com>
#

import argparse
import logging
import sys
import traceback


def main() -> int:
    """
    Check a BED file for overlapping intervals.

    Returns:
        int: 0 if no overlaps, 1 if overlaps found
    """
    parser = argparse.ArgumentParser("Overlap Checker")
    parser.add_argument("bed_file", help="BED file to check for overlaps")
    parser.add_argument(
        "-v", "--verbose", action="store_true", help="Show detailed error messages"
    )

    args = parser.parse_args()

    logging_level = logging.INFO if args.verbose else logging.WARNING
    logging.basicConfig(level=logging_level, format="%(levelname)s: %(message)s")

    try:
        with open(args.bed_file, encoding="utf-8") as f:
            last = -1
            lines = 1

            for line in f:
                chunk = line.split("\t")
                if len(chunk) > 3 and (last - 1) > int(chunk[1]):
                    print(f"Overlap at {chunk[0]}:{int(chunk[1])} (line {lines})")
                    return 1
                elif len(chunk) > 3:
                    last = int(chunk[2])
                lines += 1

        return 0
    except Exception as e:
        logging.error(f"Error: {str(e)}")
        if args.verbose:
            traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())
