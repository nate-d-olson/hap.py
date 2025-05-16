# coding=utf-8
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
# 3/9/2014
#
# Author:
#
# Peter Krusche <pkrusche@illumina.com>
#

import logging
import os
import tempfile
import traceback

from Tools import LoggingWriter
from Tools.bcftools import runBcftools


def runSCmp(vcf1, vcf2, target, args):
    """Runs scmp, which outputs a file quantify can produce counts on
    vcf1 and vcf2 must be indexed and only contain a single sample column.
    """

    try:
        if args.engine == "scmp-distance":
            cmode = "distance"
        else:
            cmode = "alleles"

        tf = tempfile.NamedTemporaryFile(delete=False)
        tf.close()
        try:
            # change GTs so we can compare them
            vargs = ["merge", "--force-samples", vcf1, vcf2, "-o", tf.name]
            runBcftools(*vargs)
            vargs = [
                "view",
                tf.name,
                "|",
                "scmp",
                "-M",
                cmode,
                "-",
                "-r",
                args.ref,
                "--threads",
                # Ensure str() conversion is safe for Python 3
                (
                    str(args.threads)
                    if not isinstance(args.threads, bytes)
                    else args.threads.decode("utf-8")
                ),
                "-o",
                target,
            ]
            if args.roc:
                vargs += ["--q", args.roc]

            vargs += [
                "--distance-maxdist",
                (
                    str(args.engine_scmp_distance)
                    if not isinstance(args.engine_scmp_distance, bytes)
                    else args.engine_scmp_distance.decode("utf-8")
                ),
            ]
            runBcftools(*vargs)
        finally:
            os.remove(tf.name)

        if target.endswith(".vcf.gz"):
            runBcftools("index", "-t", target)
            return [target, target + ".tbi"]
        else:
            runBcftools("index", target)
            return [target, target + ".csi"]
    except Exception as e:
        # Handle case where error message might be bytes in Python 3
        error_msg = e.decode("utf-8") if isinstance(e, bytes) else e.decode('utf-8') if isinstance(e, bytes) else str(e)
        logging.error("Exception when running scmp: %s" % error_msg)
        logging.error("-" * 60)
        traceback.print_exc(file=LoggingWriter(logging.ERROR))
        logging.error("-" * 60)
        raise
    except BaseException as e:
        # Handle case where error message might be bytes in Python 3
        error_msg = e.decode("utf-8") if isinstance(e, bytes) else e.decode('utf-8') if isinstance(e, bytes) else str(e)
        logging.error("Exception when running scmp: %s" % error_msg)
        logging.error("-" * 60)
        traceback.print_exc(file=LoggingWriter(logging.ERROR))
        logging.error("-" * 60)
        raise
