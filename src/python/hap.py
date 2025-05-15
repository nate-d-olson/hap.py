#!/usr/bin/env python3
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
# 9/9/2014
#
# Diploid VCF File Comparison
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
import subprocess
import sys
import tempfile
import traceback
from typing import Any, Dict, List, Optional, Set, Tuple, Union

scriptDir = os.path.abspath(os.path.dirname(__file__))
sys.path.append(os.path.abspath(os.path.join(scriptDir, "..", "lib", "python3")))

import Tools
from Tools import vcfextract


def parseArgs():
    """Parse command line arguments"""

    parser = argparse.ArgumentParser("Haplotype Comparison")

    parser.add_argument(
        "ref", help="Reference sequence file (FASTA)", default=Tools.defaultReference()
    )

    parser.add_argument("query", help="Query VCF file")

    parser.add_argument("truth", help="Truth VCF file")

    parser.add_argument(
        "-o",
        "--reports-prefix",
        dest="prefix",
        help="Prefix for report output files",
        default=None,
    )

    parser.add_argument(
        "-V",
        "--verbose",
        dest="verbose",
        help="Increase verbosity of output",
        action="count",
        default=0,
    )

    parser.add_argument(
        "--logfile",
        dest="logfile",
        help="Write logging information into file rather than to stderr",
        default=None,
    )

    parser.add_argument(
        "--no-write-counts",
        dest="writeCounts",
        help="Do not write count metrics for feature stratification",
        action="store_false",
        default=True,
    )

    parser.add_argument(
        "-r",
        "--regions",
        dest="regions",
        help="Regions to restrict analysis to, must have corresponding .bed file. "
        ", separated. e.g. chr1,chr2.",
        default=None,
    )

    parser.add_argument(
        "-R",
        "--regions-file",
        dest="regions_file",
        help="Restrict analysis to given (sparse) regions. The file must be a BED file.",
        default=None,
    )

    parser.add_argument(
        "-t",
        "--type",
        dest="type",
        help="Variant types and regions to include, e.g. INDEL, SNP, ALL, NOCALL. "
        "Specify multiple types with , or use ALL to include everything",
        default="ALL",
    )

    parser.add_argument(
        "-f",
        "--feature-table",
        nargs="+",
        dest="feature_table",
        help="Features to annotate variants with and stratify on.",
        default=[],
    )

    parser.add_argument(
        "-m",
        "--stratification",
        dest="stratification",
        help="List of stratifications, comma-separated (super- and subfeatures "
        "are joined by :). For example: CONF,CONF:HET.",
    )

    parser.add_argument(
        "--confident-regions-truth",
        dest="conf_truth",
        help="Confident call regions in truth (e.g. from GIAB), a BED file.",
    )

    parser.add_argument(
        "--confident-regions-query",
        dest="conf_query",
        help="Confident call regions in query (e.g. from GIAB), a BED file.",
    )

    parser.add_argument(
        "--no-optimize",
        dest="optimize",
        help="Do not use DP optimization, ensures optimal Answer at the expense "
        "of runtime for complex regions.",
        action="store_false",
        default=True,
    )

    parser.add_argument(
        "-X",
        "--roc",
        dest="roc",
        help="Enable ROC computation (use QUAL for roc-feature to take the VCF QUAL)",
        default=None,
        action="store_true",
    )

    parser.add_argument(
        "--roc-filter",
        dest="roc_filter",
        help="When using -X, info field to use, e.g. GQX",
        default=None,
    )

    parser.add_argument(
        "--roc-delta",
        dest="roc_delta",
        help="When using -X, precision delta for ROC (for faster ROC calculation)",
        default=0.1,
        type=float,
    )

    parser.add_argument(
        "-P",
        "--preprocessed",
        dest="preprocessing",
        help=(
            "Don't preprocess, use these files (comma-separated: "
            "truth.vcf.gz,query.vcf.gz[,regions.bed.gz])"
        ),
        default=None,
    )

    parser.add_argument(
        "--usefiltered-truth",
        dest="usefiltered_truth",
        help="Use filtered variant calls in truth file (by default, only PASS is used)",
        default=False,
        action="store_true",
    )

    parser.add_argument(
        "--usefiltered-query",
        dest="usefiltered_query",
        help="Use filtered variant calls in query file (by default, only PASS is used)",
        default=False,
        action="store_true",
    )

    parser.add_argument(
        "--no-leftshift",
        dest="leftshift",
        help="Don't left-shift variants before comparison (default is to left-shift)",
        action="store_false",
        default=True,
    )

    parser.add_argument(
        "--location-features",
        dest="location_features",
        help=(
            "Create a 3-column feature table with a feature name / contig / region; "
            "this can be used to create stratifications for each contig / region combination"
        ),
        default=None,
    )

    parser.add_argument(
        "--no-decompose",
        dest="decompose",
        help="Don't decompose complex variants (i.e. MNPs, \
                           complex substitutions with length > 1)",
        action="store_false",
        default=True,
    )

    parser.add_argument(
        "--bcftools-norm",
        dest="bcftools_norm",
        help="Use bcftools norm instead of our internal implementation",
        action="store_true",
        default=False,
    )

    parser.add_argument(
        "--python-only",
        dest="pythononly",
        help="Use hap.py Python implementation (not the quantify executable)",
        action="store_true",
        default=False,
    )

    parser.add_argument(
        "--adjust-conf-regions",
        dest="adjust_conf_regions",
        help="Automatically adjust confident regions based on distance \
                           from global / filtered variants",
        action="store_true",
        default=False,
    )

    parser.add_argument(
        "--unhappy",
        dest="unhappy",
        help="Use legacy bcftools-based comparison for platform testing purposes",
        action="store_true",
        default=False,
    )

    parser.add_argument(
        "--lose",
        dest="lose",
        help="Set this to enable partial credit for heterozygous calls (diploid test)",
        default=0,
        type=float,
    )

    parser.add_argument(
        "--engine",
        dest="engine",
        help="Comparison engine to use.",
        default="xcmp",
        choices=["xcmp", "vcfeval", "scmp-somatic", "scmp-distance"],
    )

    parser.add_argument(
        "--engine-vcfeval-path",
        dest="vcfeval_path",
        help="Path to RTG Tools vcfeval executable",
        default=None,
    )

    parser.add_argument(
        "--engine-vcfeval-template",
        dest="vcfeval_template",
        help="Path to RTG Tools SDF template",
        default=None,
    )

    parser.add_argument(
        "--preserveAllVariants",
        dest="preserve_all_variants",
        help="Use --preserve-all-variants VCFeval setting (take the first truth match "
        "for each query, using the match number as allele count.)",
        action="store_true",
        default=False,
    )

    parser.add_argument(
        "--threads",
        dest="threads",
        help="Number of threads to use",
        default=multiprocessing.cpu_count(),
    )

    parser.add_argument(
        "--scratch-prefix",
        dest="scratch_prefix",
        help="Directory for scratch files",
        default=None,
    )

    parser.add_argument(
        "--output-vtc",
        dest="output_vtc",
        help="Write variant compare outputs in VTC format if prefix is specified",
        action="store_true",
        default=False,
    )

    parser.add_argument(
        "--gender",
        dest="gender",
        help="Gender, determines how to handle X/Y",
        default="female",
        choices=["male", "female", "auto", "none"],
    )

    parser.add_argument(
        "--window-size",
        dest="window",
        help="Preprocessing window size. Set this to 0 to disable windowing.",
        default=20000,
        type=int,
    )

    parser.add_argument(
        "--xcmp-enumeration-threshold",
        dest="enum_threshold",
        help=(
            "Maximum number of enumerations per window for xcmp (use > 0 to disable "
            "window-based decomposition of difficult regions)."
        ),
        default=10,
        type=int,
    )

    parser.add_argument(
        "--pass-only",
        dest="pass_only",
        help="Keep only PASS variants",
        action="store_true",
        default=False,
    )

    parser.add_argument(
        "--fixchr",
        dest="fixchr",
        help="Add chr prefix to VCF records where needed",
        action="store_true",
        default=False,
    )

    parser.add_argument(
        "--no-fixchr",
        dest="fixchr",
        help="Do not fix chr prefix",
        action="store_false",
        default=False,
    )

    parser.add_argument(
        "--set-globals",
        dest="globals",
        help="Set global flagss for preprocessing",
        default=None,
    )

    parser.add_argument(
        "--force-interactive",
        dest="forceinteractive",
        help="Force to run interactively, i.e. do not submit to SGE",
        action="store_true",
        default=False,
    )

    parser.add_argument(
        "--force-overwrite",
        dest="forceoverwrite",
        help="Force to overwrite existing result files.",
        action="store_true",
        default=False,
    )

    parser.add_argument(
        "--preprocess-truth",
        dest="preprocess_truth",
        help=(
            "Preprocessing flag string for truth, e.g. DECOMPOSE,\
                            NO_COMPLEX_BRANCH,NO_LEFTSHIFT"
        ),
        default=None,
    )

    parser.add_argument(
        "--preprocess-query",
        dest="preprocess_query",
        help=(
            "Preprocessing flag string for query, e.g. DECOMPOSE,\
                            NO_COMPLEX_BRANCH,NO_LEFTSHIFT"
        ),
        default=None,
    )

    parser.add_argument(
        "--preprocess-window",
        dest="preprocess_window",
        help="Window size for preprocessing",
        type=int,
        default=1000,
    )

    parser.add_argument(
        "--write-vcf",
        dest="write_vcf",
        help="Write annotated VCF file for debugging",
        action="store_true",
        default=False,
    )

    parser.add_argument(
        "--output-vtc-max-size",
        dest="output_vtc_max_size",
        help="Maximum size for output VTC file",
        default=2048,
        type=int,
    )

    parser.add_argument(
        "--write-counts",
        dest="write_counts",
        help="Write count metrics for feature stratification",
        action="store_true",
    )

    parser.add_argument(
        "--false-positives",
        dest="false_positives",
        help="Suppress false positives from output, taking the given BED file as truth",
        default=None,
    )

    parser.add_argument(
        "--qq-plot",
        dest="qq_plot",
        help=(
            "Make Q-Q plots for all stratifications; "
            "this only works when false_positives is given."
        ),
        action="store_true",
        default=False,
    )

    parser.add_argument(
        "--no-fixchr-truth",
        dest="fixchr_truth",
        help="Do not fix chr prefixes in truth",
        action="store_false",
        default=True,
    )

    parser.add_argument(
        "--no-fixchr-query",
        dest="fixchr_query",
        help="Do not fix chr prefixes in query",
        action="store_false",
        default=True,
    )

    result = parser.parse_args()

    # check engine options

    if result.engine == "vcfeval":
        if not result.vcfeval_path:
            parser.error(
                "Path to vcfeval (--engine-vcfeval-path) required when engine=vcfeval"
            )
        if not os.path.exists(result.vcfeval_path):
            parser.error("RTG Tools vcfeval not found at %s" % result.vcfeval_path)

        if not result.vcfeval_template:
            parser.error(
                "Path to vcfeval template (--engine-vcfeval-template) required when engine=vcfeval"
            )
        if not os.path.exists(result.vcfeval_template):
            parser.error(
                "RTG Tools vcfeval template not found at %s" % result.vcfeval_template
            )

    elif result.engine != "xcmp" and not result.unhappy:
        parser.error(
            "Engine %s is not supported in this version of hap.py" % result.engine
        )

    if result.unhappy and result.engine != "xcmp":
        parser.error("--unhappy can only be used with the xcmp engine")

    if result.threads:
        result.threads = int(result.threads)

    if result.type:
        result.type = result.type.upper()

    # if no regions were given, this will look at all of them.
    if result.regions:
        # regions are specified as VCF chrom IDs
        result.regions = result.regions.split(",")

    return result


def main():
    """Main method"""
    try:
        args = parseArgs()

        # Default threshold for ROC
        flt_threshold = -1

        # how verbose?
        if args.verbose == 0:
            loglevel = logging.WARN
        elif args.verbose == 1:
            loglevel = logging.INFO
        else:
            loglevel = logging.DEBUG

        if args.logfile:
            logging.basicConfig(filename=args.logfile, level=loglevel)
        else:
            logging.basicConfig(level=loglevel)

        logging.getLogger("requests").setLevel(logging.WARNING)

        if (
            args.prefix is not None
            and os.path.exists(args.prefix + ".summary.csv")
            and not args.forceoverwrite
        ):
            # check if the output file exists already
            logging.error(
                "Results file %s exists already. Use --force-overwrite to overwrite."
                % (args.prefix + ".summary.csv")
            )
            return 1

        if args.engine == "xcmp":
            # Use hap.py
            logging.info("Using haplotype comparison engine.")
            import Haplo.haplotypes
            import Haplo.happyroc
            import Haplo.quantify

            Haplo.quantify.run(args)

            if args.roc:
                if not args.roc_filter:
                    logging.info("Creating ROC curve using QUAL")
                    filter_key = "QUAL"
                else:
                    logging.info("Creating ROC curve using %s" % args.roc_filter)
                    filter_key = args.roc_filter
                try:
                    Haplo.happyroc.makeRocs(
                        args.prefix, filter_key, args.roc_delta, args.type
                    )
                except Exception as e:
                    logging.error(str(e))
                    logging.info(
                        "Could not create ROC file, possibly because there were not enough TP/FP variants."
                    )
        elif args.engine == "vcfeval":
            # Use RTG Tools vcfeval
            logging.info("Using RTG vcfeval engine.")

            # Long import...
            from Tools import rtgtools

            rtgtools.runRTGTools(args)

            # Make a happy file
            vcfextract.makeHappyVCF(args)

            if args.roc:
                if not args.roc_filter:
                    logging.info("Creating ROC curve using QUAL")
                    filter_key = "QUAL"
                else:
                    logging.info("Creating ROC curve using %s" % args.roc_filter)
                    filter_key = args.roc_filter
                try:
                    import Haplo.happyroc

                    Haplo.happyroc.makeRocs(
                        args.prefix, filter_key, args.roc_delta, args.type
                    )
                except Exception as e:
                    logging.error(str(e))
                    logging.info(
                        "Could not create ROC file, possibly because there were not enough TP/FP variants."
                    )
        else:
            logging.error("Unsupported engine: %s" % args.engine)

        return 0
    except Exception as e:
        print(str(e))
        traceback.print_exc(file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
