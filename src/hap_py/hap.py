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
import multiprocessing
import os
import sys
import tempfile
import time
import traceback
from pathlib import Path

import contextlib

# Modern imports using the new package structure
from . import pre, qfy
from .haplo import gvcf2bed, partialcredit, quantify, vcfeval
from .tools import bcftools, vcfextract
from .tools.bcftools import bedOverlapCheck
from .tools.fastasize import fastaContigLengths
from .tools.parallel import getPool
from .tools.sessioninfo import sessionInfo
from .tools import version


def main() -> int:
    parser = argparse.ArgumentParser("Haplotype Comparison")

    # input
    parser.add_argument(
        "-v",
        "--version",
        dest="version",
        action="store_true",
        help="Show version number and exit.",
    )

    parser.add_argument(
        "-r", "--reference", dest="ref", default=None, help="Specify a reference file."
    )

    # output
    parser.add_argument(
        "-o",
        "--report-prefix",
        dest="reports_prefix",
        default=None,
        help="Filename prefix for report output.",
    )
    parser.add_argument(
        "--scratch-prefix",
        dest="scratch_prefix",
        default=None,
        help="Directory for scratch files.",
    )
    parser.add_argument(
        "--keep-scratch",
        dest="delete_scratch",
        default=True,
        action="store_false",
        help="Filename prefix for scratch report output.",
    )

    # add quantification args
    qfy.updateArgs(parser)

    # control preprocessing
    pre.updateArgs(parser)
    parser.add_argument(
        "--convert-gvcf-truth",
        dest="convert_gvcf_truth",
        action="store_true",
        default=False,
        help="Convert the truth set from genome VCF format to a VCF before processing.",
    )
    parser.add_argument(
        "--convert-gvcf-query",
        dest="convert_gvcf_query",
        action="store_true",
        default=False,
        help="Convert the query set from genome VCF format to a VCF before processing.",
    )
    parser.add_argument(
        "--preprocess-truth",
        dest="preprocessing_truth",
        action="store_true",
        default=False,
        help="Preprocess truth file with same settings as query (default is to accept truth in original format).",
    )
    parser.add_argument(
        "--usefiltered-truth",
        dest="usefiltered_truth",
        action="store_true",
        default=False,
        help="Use filtered variant calls in truth file (by default, only PASS calls in the truth file are used)",
    )
    parser.add_argument(
        "--preprocessing-window-size",
        dest="preprocess_window",
        default=10000,
        type=int,
        help="Preprocessing window size (variants further apart than that size are not expected to interfere).",
    )
    parser.add_argument(
        "--adjust-conf-regions",
        dest="preprocessing_truth_confregions",
        action="store_true",
        default=True,
        help="Adjust confident regions to include variant locations. Note this will only include variants "
        "that are included in the CONF regions already when viewing with bcftools; this option only "
        "makes sure insertions are padded correctly in the CONF regions (to capture these, both the "
        "base before and after must be contained in the bed file).",
    )
    parser.add_argument(
        "--no-adjust-conf-regions",
        dest="preprocessing_truth_confregions",
        action="store_false",
        help="Do not adjust confident regions for insertions.",
    )

    # detailed control of comparison
    parser.add_argument(
        "--unhappy",
        "--no-haplotype-comparison",
        dest="no_hc",
        action="store_true",
        default=False,
        help="Disable haplotype comparison (only count direct GT matches as TP).",
    )

    parser.add_argument(
        "-w",
        "--window-size",
        dest="window",
        default=50,
        type=int,
        help="Minimum distance between variants such that they fall into the same superlocus.",
    )
    parser.add_argument(
        "--threads",
        dest="threads",
        default=multiprocessing.cpu_count(),
        type=int,
        help="Number of threads to use.",
    )

    parser.add_argument(
        "--engine",
        dest="engine",
        default="vcfeval",
        choices=["vcfeval"],
        help="Comparison engine to use. Only vcfeval is supported in Python 3.",
    )

    parser.add_argument(
        "--engine-vcfeval-path",
        dest="engine_vcfeval",
        required=False,
        default=vcfeval.findVCFEval(),
        help='This parameter should give the path to the "rtg" executable. '
        "The default is %s" % vcfeval.findVCFEval(),
    )

    parser.add_argument(
        "--engine-vcfeval-template",
        dest="engine_vcfeval_template",
        required=False,
        help="Vcfeval needs the reference sequence formatted in its own file format "
        "(SDF -- run rtg format -o ref.SDF ref.fa). You can specify this here "
        "to save time when running hap.py with vcfeval. If no SDF folder is "
        "specified, hap.py will create a temporary one.",
    )

    # Remove SGE dependency - modern systems don't typically use SGE
    has_sge = False
    if has_sge:
        parser.add_argument(
            "--force-interactive",
            dest="force_interactive",
            default=False,
            action="store_true",
            help="Force running interactively (i.e. when JOB_ID is not in the environment)",
        )

    parser.add_argument("_vcfs", help="Two VCF files.", default=[], nargs="*")

    parser.add_argument(
        "--logfile",
        dest="logfile",
        default=None,
        help="Write logging information into file rather than to stderr",
    )

    verbosity_options = parser.add_mutually_exclusive_group(required=False)

    verbosity_options.add_argument(
        "--verbose",
        dest="verbose",
        default=False,
        action="store_true",
        help="Raise logging level from warning to info.",
    )

    verbosity_options.add_argument(
        "--quiet",
        dest="quiet",
        default=False,
        action="store_true",
        help="Set logging level to output errors only.",
    )

    args, unknown_args = parser.parse_known_args()

    if not has_sge:
        args.force_interactive = True

    if args.verbose:
        loglevel = logging.INFO
    elif args.quiet:
        loglevel = logging.ERROR
    else:
        loglevel = logging.WARNING

    # reinitialize logging
    for handler in logging.root.handlers[:]:
        logging.root.removeHandler(handler)
    logging.basicConfig(
        filename=args.logfile,
        format="%(asctime)s %(levelname)-8s %(message)s",
        level=loglevel,
    )

    # remove some safe unknown args
    unknown_args = [x for x in unknown_args if x not in ["--force-interactive"]]
    if len(sys.argv) < 2 or len(unknown_args) > 0:
        if unknown_args:
            logging.error(f"Unknown arguments specified: {unknown_args}")
        parser.print_help()
        exit(1)

    print(f"Hap.py {version.version}")
    if args.version:
        exit(0)

    if args.roc:
        args.write_vcf = True

    # sanity-check regions bed file (HAP-57)
    if args.regions_bedfile:
        logging.info("Checking input regions.")
        if bedOverlapCheck(args.regions_bedfile):
            raise ValueError(
                "The regions bed file (specified using -R) has overlaps. "
                "You can either use -T, or run the file through bedtools merge"
            )

    if args.fp_bedfile and not os.path.exists(args.fp_bedfile):
        raise FileNotFoundError("FP/confident call region bed file does not exist.")

    if not args.force_interactive and "JOB_ID" not in os.environ:
        parser.print_help()
        raise RuntimeError(
            "Please qsub me so I get approximately 1 GB of RAM per thread."
        )

    if not args.ref:
        args.ref = None  # Remove default reference dependency

    if not args.ref or not os.path.exists(args.ref):
        raise FileNotFoundError("Please specify a valid reference path using -r.")

    if not args.reports_prefix:
        raise ValueError("Please specify an output prefix using -o")

    if not os.path.exists(os.path.dirname(os.path.abspath(args.reports_prefix))):
        raise FileNotFoundError(
            "The output path does not exist. Please specify a valid output path and prefix using -o"
        )

    if os.path.basename(args.reports_prefix) == "" or os.path.isdir(
        args.reports_prefix
    ):
        raise ValueError(
            "The output path should specify a file name prefix. Please specify a valid output path "
            "and prefix using -o. For example, -o /tmp/test will create files named /tmp/test* ."
        )

    # noinspection PyProtectedMember
    if not args._vcfs or len(args._vcfs) != 2:
        raise ValueError("Please specify exactly two input VCFs.")

    # noinspection PyProtectedMember
    args.vcf1 = args._vcfs[0]
    # noinspection PyProtectedMember
    args.vcf2 = args._vcfs[1]

    if not os.path.exists(args.vcf1):
        raise FileNotFoundError(f"Input file {args.vcf1} does not exist.")
    if not os.path.exists(args.vcf2):
        raise FileNotFoundError(f"Input file {args.vcf2} does not exist.")

    tempfiles = []

    # Define internal format suffix
    internal_format_suffix = ".vcf.gz"

    # write session info and args file
    session = sessionInfo()
    session["final_args"] = args.__dict__
    with open(
        args.reports_prefix + ".runinfo.json", "w", encoding="utf-8"
    ) as sessionfile:
        json.dump(session, sessionfile)

    try:
        logging.info(f"Comparing {args.vcf1} and {args.vcf2}")

        logging.info("Preprocessing truth: %s" % args.vcf1)
        starttime = time.time()

        ttf = tempfile.NamedTemporaryFile(
            delete=False,
            dir=args.scratch_prefix,
            prefix="truth.pp",
            suffix=internal_format_suffix,
        )
        ttf.close()

        if args.preprocessing_truth and args.filter_nonref:
            logging.info("Filtering out any variants genotyped as <NON_REF>")

        ## Only converting truth gvcf to vcf if both arguments are true
        convert_gvcf_truth = False
        if args.convert_gvcf_truth or args.convert_gvcf_to_vcf:
            logging.info("Converting genome VCF to VCF")
            convert_gvcf_truth = True

        tempfiles.append(ttf.name)
        tempfiles.append(ttf.name + ".csi")
        tempfiles.append(ttf.name + ".tbi")
        args.gender = pre.preprocess(
            args.vcf1,
            ttf.name,
            args.ref,
            args.locations,
            None if args.usefiltered_truth else "*",  # filters
            args.fixchr,
            args.regions_bedfile,
            args.targets_bedfile,
            args.preprocessing_leftshift if args.preprocessing_truth else False,
            args.preprocessing_decompose if args.preprocessing_truth else False,
            args.preprocessing_norm if args.preprocessing_truth else False,
            args.preprocess_window,
            args.threads,
            args.gender,
            False,
            "TRUTH",
            filter_nonref=args.filter_nonref if args.preprocessing_truth else False,
            convert_gvcf_to_vcf=convert_gvcf_truth,
        )

        args.vcf1 = ttf.name

        if args.fp_bedfile and args.preprocessing_truth_confregions:
            conf_temp = gvcf2bed.gvcf2bed(
                args.vcf1, args.ref, args.fp_bedfile, args.scratch_prefix
            )
            tempfiles.append(conf_temp)
            args.strat_regions.append(f"CONF_VARS:{conf_temp}")

        h1 = vcfextract.extractHeadersJSON(args.vcf1)

        elapsed = time.time() - starttime
        logging.info(f"preprocess for {args.vcf1} -- time taken {elapsed:.2f}")

        # once we have preprocessed the truth file we can resolve the locations
        # doing this here improves the time for query preprocessing below
        reference_contigs = set(fastaContigLengths(args.ref).keys())

        if not args.locations:
            # default set of locations is the overlap between truth and reference
            args.locations = list(reference_contigs & set(h1["tabix"]["chromosomes"]))
            if not args.locations:
                raise ValueError("Truth and reference have no chromosomes in common!")
        elif type(args.locations) is not list:
            args.locations = args.locations.split(",")

        args.locations = sorted(args.locations)

        logging.info("Preprocessing query: %s" % args.vcf2)
        if args.filter_nonref:
            logging.info("Filtering out any variants genotyped as <NON_REF>")

        ## Only converting truth gvcf to vcf if both arguments are true
        convert_gvcf_query = False
        if args.convert_gvcf_query or args.convert_gvcf_to_vcf:
            logging.info("Converting genome VCF to VCF")
            convert_gvcf_query = True

        starttime = time.time()

        filtering = "*" if args.pass_only else args.filters_only

        qtf = tempfile.NamedTemporaryFile(
            delete=False,
            dir=args.scratch_prefix,
            prefix="query.pp",
            suffix=internal_format_suffix,
        )
        qtf.close()
        tempfiles.append(qtf.name)
        tempfiles.append(qtf.name + ".csi")
        tempfiles.append(qtf.name + ".tbi")

        pre.preprocess(
            args.vcf2,
            qtf.name,
            args.ref,
            str(",".join(args.locations)),
            filtering,
            args.fixchr,
            args.regions_bedfile,
            args.targets_bedfile,
            args.preprocessing_leftshift,
            args.preprocessing_decompose,
            args.preprocessing_norm,
            args.preprocess_window,
            args.threads,
            args.gender,  # same gender as truth above
            False,
            "QUERY",
            filter_nonref=args.filter_nonref,
            convert_gvcf_to_vcf=convert_gvcf_query,
        )

        args.vcf2 = qtf.name
        h2 = vcfextract.extractHeadersJSON(args.vcf2)

        elapsed = time.time() - starttime
        logging.info(f"preprocess for {args.vcf2} -- time taken {elapsed:.2f}")

        if not h1["tabix"]:
            raise RuntimeError("Truth file is not indexed after preprocesing.")

        if not h2["tabix"]:
            raise RuntimeError("Query file is not indexed after preprocessing.")

        for _xc in args.locations:
            if _xc not in h2["tabix"]["chromosomes"]:
                logging.warning(f"No calls for location {_xc} in query!")

        getPool(args.threads)

        # count variants before normalisation
        if "samples" not in h1 or not h1["samples"]:
            raise ValueError("Cannot read sample names from truth VCF file")

        if "samples" not in h2 or not h2["samples"]:
            raise ValueError("Cannot read sample names from query VCF file")

        tf = tempfile.NamedTemporaryFile(
            delete=False,
            dir=args.scratch_prefix,
            prefix="hap.py.result.",
            suffix=internal_format_suffix,
        )
        tf.close()
        tempfiles.append(tf.name)
        tempfiles.append(tf.name + ".tbi")
        tempfiles.append(tf.name + ".csi")
        output_name = tf.name

        if args.engine == "vcfeval":
            tempfiles += vcfeval.runVCFEval(
                args.vcf1, args.vcf2, output_name, args
            )
            # passed to quantify
            args.type = "ga4gh"
        else:
            raise ValueError(f"Unknown comparison engine: {args.engine}")

        if args.preserve_info and args.engine == "vcfeval":
            # if we use vcfeval we need to merge the INFO fields back in.
            # Use context manager for better resource management
            with tempfile.NamedTemporaryFile(suffix=".txt", delete=False) as tf:
                tempfiles.append(tf.name)
                print("TRUTH_IN", file=tf)
                print("QUERY_IN", file=tf)

            # Use context manager for better resource management
            with tempfile.NamedTemporaryFile(
                suffix=".vcf.gz", delete=False
            ) as info_file:
                tempfiles.append(info_file.name)

            bcftools.runBcftools(
                "merge",
                args.vcf1,
                args.vcf2,
                "--force-samples",
                "-m",
                "all",
                "|",
                "bcftools",
                "reheader",
                "-s",
                tf.name,
                "|",
                "bcftools",
                "view",
                "-o",
                info_file.name,
                "-O",
                "z",
            )
            bcftools.runBcftools("index", info_file.name)

            # Use context manager for better resource management
            merged_info_file_path = ""
            with tempfile.NamedTemporaryFile(
                suffix=".vcf.gz", delete=False
            ) as merged_info_file:
                merged_info_file_path = merged_info_file.name
                tempfiles.append(merged_info_file_path)

            bcftools.runBcftools(
                "merge",
                output_name,
                info_file.name,
                "-m",
                "all",
                "|",
                "bcftools",
                "view",
                "-s",
                "^TRUTH_IN,QUERY_IN",
                "-X",
                "-U",
                "-o",
                merged_info_file_path,
                "-O",
                "z",
            )
            output_name = merged_info_file_path

        args.in_vcf = [output_name]
        args.runner = "hap.py"
        qfy.quantify(args)
        return 0

    finally:
        if args.delete_scratch:
            for x in tempfiles:
                with contextlib.suppress(Exception):
                    os.remove(x)

        else:
            logging.info(f"Scratch files kept: {tempfiles}")


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as e:
        logging.error(str(e))
        traceback.print_exc()
        sys.exit(1)
