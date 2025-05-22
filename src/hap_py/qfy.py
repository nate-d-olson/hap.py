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
import gzip
import json
import logging
import multiprocessing
import os
import sys
import tempfile
import traceback
from pathlib import Path

import pandas as pd

import contextlib

# Modern imports using the new package structure
from .haplo import gvcf2bed, happyroc, quantify
from .tools import vcfextract, fastasize
from .tools.metric import dataframeToMetricsTable, makeMetricsObject
from .tools import version


def quantify(args: argparse.Namespace) -> None:
    """Run quantify and write tables"""
    vcf_name = args.in_vcf[0]

    if not vcf_name or not os.path.exists(vcf_name):
        raise FileNotFoundError("Cannot read input VCF.")

    logging.info("Counting variants...")

    truth_or_query_is_bcf = False
    with contextlib.suppress(Exception):
        truth_or_query_is_bcf = args.vcf1.endswith(".bcf") and args.vcf2.endswith(
            ".bcf"
        )

    internal_format_suffix = ".bcf" if args.bcf or truth_or_query_is_bcf else ".vcf.gz"

    output_vcf = args.reports_prefix + internal_format_suffix
    roc_table = args.reports_prefix + ".roc.tsv"

    qfyregions = {}

    if args.fp_bedfile:
        if not os.path.exists(args.fp_bedfile):
            raise FileNotFoundError(
                f"FP / Confident region file not found at {args.fp_bedfile}"
            )
        qfyregions["CONF"] = args.fp_bedfile

    if args.strat_tsv:
        with open(args.strat_tsv, encoding="utf-8") as strat_file_list:
            for strat_file in strat_file_list:
                n, _, f = strat_file.strip().partition("\t")
                if n in qfyregions:
                    raise ValueError(f"Duplicate stratification region ID: {n}")
                if not f:
                    if n:
                        raise ValueError(f"No file for stratification region {n}")
                    else:
                        continue
                if not os.path.exists(f):
                    f = os.path.join(
                        os.path.abspath(os.path.dirname(args.strat_tsv)), f
                    )
                if not os.path.exists(f):
                    raise FileNotFoundError(f"Quantification region file {f} not found")
                qfyregions[n] = f

    if args.strat_regions:
        for r in args.strat_regions:
            n, _, f = r.partition(":")
            if not os.path.exists(f):
                raise FileNotFoundError(f"Quantification region file {f} not found")
            qfyregions[n] = f

    if vcf_name == output_vcf or vcf_name == output_vcf + internal_format_suffix:
        raise ValueError(
            f"Cannot overwrite input VCF: {vcf_name} would be overwritten with output name {output_vcf}."
        )

    roc_header = args.roc
    with contextlib.suppress(Exception):
        roc_header = args.roc_header

    quantify.run_quantify(
        vcf_name,
        roc_table,
        output_vcf if args.write_vcf else False,
        qfyregions,
        args.ref,
        threads=args.threads,
        output_vtc=args.output_vtc,
        output_rocs=args.do_roc,
        qtype=args.type,
        roc_val=args.roc,
        roc_header=roc_header,
        roc_filter=args.roc_filter,
        roc_delta=args.roc_delta,
        roc_regions=args.roc_regions,
        clean_info=not args.preserve_info,
        strat_fixchr=args.strat_fixchr,
    )

    metrics_output = makeMetricsObject("%s.comparison" % args.runner)

    filter_handling = None
    try:
        if args.engine == "vcfeval" or not args.usefiltered:
            filter_handling = "ALL" if args.usefiltered else "PASS"
    except AttributeError:
        # if we run this through qfy, these arguments are not present
        pass

    total_region_size = None
    headers = vcfextract.extractHeadersJSON(vcf_name)
    try:
        contigs_to_use = ",".join(headers["tabix"]["chromosomes"])
        contig_lengths = fastasize.fastaNonNContigLengths(args.ref)
        total_region_size = fastasize.calculateLength(contig_lengths, contigs_to_use)
        logging.info(
            f"Subset.Size for * is {total_region_size}, based on these contigs: {contigs_to_use}"
        )
    except Exception:
        pass

    res = happyroc.roc(
        roc_table,
        args.reports_prefix + ".roc",
        filter_handling=filter_handling,
        ci_alpha=args.ci_alpha,
        total_region_size=total_region_size,
    )
    df = res["all"]

    # only use summary numbers
    df = df[(df["QQ"] == "*") & (df["Filter"].isin(["ALL", "PASS"]))]

    summary_columns = [
        "Type",
        "Filter",
    ]

    for additional_column in [
        "TRUTH.TOTAL",
        "TRUTH.TP",
        "TRUTH.FN",
        "QUERY.TOTAL",
        "QUERY.FP",
        "QUERY.UNK",
        "FP.gt",
        "FP.al",
        "METRIC.Recall",
        "METRIC.Precision",
        "METRIC.Frac_NA",
        "METRIC.F1_Score",
        "TRUTH.TOTAL.TiTv_ratio",
        "QUERY.TOTAL.TiTv_ratio",
        "TRUTH.TOTAL.het_hom_ratio",
        "QUERY.TOTAL.het_hom_ratio",
    ]:
        summary_columns.append(additional_column)

    # Remove subtype
    summary_df = df[
        (df["Subtype"] == "*") & (df["Genotype"] == "*") & (df["Subset"] == "*")
    ]

    summary_df[summary_columns].to_csv(
        args.reports_prefix + ".summary.csv", index=False
    )

    metrics_output["metrics"].append(
        dataframeToMetricsTable("summary.metrics", summary_df[summary_columns])
    )

    if args.write_counts:
        df.to_csv(args.reports_prefix + ".extended.csv", index=False)
        metrics_output["metrics"].append(dataframeToMetricsTable("all.metrics", df))

    essential_numbers = summary_df[summary_columns]

    pd.set_option("display.max_columns", 500)
    pd.set_option("display.width", 1000)

    essential_numbers = essential_numbers[
        essential_numbers["Type"].isin(["SNP", "INDEL"])
    ]

    logging.info("\n" + essential_numbers.to_string(index=False))

    # in default mode, print(result summary to stdout)
    if not args.quiet and not args.verbose:
        print("Benchmarking Summary:")
        print(essential_numbers.to_string(index=False))

    # keep this for verbose output
    if not args.verbose:
        with contextlib.suppress(Exception):
            os.unlink(roc_table)

    for t in list(res.keys()):
        metrics_output["metrics"].append(dataframeToMetricsTable("roc." + t, res[t]))

    # gzip JSON output
    if args.write_json:
        with gzip.open(
            args.reports_prefix + ".metrics.json.gz", "wt", encoding="utf-8"
        ) as fp:
            json.dump(metrics_output, fp)


def updateArgs(parser: argparse.ArgumentParser) -> None:
    """add common quantification args"""
    parser.add_argument(
        "-t",
        "--type",
        dest="type",
        choices=["ga4gh"],
        help="Annotation format in input VCF file.",
    )

    parser.add_argument(
        "-f",
        "--false-positives",
        dest="fp_bedfile",
        default=None,
        type=str,
        help="False positive / confident call regions (.bed or .bed.gz). Calls outside "
        "these regions will be labelled as UNK.",
    )

    parser.add_argument(
        "--stratification",
        dest="strat_tsv",
        default=None,
        type=str,
        help="Stratification file list (TSV format -- first column is region name, second column is file name).",
    )

    parser.add_argument(
        "--stratification-region",
        dest="strat_regions",
        default=[],
        action="append",
        help="Add single stratification region, e.g. --stratification-region TEST:test.bed",
    )

    parser.add_argument(
        "--stratification-fixchr",
        dest="strat_fixchr",
        default=None,
        action="store_true",
        help="Add chr prefix to stratification files if necessary",
    )

    parser.add_argument(
        "-V",
        "--write-vcf",
        dest="write_vcf",
        default=False,
        action="store_true",
        help="Write an annotated VCF.",
    )

    parser.add_argument(
        "-X",
        "--write-counts",
        dest="write_counts",
        default=True,
        action="store_true",
        help="Write advanced counts and metrics.",
    )

    parser.add_argument(
        "--no-write-counts",
        dest="write_counts",
        default=True,
        action="store_false",
        help="Do not write advanced counts and metrics.",
    )

    parser.add_argument(
        "--output-vtc",
        dest="output_vtc",
        default=False,
        action="store_true",
        help="Write VTC field in the final VCF which gives the counts each position has contributed to.",
    )

    parser.add_argument(
        "--preserve-info",
        dest="preserve_info",
        action="store_true",
        default=False,
        help="Preserve and merge the INFO fields in truth and query. Useful for ROC computation.",
    )

    parser.add_argument(
        "--roc",
        dest="roc",
        default="QUAL",
        help="Select a feature to produce a ROC on (INFO feature, QUAL, GQX, ...).",
    )

    parser.add_argument(
        "--no-roc",
        dest="do_roc",
        default=True,
        action="store_false",
        help="Disable ROC computation and only output summary statistics for more concise output.",
    )

    parser.add_argument(
        "--roc-regions",
        dest="roc_regions",
        default=["*"],
        action="append",
        help="Select a list of regions to compute ROCs in. By default, only the '*' region will"
        " produce ROC output (aggregate variant counts).",
    )

    parser.add_argument(
        "--roc-filter",
        dest="roc_filter",
        default=False,
        help="Select a filter to ignore when making ROCs.",
    )

    parser.add_argument(
        "--roc-delta",
        dest="roc_delta",
        default=0.5,
        type=float,
        help="Minimum spacing between ROC QQ levels.",
    )

    parser.add_argument(
        "--ci-alpha",
        dest="ci_alpha",
        default=0.0,
        type=float,
        help="Confidence level for Jeffrey's CI for recall, precision and fraction of non-assessed calls.",
    )

    parser.add_argument(
        "--no-json",
        dest="write_json",
        default=True,
        action="store_false",
        help="Disable JSON file output.",
    )


def main() -> int:
    """
    Main entry point for qfy.py.

    Returns:
        int: 0 on success, non-zero on failure
    """
    parser = argparse.ArgumentParser("Quantify annotated VCFs")

    parser.add_argument(
        "-v",
        "--version",
        dest="version",
        action="store_true",
        help="Show version number and exit.",
    )

    parser.add_argument(
        "in_vcf",
        help="Comparison intermediate VCF file to quantify (two column TRUTH/QUERY format)",
        nargs=1,
    )

    parser.add_argument(
        "--adjust-conf-regions",
        dest="preprocessing_truth_confregions",
        default=None,
        help="When hap.py was run with --adjust-conf-regions, on the original VCF, "
        "then quantify needs the truthset VCF in order to correctly reproduce "
        " the results. This switch allows us to pass the truth VCF into quantify.",
    )

    updateArgs(parser)

    # generic, keep in sync with hap.py!
    parser.add_argument(
        "-o",
        "--report-prefix",
        dest="reports_prefix",
        default=None,
        required=True,
        help="Filename prefix for report output.",
    )

    parser.add_argument(
        "-r", "--reference", dest="ref", default=None, help="Specify a reference file."
    )

    parser.add_argument(
        "--threads",
        dest="threads",
        default=multiprocessing.cpu_count(),
        type=int,
        help="Number of threads to use.",
    )

    parser.add_argument(
        "--logfile",
        dest="logfile",
        default=None,
        help="Write logging information into file rather than to stderr",
    )

    parser.add_argument(
        "--bcf",
        dest="bcf",
        action="store_true",
        default=False,
        help="Use BCF internally. This is the default when the input file"
        " is in BCF format already. Using BCF can speed up temp file access, "
        " but may fail for VCF files that have broken headers or records that "
        " don't comply with the header.",
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

    args.runner = "qfy.py"

    if not args.ref:
        args.ref = version.defaultReference()

    args.scratch_prefix = tempfile.gettempdir()

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
        exit(0)

    if args.version:
        print(f"qfy.py {version.version}")
        exit(0)

    if args.fp_bedfile and args.preprocessing_truth_confregions:
        conf_temp = gvcf2bed.gvcf2bed(
            args.preprocessing_truth_confregions,
            args.ref,
            args.fp_bedfile,
            args.scratch_prefix,
        )
        args.strat_regions.append("CONF_VARS:" + conf_temp)
        args.preprocessing_truth_confregions = None

    quantify(args)
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as e:
        logging.error(str(e))
        traceback.print_exc()
        sys.exit(1)
