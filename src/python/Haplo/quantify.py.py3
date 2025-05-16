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
# 01/04/2015
#
# Process raw counts coming out of quantify

import os
import tempfile
import subprocess
import copy
import json
import logging
import Tools
import pipes
from typing import Dict, List, Union, Optional, Any, Tuple, Set, TextIO

from Tools.bcftools import runBcftools


def _locations_tmp_bed_file(locations: Union[str, List[str]]) -> str:
    """ Turn a list of locations into a bed file

    Args:
        locations: List of locations as strings or comma-separated string

    Returns:
        Path to temporary BED file

    Raises:
        Exception: For invalid location formats
    """
    if isinstance(locations, str):
        locations = locations.split(",")
    if not isinstance(locations, list):
        raise Exception(f"Invalid list of locations (must be str or list): {str(locations)}")

    llocations = []

    for l in locations:
        xchr, _, _pos = l.partition(":")
        start, _, end = _pos.partition("-")
        if not xchr:
            raise Exception(f"Invalid chromosome name in {str(l)}")
        try:
            start = int(start)
        except:
            start = 0

        try:
            end = int(end)
        except:
            end = 2 ** 31 - 1

        llocations.append([xchr, start, end])

    # setup temporary file for locations
    fd, tpath = tempfile.mkstemp(suffix=".bed")
    os.close(fd)

    with open(tpath, "w") as f:
        for l in llocations:
            f.write("%s\t%i\t%i\n" % tuple(l))

    return tpath


def run(args: Any) -> None:
    """Run comparison and create summary statistics

    Args:
        args: Parsed command line arguments
    """
    import Haplo.xcmp
    import Haplo.vcfeval
    import Haplo.scmp

    outfiles = {}

    outprefix = args.prefix
    if not outprefix:
        fd, outprefix = tempfile.mkstemp()
        os.close(fd)

    if args.type != "ALL":
        typelist = args.type.split(",")
    else:
        typelist = ["INDEL", "SNP", "COMPLEX"]

    logging.info("Variant types to process: %s" % str(typelist))

    outvcfs = []
    for t in typelist:
        t_outprefix = outprefix + "." + t
        if args.unhappy:
            outfiles[t] = u_unhappy(args.truth, args.query, args.ref,
                                  args.regions, args.regions_file,
                                  t_outprefix, t.lower(),
                                  args.usefiltered_truth, args.usefiltered_query)
        else:
            if args.gender == "auto" or args.gender == "none":
                logging.warning("Auto / none for gender selection are not supported. Using female.")
                is_male = False
            else:
                is_male = args.gender.lower() == "male"

            if args.engine == "xcmp":
                # Use C++ haplotype comparison
                outfiles[t] = u_happyc(args.truth, args.query, args.ref,
                                    args.regions, args.regions_file,
                                    outprefix, t, args.preprocessing,
                                    args.window, args.fixchr_truth, args.fixchr_query,
                                    args.scratch_prefix, args.feature_table,
                                    args.usefiltered_truth, args.usefiltered_query,
                                    args.leftshift, args.decompose,
                                    args.bcftools_norm,
                                    args.threads, args.engine, args.preserve_all_variants,
                                    args.write_vcf, args.output_vtc, args.output_vtc_max_size,
                                    args.lose, args.xcmp_enumeration_threshold, is_male,
                                    pass_only=args.pass_only, conf_truth=args.conf_truth,
                                    conf_query=args.conf_query, optimize=args.optimize,
                                    preprocess_truth=args.preprocess_truth,
                                    preprocess_query=args.preprocess_query,
                                    preprocess_window=args.preprocess_window,
                                    location_features=args.location_features,
                                    adjust_conf_regions=args.adjust_conf_regions,
                                    fp_bedfile=args.false_positives)
            elif args.engine == "vcfeval":
                # Use RTG vcfeval
                outfiles[t] = v_vcfeval(args.truth, args.query, args.ref,
                                     args.regions, args.regions_file,
                                     outprefix, t, args.preprocessing,
                                     args.window, args.fixchr_truth, args.fixchr_query,
                                     args.scratch_prefix,
                                     args.usefiltered_truth, args.usefiltered_query,
                                     args.threads, args.vcfeval_path, args.vcfeval_template,
                                     args.preserve_all_variants,
                                     args.write_vcf, args.output_vtc, args.output_vtc_max_size,
                                     feature_table=args.feature_table)
            # elif args.engine == "scmp-somatic" or \
            #         (args.engine == "scmp-distance"):
            #     pass
            else:
                raise Exception(f"Invalid engine name: {args.engine}")

        if args.write_vcf and os.path.exists(t_outprefix + ".vcf.gz"):
            outvcfs.append(t_outprefix + ".vcf.gz")

    _write_outfiles(outfiles, outprefix, typelist, args.writeCounts)
    if outvcfs:
        xvcf = outprefix + ".vcf.gz"
        # open and pipe to bgzip
        _merge_vcfs(outvcfs, xvcf)


def _make_cmdline(args: List[str]) -> str:
    """Make a command line from arguments

    Args:
        args: List of command line arguments

    Returns:
        Formatted command line string
    """
    qargs = []
    for a in args:
        if a.strip() != "|":
            qargs.append(pipes.quote(a))
        else:
            qargs.append("|")
    return " ".join(qargs)


def _merge_vcfs(vcfs: List[str], outvcf: str) -> None:
    """Merge VCFs

    Args:
        vcfs: List of VCF files to merge
        outvcf: Output VCF file path
    """
    if os.path.exists(outvcf):
        try:
            os.unlink(outvcf)
        except:
            pass

    cmd_line = ["bcftools", "concat", "-a"]
    cmd_line.extend(vcfs)
    cmd_line.extend(["-o", outvcf])

    cmd_line_str = _make_cmdline(cmd_line)
    logging.info(cmd_line_str)

    po = subprocess.Popen(cmd_line_str,
                          shell=True,
                          stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE,
                          universal_newlines=True)

    stdout, stderr = po.communicate()

    po.wait()

    return_code = po.returncode

    if return_code != 0:
        logging.error(f"bcftools concat error: {stderr}")
        raise Exception(f"Failed to concatenate {str(vcfs)}")

    # index vcf
    cmd_line = ["bcftools", "index", outvcf]

    cmd_line_str = _make_cmdline(cmd_line)
    logging.info(cmd_line_str)

    po = subprocess.Popen(cmd_line_str,
                          shell=True,
                          stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE,
                          universal_newlines=True)

    stdout, stderr = po.communicate()

    po.wait()

    return_code = po.returncode

    if return_code != 0:
        logging.error(f"bcftools index error: {stderr}")
        logging.warning(f"Failed to index {outvcf}")


def _write_outfiles(outfiles: Dict[str, Dict[str, Any]],
                   outprefix: str,
                   typelist: List[str],
                   writeCounts: bool) -> None:
    """Write output files

    Args:
        outfiles: Dictionary of output files by variant type
        outprefix: Output file prefix
        typelist: List of variant types
        writeCounts: Whether to write count metrics
    """
    # write CSV outputs
    of_summary = open(outprefix + ".summary.csv", "w")
    of_extended = None
    of_metrics = None

    file_header = False

    header_lines = []
    data = {}

    for t in typelist:
        if t in outfiles:
            try:
                if not file_header:
                    header_lines.append("#" + outfiles[t]["summary_header"])
                    file_header = True
                data[t] = outfiles[t]["summary_csv"]
            except:
                pass  # might not have all outputs

    for h in header_lines:
        of_summary.write(h + "\n")

    for t in sorted(data.keys()):
        for l in data[t].splitlines():
            if not l.startswith("Type"):
                of_summary.write(l + "\n")

    for t in typelist:
        try:
            if of_extended is None and writeCounts and "extended_csv" in outfiles[t]:
                of_extended = open(outprefix + ".extended.csv", "w")
                of_extended.write(outfiles[t]["extended_header"] + "\n")

            if of_extended and "extended_csv" in outfiles[t]:
                for l in outfiles[t]["extended_csv"].splitlines():
                    if not l.startswith("#"):
                        of_extended.write(l + "\n")

            if of_metrics is None and writeCounts and "metrics" in outfiles[t]:
                of_metrics = open(outprefix + ".metrics.json.gz", "wb")
                import gzip
                of_metrics = gzip.GzipFile(fileobj=of_metrics)

            if of_metrics and "metrics" in outfiles[t]:
                of_metrics.write(json.dumps(outfiles[t]["metrics"]).encode('utf-8'))
        except:
            pass  # might not have all outputs
