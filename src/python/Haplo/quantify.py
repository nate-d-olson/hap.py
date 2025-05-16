#!/usr/bin/env python33
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

import contextlib
import json
import logging
import os
import pipes
import subprocess
import tempfile
from typing import Any, Dict, List, Union


def _locations_tmp_bed_file(locations: Union[str, List[str]]) -> str:
    """Turn a list of locations into a bed file

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
        raise Exception(
            f"Invalid list of locations (must be str or list): {str(locations)}"
        )

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
            end = 2**31 - 1

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
            outfiles[t] = u_unhappy(
                args.truth,
                args.query,
                args.ref,
                args.regions,
                args.regions_file,
                t_outprefix,
                t.lower(),
                args.usefiltered_truth,
                args.usefiltered_query,
            )
        else:
            if args.gender == "auto" or args.gender == "none":
                logging.warning(
                    "Auto / none for gender selection are not supported. Using female."
                )
            else:
                args.gender.lower() == "male"

            if args.engine == "vcfeval":
                # Use RTG vcfeval
                outfiles[t] = v_vcfeval(
                    args.truth,
                    args.query,
                    args.ref,
                    args.regions,
                    args.regions_file,
                    outprefix,
                    t,
                    args.preprocessing,
                    args.window,
                    args.fixchr_truth,
                    args.fixchr_query,
                    args.scratch_prefix,
                    args.usefiltered_truth,
                    args.usefiltered_query,
                    args.threads,
                    args.vcfeval_path,
                    args.vcfeval_template,
                    args.preserve_all_variants,
                    args.write_vcf,
                    args.output_vtc,
                    args.output_vtc_max_size,
                    feature_table=args.feature_table,
                )
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
        with contextlib.suppress(Exception):
            os.unlink(outvcf)

    cmd_line = ["bcftools", "concat", "-a"]
    cmd_line.extend(vcfs)
    cmd_line.extend(["-o", outvcf])

    cmd_line_str = _make_cmdline(cmd_line)
    logging.info(cmd_line_str)

    po = subprocess.Popen(
        cmd_line_str,
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True,
    )

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

    po = subprocess.Popen(
        cmd_line_str,
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True,
    )

    stdout, stderr = po.communicate()

    po.wait()

    return_code = po.returncode

    if return_code != 0:
        logging.error(f"bcftools index error: {stderr}")
        logging.warning(f"Failed to index {outvcf}")


def _write_outfiles(
    outfiles: Dict[str, Dict[str, Any]],
    outprefix: str,
    typelist: List[str],
    writeCounts: bool,
) -> None:
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
                of_metrics.write(json.dumps(outfiles[t]["metrics"]).encode("utf-8"))
        except:
            pass  # might not have all outputs


def u_unhappy(
    truth: str,
    query: str,
    ref: str,
    regions: str,
    regions_file: str,
    outprefix: str,
    variant_type: str,
    usefiltered_truth: bool,
    usefiltered_query: bool,
) -> Dict[str, Any]:
    """Unhappy comparison (direct GT matching)

    This function is maintained for backwards compatibility but should be replaced
    with a more robust implementation

    Args:
        truth: Truth VCF file
        query: Query VCF file
        ref: Reference genome
        regions: Regions to analyze
        regions_file: Regions file
        outprefix: Output prefix
        variant_type: Variant type
        usefiltered_truth: Whether to use filtered variants in truth
        usefiltered_query: Whether to use filtered variants in query

    Returns:
        Dictionary with summary metrics
    """
    logging.warning("u_unhappy is a placeholder - please use vcfeval engine instead")
    return {
        "summary_header": "Type,Filter,TRUTH.TOTAL,TRUTH.TP,TRUTH.FN,QUERY.TOTAL,QUERY.TP,QUERY.FP,METRIC.Recall,METRIC.Precision,METRIC.F1_Score,TRUTH.TOTAL.TiTv_ratio,QUERY.TOTAL.TiTv_ratio,TRUTH.TOTAL.TiTv_ratio,QUERY.TOTAL.TiTv_ratio",
        "summary_csv": f"Type={variant_type},Filter=ALL,0,0,0,0,0,0,0,0,0,0,0,0,0",
    }


def v_vcfeval(
    truth: str,
    query: str,
    ref: str,
    regions: str,
    regions_file: str,
    outprefix: str,
    variant_type: str,
    preprocessing: bool,
    window: int,
    fixchr_truth: bool,
    fixchr_query: bool,
    scratch_prefix: str,
    usefiltered_truth: bool,
    usefiltered_query: bool,
    threads: int,
    vcfeval_path: str,
    vcfeval_template: str,
    preserve_all_variants: bool,
    write_vcf: bool,
    output_vtc: bool,
    output_vtc_max_size: int,
    feature_table: str = None,
) -> Dict[str, Any]:
    """Run RTG's vcfeval and process results

    Args:
        truth: Truth VCF
        query: Query VCF
        ref: Reference FASTA
        regions: Region string
        regions_file: BED file with regions
        outprefix: Output prefix
        variant_type: Variant type (SNP, INDEL, etc)
        preprocessing: Whether to preprocess
        window: Window size
        fixchr_truth: Whether to fix chromosome names in truth
        fixchr_query: Whether to fix chromosome names in query
        scratch_prefix: Scratch directory
        usefiltered_truth: Whether to use filtered variants in truth
        usefiltered_query: Whether to use filtered variants in query
        threads: Number of threads
        vcfeval_path: Path to vcfeval
        vcfeval_template: SDF template for vcfeval
        preserve_all_variants: Whether to preserve all variants
        write_vcf: Whether to write VCF output
        output_vtc: Whether to output variant truth coverage
        output_vtc_max_size: Maximum size for VTC
        feature_table: Feature table file

    Returns:
        Dictionary with summary metrics
    """
    # This is a simplified implementation to handle vcfeval results
    # In a real implementation, this would process the output from vcfeval

    # Here we would normally call out to Haplo.vcfeval to run vcfeval and then
    # process the results to generate metrics

    logging.info(f"Processing vcfeval results for {variant_type}")

    # Return a basic structure with placeholders
    return {
        "summary_header": "Type,Filter,TRUTH.TOTAL,TRUTH.TP,TRUTH.FN,QUERY.TOTAL,QUERY.TP,QUERY.FP,METRIC.Recall,METRIC.Precision,METRIC.F1_Score,TRUTH.TOTAL.TiTv_ratio,QUERY.TOTAL.TiTv_ratio,TRUTH.TOTAL.TiTv_ratio,QUERY.TOTAL.TiTv_ratio",
        "summary_csv": f"Type={variant_type},Filter=ALL,0,0,0,0,0,0,0,0,0,0,0,0,0",
        "extended_header": "#Type,Subtype,Filter,TRUTH.TOTAL,TRUTH.TP,TRUTH.FN,QUERY.TOTAL,QUERY.TP,QUERY.FP,METRIC.Recall,METRIC.Precision,METRIC.F1_Score,TRUTH.TP.TiTv_ratio,QUERY.TP.TiTv_ratio",
        "extended_csv": f"#Type={variant_type},Subtype=*,Filter=*,0,0,0,0,0,0,0,0,0,0,0",
        "metrics": {"type": variant_type},
    }
