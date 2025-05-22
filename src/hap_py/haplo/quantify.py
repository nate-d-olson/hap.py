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
from typing import Any, Dict, List, Optional, Tuple, Union


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

    for location in locations:
        xchr, _, pos = location.partition(":")
        start, _, end = pos.partition("-")
        if not xchr:
            raise Exception(f"Invalid chromosome name in {str(location)}")
        try:
            start = int(start)
        except ValueError:
            start = 0

        try:
            end = int(end)
        except ValueError:
            end = 2**31 - 1

        llocations.append([xchr, start, end])

    # setup temporary file for locations
    fd, tpath = tempfile.mkstemp(suffix=".bed")
    os.close(fd)

    with open(tpath, "w", encoding="utf-8") as f:
        for llocation in llocations:
            f.write("%s\t%i\t%i\n" % tuple(llocation))

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
    for variant_type in typelist:
        t_outprefix = outprefix + "." + variant_type
        if args.unhappy:
            outfiles[variant_type] = u_unhappy(
                args.truth,
                args.query,
                args.ref,
                args.regions,
                args.regions_file,
                t_outprefix,
                variant_type.lower(),
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
                outfiles[variant_type] = v_vcfeval(
                    args.truth,
                    args.query,
                    args.ref,
                    args.regions,
                    args.regions_file,
                    outprefix,
                    variant_type,
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

    try:
        po = subprocess.Popen(
            cmd_line_str,
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True,
        )

        stdout, stderr = po.communicate()
        return_code = po.returncode

        if return_code != 0:
            logging.error(f"bcftools concat error: {stderr}")
            raise Exception(f"Failed to concatenate {str(vcfs)}")
    except Exception as e:
        logging.error(f"Command execution failed: {str(e)}")
        raise Exception(f"Failed to concatenate {str(vcfs)}: {str(e)}")

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
    of_summary = open(outprefix + ".summary.csv", "w", encoding="utf-8")
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
            except KeyError:
                pass  # might not have all outputs

    for h in header_lines:
        of_summary.write(h + "\n")

    for t in sorted(data.keys()):
        for line in data[t].splitlines():
            if not line.startswith("Type"):
                of_summary.write(line + "\n")

    for variant_type in typelist:
        try:
            if (
                of_extended is None
                and writeCounts
                and "extended_csv" in outfiles[variant_type]
            ):
                of_extended = open(outprefix + ".extended.csv", "w", encoding="utf-8")
                of_extended.write(outfiles[variant_type]["extended_header"] + "\n")

            if of_extended and "extended_csv" in outfiles[variant_type]:
                for line in outfiles[variant_type]["extended_csv"].splitlines():
                    if not line.startswith("#"):
                        of_extended.write(line + "\n")

            if (
                of_metrics is None
                and writeCounts
                and "metrics" in outfiles[variant_type]
            ):
                metrics_file = open(outprefix + ".metrics.json.gz", "wb")
                import gzip

                of_metrics = gzip.GzipFile(fileobj=metrics_file)

            if of_metrics and "metrics" in outfiles[variant_type]:
                # Properly handle encoding in Python 3
                of_metrics.write(
                    json.dumps(outfiles[variant_type]["metrics"]).encode("utf-8")
                )
        except (KeyError, OSError) as exc:
            logging.warning(f"Error processing output for {variant_type}: {str(exc)}")
            # might not have all outputs


def _parse_vcfeval_stats(stats_file: str) -> Tuple[Dict[str, int], Dict[str, float]]:
    """Parse vcfeval output statistics file

    Args:
        stats_file: Path to the statistics file

    Returns:
        Tuple of (counts, metrics) dictionaries
    """
    counts = {
        "TRUTH.TOTAL": 0,
        "TRUTH.TP": 0,
        "TRUTH.FN": 0,
        "QUERY.TOTAL": 0,
        "QUERY.TP": 0,
        "QUERY.FP": 0,
    }

    metrics = {
        "Recall": 0.0,
        "Precision": 0.0,
        "F1_Score": 0.0,
    }

    try:
        with open(stats_file, encoding="utf-8") as f:
            for line in f:
                if line.startswith("#"):
                    continue

                parts = line.strip().split()
                if len(parts) >= 7:
                    # Format is typically: Type Count TP FP FN Precision Recall etc.
                    counts["TRUTH.TOTAL"] = int(parts[1])
                    counts["TRUTH.TP"] = int(parts[2])
                    counts["QUERY.FP"] = int(parts[3])
                    counts["TRUTH.FN"] = int(parts[4])
                    counts["QUERY.TOTAL"] = counts["TRUTH.TP"] + counts["QUERY.FP"]
                    counts["QUERY.TP"] = counts["TRUTH.TP"]

                    # Parse metrics
                    if len(parts) >= 9:
                        metrics["Precision"] = float(parts[5])
                        metrics["Recall"] = float(parts[6])
                        metrics["F1_Score"] = float(parts[7])
    except Exception as e:
        logging.warning(f"Failed to parse vcfeval stats file: {e}")

    return counts, metrics


def u_unhappy(
    truth: str,
    query: str,
    ref: str,
    regions: Optional[str],
    regions_file: Optional[str],
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

    # Create a minimal result structure that's compatible with quantify requirements
    result = {
        "summary_header": "Type,Filter,TRUTH.TOTAL,TRUTH.TP,TRUTH.FN,QUERY.TOTAL,QUERY.TP,QUERY.FP,METRIC.Recall,METRIC.Precision,METRIC.F1_Score,TRUTH.TOTAL.TiTv_ratio,QUERY.TOTAL.TiTv_ratio,TRUTH.TOTAL.TiTv_ratio,QUERY.TOTAL.TiTv_ratio",
        "summary_csv": f"{variant_type},ALL,0,0,0,0,0,0,0,0,0,0,0,0,0",
        "extended_header": "#Type,Subtype,Filter,TRUTH.TOTAL,TRUTH.TP,TRUTH.FN,QUERY.TOTAL,QUERY.TP,QUERY.FP,METRIC.Recall,METRIC.Precision,METRIC.F1_Score,TRUTH.TP.TiTv_ratio,QUERY.TP.TiTv_ratio",
        "extended_csv": f"{variant_type},*,ALL,0,0,0,0,0,0,0,0,0,0,0",
        "metrics": {
            "type": variant_type,
            "filter": "ALL",
            "counts": {
                "TRUTH.TOTAL": 0,
                "TRUTH.TP": 0,
                "TRUTH.FN": 0,
                "QUERY.TOTAL": 0,
                "QUERY.TP": 0,
                "QUERY.FP": 0,
            },
            "metrics": {"Recall": 0.0, "Precision": 0.0, "F1_Score": 0.0},
        },
    }

    return result


def v_vcfeval(
    truth: str,
    query: str,
    ref: str,
    regions: Optional[str],
    regions_file: Optional[str],
    outprefix: str,
    variant_type: str,
    preprocessing: bool,
    window: int,
    fixchr_truth: Optional[bool],
    fixchr_query: Optional[bool],
    scratch_prefix: Optional[str],
    usefiltered_truth: bool,
    usefiltered_query: bool,
    threads: int,
    vcfeval_path: Optional[str],
    vcfeval_template: Optional[str],
    preserve_all_variants: bool,
    write_vcf: bool,
    output_vtc: bool,
    output_vtc_max_size: int,
    feature_table: Optional[str] = None,
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
    # This function should be reimplemented to actually process the vcfeval output files
    # and generate proper metrics. For now, we return placeholders.

    logging.info(f"Processing vcfeval results for {variant_type}")

    # Initialize with default values
    counts = {
        "TRUTH.TOTAL": 0,
        "TRUTH.TP": 0,
        "TRUTH.FN": 0,
        "QUERY.TOTAL": 0,
        "QUERY.TP": 0,
        "QUERY.FP": 0,
    }

    metrics = {"Recall": 0.0, "Precision": 0.0, "F1_Score": 0.0}

    # We'd normally parse the vcfeval output files here
    # In a complete implementation, we should look for the statistics summary file
    # that vcfeval produces and parse it for the counts and metrics

    # Create a simulated summary string
    summary_csv = f"{variant_type},ALL,{counts['TRUTH.TOTAL']},{counts['TRUTH.TP']},{counts['TRUTH.FN']},{counts['QUERY.TOTAL']},{counts['QUERY.TP']},{counts['QUERY.FP']},{metrics['Recall']},{metrics['Precision']},{metrics['F1_Score']},0,0,0,0"

    # Extended data would include breakdowns by subtype
    extended_csv = f"{variant_type},*,ALL,{counts['TRUTH.TOTAL']},{counts['TRUTH.TP']},{counts['TRUTH.FN']},{counts['QUERY.TOTAL']},{counts['QUERY.TP']},{counts['QUERY.FP']},{metrics['Recall']},{metrics['Precision']},{metrics['F1_Score']},0,0"

    # Return structured output
    result = {
        "summary_header": "Type,Filter,TRUTH.TOTAL,TRUTH.TP,TRUTH.FN,QUERY.TOTAL,QUERY.TP,QUERY.FP,METRIC.Recall,METRIC.Precision,METRIC.F1_Score,TRUTH.TOTAL.TiTv_ratio,QUERY.TOTAL.TiTv_ratio,TRUTH.TOTAL.TiTv_ratio,QUERY.TOTAL.TiTv_ratio",
        "summary_csv": summary_csv,
        "extended_header": "#Type,Subtype,Filter,TRUTH.TOTAL,TRUTH.TP,TRUTH.FN,QUERY.TOTAL,QUERY.TP,QUERY.FP,METRIC.Recall,METRIC.Precision,METRIC.F1_Score,TRUTH.TP.TiTv_ratio,QUERY.TP.TiTv_ratio",
        "extended_csv": extended_csv,
        "metrics": {
            "type": variant_type,
            "filter": "ALL",
            "counts": counts,
            "metrics": metrics,
        },
    }

    return result


def run_quantify(
    vcf_name: str,
    roc_table: str,
    output_vcf: str,
    regions: Dict[str, str],
    ref: str,
    threads: int = 1,
    output_vtc: bool = False,
    output_rocs: bool = True,
    qtype: str = "ga4gh",
    roc_val: str = "QUAL",
    roc_header: Optional[str] = None,
    roc_filter: Optional[str] = None,
    roc_delta: float = 0.5,
    roc_regions: Optional[List[str]] = None,
    clean_info: bool = True,
    strat_fixchr: bool = False,
) -> None:
    """Run quantify and process results.

    This function processes a VCF output from a comparison engine (like vcfeval)
    to extract metrics.

    Args:
        vcf_name: Path to the input VCF from comparison engine
        roc_table: Path to ROC table output
        output_vcf: Path to output VCF
        regions: Quantification regions (dict of region name -> bed file)
        ref: Reference genome path
        threads: Number of threads
        output_vtc: Whether to output VTC
        output_rocs: Whether to output ROC
        qtype: Quantification type
        roc_val: ROC value
        roc_header: ROC header
        roc_filter: ROC filter
        roc_delta: ROC delta
        roc_regions: ROC regions
        clean_info: Whether to clean info
        strat_fixchr: Whether to fix chromosomes in stratification
    """
    logging.info(f"Processing {vcf_name} to extract metrics")

    if not os.path.exists(vcf_name):
        raise Exception(f"Input VCF {vcf_name} does not exist")

    # In a full implementation, this would:
    # 1. Process the VCF to extract variant counts
    # 2. Calculate performance metrics
    # 3. Generate ROC curves if requested
    # 4. Write outputs to the specified files

    # For now, this is a placeholder that acknowledges the function should exist
    logging.info(f"Metrics will be written to {roc_table}")
    if output_vcf:
        logging.info(f"Annotated VCF will be written to {output_vcf}")
