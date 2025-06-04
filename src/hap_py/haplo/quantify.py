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
import shlex
import subprocess
import tempfile
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

import numpy as np
import pandas as pd

from .metrics_calculator import MetricsCalculator

# Set up logging
logger = logging.getLogger(__name__)


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
            qargs.append(shlex.quote(a))
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


def _parse_vcfeval_stats(
    stats_file: str,
) -> Tuple[Dict[str, int], Dict[str, float]]:
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
    from pathlib import Path

    logging.info(f"Processing vcfeval results for {variant_type}")

    # First, run vcfeval to generate the comparison files
    from .vcfeval import runVCFEval

    # Set up vcfeval arguments
    class VCFEvalArgs:
        def __init__(self):
            self.truth = truth
            self.query = query
            self.reference = ref
            self.regions = regions
            self.regions_file = regions_file
            self.scratch_prefix = scratch_prefix or "/tmp"
            self.output = outprefix
            self.preprocess_truth = preprocessing
            self.preprocess_query = preprocessing
            self.window = window
            self.fixchr_truth = fixchr_truth
            self.fixchr_query = fixchr_query
            self.usefiltered_truth = usefiltered_truth
            self.usefiltered_query = usefiltered_query
            self.threads = threads
            self.preserve_info = preserve_all_variants
            self.write_vcf = write_vcf
            self.write_counts = True
            self.engine_vcfeval_path = vcfeval_path
            self.engine_vcfeval_template = vcfeval_template

    vcfeval_args = VCFEvalArgs()

    # Run vcfeval
    try:
        result = runVCFEval(vcfeval_args)
        if result is None:
            raise RuntimeError("vcfeval failed to produce results")
    except Exception as e:
        logging.error(f"vcfeval failed: {e}")
        # Return empty results if vcfeval fails
        return _create_empty_vcfeval_result(variant_type)

    # Process vcfeval output files to generate quantify metrics
    output_path = Path(outprefix)
    tp_vcf = output_path.with_suffix(".vcf.gz")
    fn_vcf = str(output_path) + ".fn.vcf.gz"
    fp_vcf = str(output_path) + ".fp.vcf.gz"

    # Parse vcfeval summary file if it exists
    summary_file = Path(
        str(output_path).replace(".vcf", "") + "_vcfeval" + "/summary.txt"
    )
    counts = _parse_vcfeval_summary(summary_file) if summary_file.exists() else {}

    # If summary parsing failed, count variants manually
    if not counts:
        counts = _count_variants_from_vcfs(tp_vcf, fn_vcf, fp_vcf)

    # Calculate metrics using the MetricsCalculator
    metrics_calc = MetricsCalculator.calculate_basic_metrics(
        tp_count=counts.get("TRUTH.TP", 0),
        fp_count=counts.get("QUERY.FP", 0),
        fn_count=counts.get("TRUTH.FN", 0),
        total_truth=counts.get("TRUTH.TOTAL", None),
        total_query=counts.get("QUERY.TOTAL", None),
    )

    # Generate ROC table
    roc_table_path = str(output_path) + ".roc.tsv"
    _generate_roc_table(tp_vcf, roc_table_path, variant_type)

    # Create CSV output strings
    summary_csv = f"{variant_type},ALL,{counts.get('TRUTH.TOTAL', 0)},{counts.get('TRUTH.TP', 0)},{counts.get('TRUTH.FN', 0)},{counts.get('QUERY.TOTAL', 0)},{counts.get('QUERY.TP', 0)},{counts.get('QUERY.FP', 0)},{metrics_calc.recall:.6f},{metrics_calc.precision:.6f},{metrics_calc.f1_score:.6f},0,0,0,0"

    extended_csv = f"{variant_type},*,ALL,{counts.get('TRUTH.TOTAL', 0)},{counts.get('TRUTH.TP', 0)},{counts.get('TRUTH.FN', 0)},{counts.get('QUERY.TOTAL', 0)},{counts.get('QUERY.TP', 0)},{counts.get('QUERY.FP', 0)},{metrics_calc.recall:.6f},{metrics_calc.precision:.6f},{metrics_calc.f1_score:.6f},0,0"

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
            "metrics": {
                "Recall": metrics_calc.recall,
                "Precision": metrics_calc.precision,
                "F1_Score": metrics_calc.f1_score,
            },
        },
    }

    return result


def _create_empty_vcfeval_result(variant_type: str) -> Dict[str, Any]:
    """Create empty result when vcfeval fails."""
    counts = {
        "TRUTH.TOTAL": 0,
        "TRUTH.TP": 0,
        "TRUTH.FN": 0,
        "QUERY.TOTAL": 0,
        "QUERY.TP": 0,
        "QUERY.FP": 0,
    }

    metrics = {"Recall": 0.0, "Precision": 0.0, "F1_Score": 0.0}

    summary_csv = f"{variant_type},ALL,0,0,0,0,0,0,0.0,0.0,0.0,0,0,0,0"
    extended_csv = f"{variant_type},*,ALL,0,0,0,0,0,0,0.0,0.0,0.0,0,0"

    return {
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


def _parse_vcfeval_summary(summary_file: Path) -> Dict[str, int]:
    """Parse vcfeval summary.txt file to extract counts."""
    counts = {}

    try:
        with open(summary_file) as f:
            for line in f:
                line = line.strip()
                if line.startswith("Threshold"):
                    continue

                # Look for lines like: "None         12345   123   567   8901   234   345   0.9567   0.8912   0.9234"
                parts = line.split()
                if len(parts) >= 10 and parts[0] in ["None", "PASS"]:
                    try:
                        # RTG vcfeval format: threshold tp_baseline fp fn_baseline tp_call fp_call precision recall f_measure
                        counts["TRUTH.TP"] = int(parts[1])
                        counts["QUERY.FP"] = int(parts[2])
                        counts["TRUTH.FN"] = int(parts[3])
                        counts["QUERY.TP"] = int(parts[4])
                        # Calculate totals
                        counts["TRUTH.TOTAL"] = counts["TRUTH.TP"] + counts["TRUTH.FN"]
                        counts["QUERY.TOTAL"] = counts["QUERY.TP"] + counts["QUERY.FP"]
                        break
                    except (ValueError, IndexError):
                        continue

    except Exception as e:
        logging.warning(f"Failed to parse vcfeval summary file {summary_file}: {e}")

    return counts


def _count_variants_from_vcfs(tp_vcf: str, fn_vcf: str, fp_vcf: str) -> Dict[str, int]:
    """Count variants from VCF files if summary parsing fails."""
    counts = {
        "TRUTH.TP": 0,
        "TRUTH.FN": 0,
        "QUERY.TP": 0,
        "QUERY.FP": 0,
    }

    # Count TP variants
    if os.path.exists(tp_vcf):
        try:
            import pysam

            with pysam.VariantFile(tp_vcf) as vcf:
                tp_count = sum(1 for _ in vcf)
                counts["TRUTH.TP"] = tp_count
                counts["QUERY.TP"] = tp_count  # TP count is same for both
        except Exception as e:
            logging.warning(f"Failed to count variants in {tp_vcf}: {e}")

    # Count FN variants
    if os.path.exists(fn_vcf):
        try:
            import pysam

            with pysam.VariantFile(fn_vcf) as vcf:
                counts["TRUTH.FN"] = sum(1 for _ in vcf)
        except Exception as e:
            logging.warning(f"Failed to count variants in {fn_vcf}: {e}")

    # Count FP variants
    if os.path.exists(fp_vcf):
        try:
            import pysam

            with pysam.VariantFile(fp_vcf) as vcf:
                counts["QUERY.FP"] = sum(1 for _ in vcf)
        except Exception as e:
            logging.warning(f"Failed to count variants in {fp_vcf}: {e}")

    # Calculate totals
    counts["TRUTH.TOTAL"] = counts["TRUTH.TP"] + counts["TRUTH.FN"]
    counts["QUERY.TOTAL"] = counts["QUERY.TP"] + counts["QUERY.FP"]

    return counts


def _generate_roc_table(tp_vcf: str, roc_table_path: str, variant_type: str) -> None:
    """Generate ROC table from TP VCF file."""
    import pysam

    try:
        # Create basic ROC table structure
        roc_data = []

        if os.path.exists(tp_vcf):
            # Read quality scores from TP VCF
            with pysam.VariantFile(tp_vcf) as vcf:
                quality_scores = []
                for record in vcf:
                    # Extract quality score
                    if record.qual is not None:
                        quality_scores.append(float(record.qual))
                    elif "GQ" in record.format and len(record.samples) > 0:
                        # Try to get genotype quality
                        sample = list(record.samples.values())[0]
                        if "GQ" in sample and sample["GQ"] is not None:
                            quality_scores.append(float(sample["GQ"]))
                    else:
                        quality_scores.append(0.0)

                # Generate ROC points at different thresholds
                if quality_scores:
                    thresholds = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
                    for threshold in thresholds:
                        tp_count = sum(1 for q in quality_scores if q >= threshold)
                        fp_count = (
                            0  # Simplified - would need FP VCF for real calculation
                        )

                        # Calculate precision/recall (simplified)
                        total_tp = len(quality_scores)
                        precision = tp_count / max(tp_count + fp_count, 1)
                        recall = tp_count / max(total_tp, 1)

                        roc_data.append(
                            {
                                "Score": threshold,
                                "Threshold": f">= {threshold}",
                                "True-pos": tp_count,
                                "False-pos": fp_count,
                                "False-neg": total_tp - tp_count,
                                "Precision": precision,
                                "Recall": recall,
                                "F-measure": 2
                                * precision
                                * recall
                                / max(precision + recall, 1e-10),
                                "Filter": "ALL",
                                "Subset": variant_type,
                            }
                        )

        # Create DataFrame and save
        if roc_data:
            roc_df = pd.DataFrame(roc_data)
        else:
            # Create empty ROC table
            roc_df = pd.DataFrame(
                {
                    "Score": [0],
                    "Threshold": [">= 0"],
                    "True-pos": [0],
                    "False-pos": [0],
                    "False-neg": [0],
                    "Precision": [0.0],
                    "Recall": [0.0],
                    "F-measure": [0.0],
                    "Filter": ["ALL"],
                    "Subset": [variant_type],
                }
            )

        # Save ROC table
        roc_df.to_csv(roc_table_path, sep="\t", index=False)
        logger.info(f"ROC table saved: {roc_table_path}")

    except Exception as e:
        logger.warning(f"Failed to generate ROC table: {e}")
        # Create minimal ROC table
        roc_df = pd.DataFrame(
            {
                "Score": [0],
                "Threshold": [">= 0"],
                "True-pos": [0],
                "False-pos": [0],
                "False-neg": [0],
                "Precision": [0.0],
                "Recall": [0.0],
                "F-measure": [0.0],
                "Filter": ["ALL"],
                "Subset": [variant_type],
            }
        )
        roc_df.to_csv(roc_table_path, sep="\t", index=False)


def run_quantify(
    vcf_name: str,
    roc_table: str,
    output_vcf: Union[str, bool] = False,
    regions: Optional[Dict[str, str]] = None,
    reference: Optional[str] = None,
    threads: int = 1,
    output_vtc: bool = False,
    output_rocs: bool = True,
    qtype: str = "xcmp",
    roc_val: str = "QUAL",
    roc_header: Optional[str] = None,
    roc_filter: Optional[str] = None,
    roc_delta: float = 0.001,
    roc_regions: Optional[List[str]] = None,
    clean_info: bool = True,
    strat_fixchr: bool = False,
) -> None:
    """
    Main quantify function that processes VCF comparison results.

    This function serves as the primary entry point for quantification analysis,
    processing VCF files to generate metrics, ROC curves, and stratified results.

    Args:
        vcf_name: Path to the input VCF file (comparison results)
        roc_table: Output path for ROC table
        output_vcf: Output VCF path or False to skip VCF output
        regions: Dictionary of region names to BED file paths for stratification
        reference: Path to reference genome FASTA file
        threads: Number of threads to use (currently not implemented)
        output_vtc: Whether to output variant truth categories
        output_rocs: Whether to generate ROC curves
        qtype: Quantification type (e.g., 'xcmp')
        roc_val: Field to use for ROC analysis (default: 'QUAL')
        roc_header: Custom header for ROC output
        roc_filter: Filter expression for ROC analysis
        roc_delta: Delta value for ROC curve generation
        roc_regions: List of regions for ROC analysis
        clean_info: Whether to clean INFO fields in output
        strat_fixchr: Whether to fix chromosome names in stratification

    Raises:
        FileNotFoundError: If input VCF file doesn't exist
        ValueError: For invalid parameters
    """
    logger.info(f"Starting quantify analysis on {vcf_name}")

    # Validate inputs
    if not os.path.exists(vcf_name):
        raise FileNotFoundError(f"Input VCF file not found: {vcf_name}")

    # Initialize regions if not provided
    if regions is None:
        regions = {}

    # Import the QuantifyEngine here to avoid circular imports
    from .python_quantify import QuantifyEngine

    try:
        # Create the quantify engine
        engine = QuantifyEngine(
            truth_vcf=vcf_name,  # For xcmp results, the input is already compared
            query_vcf=vcf_name,
            reference=reference,
            regions=None,  # Will handle regions separately
            apply_filters=not clean_info,
            output_vtc=output_vtc,
        )

        # Process the VCF file to extract variant metrics
        logger.info("Processing VCF file...")
        df = engine.process_vcf()

        if df is None or df.empty:
            logger.warning("No variants found in VCF file")
            # Create empty ROC table
            _create_empty_roc_table(roc_table)
            return

        # Apply stratification if regions are provided
        if regions:
            logger.info(f"Applying stratification with {len(regions)} region(s)")
            stratified_results = {}

            for region_name, bed_file in regions.items():
                if os.path.exists(bed_file):
                    logger.info(f"Processing region: {region_name}")
                    region_df = engine.apply_bed_stratification(df, bed_file)
                    stratified_results[region_name] = region_df
                else:
                    logger.warning(f"Region BED file not found: {bed_file}")

        # Generate ROC curves if requested
        if output_rocs and roc_table:
            logger.info("Generating ROC curves...")
            _generate_roc_curves(
                df=df,
                roc_table=roc_table,
                roc_val=roc_val,
                roc_delta=roc_delta,
                roc_filter=roc_filter,
                roc_regions=roc_regions,
                regions=regions,
            )

        # Write output VCF if requested
        if output_vcf and isinstance(output_vcf, str):
            logger.info(f"Writing output VCF to {output_vcf}")
            engine.write_output_vcf(df, output_vcf)

        logger.info("Quantify analysis completed successfully")

    except Exception as e:
        logger.error(f"Error during quantification: {e}")
        # Create empty ROC table on error to maintain expected output
        _create_empty_roc_table(roc_table)
        raise


def _create_empty_roc_table(roc_table: str) -> None:
    """Create an empty ROC table with proper headers."""
    empty_df = pd.DataFrame(
        columns=[
            "Type",
            "Filter",
            "TRUTH.TOTAL",
            "TRUTH.TP",
            "TRUTH.FN",
            "QUERY.TOTAL",
            "QUERY.TP",
            "QUERY.FP",
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
        ]
    )
    empty_df.to_csv(roc_table, sep="\t", index=False)


def _generate_roc_curves(
    df: pd.DataFrame,
    roc_table: str,
    roc_val: str = "QUAL",
    roc_delta: float = 0.001,
    roc_filter: Optional[str] = None,
    roc_regions: Optional[List[str]] = None,
    regions: Optional[Dict[str, str]] = None,
) -> None:
    """Generate ROC curves and write to table."""
    from .metrics_calculator import MetricsCalculator

    calculator = MetricsCalculator()

    # Generate thresholds for ROC curve
    if roc_val in df.columns:
        qual_values = df[roc_val].dropna()
        if len(qual_values) > 0:
            min_qual = qual_values.min()
            max_qual = qual_values.max()
            thresholds = np.arange(min_qual, max_qual + roc_delta, roc_delta)
        else:
            thresholds = [0.0]
    else:
        logger.warning(
            f"ROC field '{roc_val}' not found in VCF, using default thresholds"
        )
        thresholds = [0.0]

    roc_rows = []

    # Process different variant types
    variant_types = ["INDEL", "SNP"]
    if "Type" in df.columns:
        variant_types = df["Type"].unique().tolist()

    for variant_type in variant_types:
        type_df = (
            df[df.get("Type", "SNP") == variant_type] if "Type" in df.columns else df
        )

        if type_df.empty:
            continue

        for threshold in thresholds[
            :100
        ]:  # Limit to first 100 thresholds for performance
            # Filter by quality threshold
            filtered_df = (
                type_df[type_df.get(roc_val, 0) >= threshold]
                if roc_val in type_df.columns
                else type_df
            )

            # Calculate metrics
            metrics = calculator.calculate_summary_metrics(filtered_df)

            roc_row = {
                "Type": variant_type,
                "Filter": "ALL",
                "TRUTH.TOTAL": metrics.get("truth_total", 0),
                "TRUTH.TP": metrics.get("truth_tp", 0),
                "TRUTH.FN": metrics.get("truth_fn", 0),
                "QUERY.TOTAL": metrics.get("query_total", 0),
                "QUERY.TP": metrics.get("query_tp", 0),
                "QUERY.FP": metrics.get("query_fp", 0),
                "FP.gt": metrics.get("fp_gt", 0),
                "FP.al": metrics.get("fp_al", 0),
                "METRIC.Recall": metrics.get("recall", 0.0),
                "METRIC.Precision": metrics.get("precision", 0.0),
                "METRIC.Frac_NA": metrics.get("frac_na", 0.0),
                "METRIC.F1_Score": metrics.get("f1_score", 0.0),
                "TRUTH.TOTAL.TiTv_ratio": metrics.get("truth_titv", 0.0),
                "QUERY.TOTAL.TiTv_ratio": metrics.get("query_titv", 0.0),
                "TRUTH.TOTAL.het_hom_ratio": metrics.get("truth_het_hom", 0.0),
                "QUERY.TOTAL.het_hom_ratio": metrics.get("query_het_hom", 0.0),
            }
            roc_rows.append(roc_row)

    # Create ROC DataFrame and save
    if roc_rows:
        roc_df = pd.DataFrame(roc_rows)
        roc_df.to_csv(roc_table, sep="\t", index=False)
        logger.info(f"ROC table written to {roc_table}")
    else:
        _create_empty_roc_table(roc_table)
        logger.warning("No ROC data generated, created empty table")


# Backward compatibility function
def quantify(*args, **kwargs):
    """
    Backward compatibility wrapper for quantify function.

    This maintains compatibility with existing code that calls quantify()
    directly instead of run_quantify().
    """
    logger.warning("quantify() function is deprecated. Use run_quantify() instead.")
    return run_quantify(*args, **kwargs)
