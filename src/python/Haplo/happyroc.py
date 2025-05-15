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
# Calculate diploid SNP and Indel ROCs
#
# This is a little different from the somatic case that is handled in Tools/roc.py
# -- we don't want to read in any of the data since we have a lot more of it.
#

import os
import re
import logging
import pandas as pd
import numpy as np
from typing import Dict, List, Tuple, Optional, Any, Union, Set

from Tools import ci

RESULT_ALLCOLUMNS = ["Type",
                     "Subtype",
                     "Subset",
                     "Filter",
                     "Genotype",
                     "QQ.Field",
                     "QQ",
                     "METRIC.Recall",
                     "METRIC.Precision",
                     "METRIC.Frac_NA",
                     "METRIC.F1_Score",
                     "FP.gt",
                     "FP.al",
                     "Subset.Size",
                     "Subset.IS_CONF.Size",
                     "Subset.Level"]

RESULT_ALLDTYPES = [str,
                    str,
                    str,
                    str,
                    str,
                    str,
                    str,
                    float,
                    float,
                    float,
                    float,
                    float,
                    float,
                    float,
                    float,
                    int]

for count_type in ["TRUTH.TOTAL", "TRUTH.TP", "TRUTH.FN",
                   "QUERY.TOTAL", "QUERY.TP", "QUERY.FP",
                   "QUERY.UNK"]:
    RESULT_ALLCOLUMNS.append(count_type)
    RESULT_ALLCOLUMNS.append(count_type + ".ti")
    RESULT_ALLCOLUMNS.append(count_type + ".tv")
    RESULT_ALLCOLUMNS.append(count_type + ".het")
    RESULT_ALLCOLUMNS.append(count_type + ".homalt")
    RESULT_ALLCOLUMNS.append(count_type + ".TiTv_ratio")
    RESULT_ALLCOLUMNS.append(count_type + ".het_hom_ratio")
    RESULT_ALLDTYPES += [float] * 7


def roc(roc_table: str, output_path: str,
        filter_handling: Optional[str] = None,
        ci_alpha: float = 0.05,
        total_region_size: Optional[float] = None) -> Dict[str, str]:
    """ Calculate SNP and indel ROC.

    Return a dictionary of variant types and corresponding files.

    Args:
        roc_table: Path to the input ROC table
        output_path: Path where output files should be written
        filter_handling: Can be None, "PASS" or "ALL" to filter rows based on the "Filter" column.
                        This is necessary because vcfeval doesn't preserve filter information
                        in GA4GH output mode, so we need to remove the corresponding rows here
        ci_alpha: Jeffrey's CI confidence level for recall, precision, na
        total_region_size: Correct Subset.Size for "*" region if a subset was selected in hap.py

    Returns:
        Dictionary mapping ROC types to file paths
    """
    result: Dict[str, List[Dict[str, str]]] = {}
    header = None
    
    with open(roc_table, 'r') as rt:
        for l in rt:
            l = l.strip()
            if not header:
                header = l.split("\t")
            else:
                rec = {}
                for k, v in zip(header, l.split("\t")):
                    rec[k] = v

                if filter_handling:
                    try:
                        if rec["Filter"] != filter_handling:
                            continue
                    except KeyError:
                        pass

                try:
                    if rec["Type"] in ["SNP", "INDEL"] \
                       and rec["Filter"] == "ALL" \
                       and rec["Subset"] == "*" \
                       and rec["Genotype"] == "*" \
                       and rec["Subtype"] == "*" \
                       and rec["QQ"] != "*":  # this is the ROC score field
                        roc_key = "Locations." + rec["Type"]
                        if roc_key not in result:
                            result[roc_key] = [rec]
                        else:
                            result[roc_key].append(rec)
                except KeyError:
                    pass

                try:
                    if rec["Type"] in ["SNP", "INDEL"] \
                       and rec["Filter"] == "PASS" \
                       and rec["Subset"] == "*" \
                       and rec["Genotype"] == "*" \
                       and rec["Subtype"] == "*" \
                       and rec["QQ"] != "*":  # this is the ROC score field
                        roc_key = "Locations." + rec["Type"] + ".PASS"
                        if roc_key not in result:
                            result[roc_key] = [rec]
                        else:
                            result[roc_key].append(rec)
                except KeyError:
                    pass

                try:
                    if rec["Type"] in ["SNP", "INDEL"] \
                       and rec["Filter"] == "SEL" \
                       and rec["Subset"] == "*" \
                       and rec["Genotype"] == "*" \
                       and rec["Subtype"] == "*" \
                       and rec["QQ"] != "*":  # this is the ROC score field
                        roc_key = "Locations." + rec["Type"] + ".SEL"
                        if roc_key not in result:
                            result[roc_key] = [rec]
                        else:
                            result[roc_key].append(rec)
                except KeyError:
                    pass

    file_result: Dict[str, str] = {}
    
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    
    for roc_key, rs in result.items():
        # Sort ROC scores in order of stringency
        try:
            rs = list(sorted(rs, key=lambda x: float(x["QQ"])))
        except (ValueError, KeyError):
            logging.error(f"Failed to sort ROC records for {roc_key}")
            continue

        # Output files
        filename = os.path.join(output_path, roc_key)
        roc_filename = filename + ".roc.tsv"
        file_result[roc_key] = roc_filename
        
        with open(roc_filename, "w") as f:
            cols_to_write = ["QQ", "METRIC.Recall", "METRIC.Precision", "TRUTH.TP", "TRUTH.FN", 
                            "QUERY.FP", "QUERY.UNK", "FP.gt", "FP.al"]
            
            f.write("\t".join(cols_to_write) + "\n")
            
            for r in rs:
                values = [r.get(c, "") for c in cols_to_write]
                f.write("\t".join(values) + "\n")

        # Calculate confidence intervals if requested
        if ci_alpha > 0:
            ci_filename = filename + ".roc.conf.tsv"
            try:
                calculate_confidence_intervals(rs, ci_filename, ci_alpha, total_region_size)
                file_result[roc_key + ".conf"] = ci_filename
            except Exception as e:
                logging.error(f"Failed to calculate confidence intervals for {roc_key}: {e}")

    return file_result


def calculate_confidence_intervals(
    rs: List[Dict[str, str]], 
    output_filename: str, 
    ci_alpha: float = 0.05, 
    total_region_size: Optional[float] = None
) -> None:
    """Calculate confidence intervals for ROC metrics
    
    Args:
        rs: List of ROC records
        output_filename: File to write results to
        ci_alpha: Jeffrey's CI confidence level
        total_region_size: Correct Subset.Size for "*" region if a subset was selected
    """
    with open(output_filename, "w") as f:
        # Write header
        cols = ["QQ", "METRIC.Recall", "METRIC.Precision", "METRIC.Frac_NA", 
                "RECALL.CONF.L", "RECALL.CONF.H", "PREC.CONF.L", "PREC.CONF.H",
                "NA.CONF.L", "NA.CONF.H"]
        f.write("\t".join(cols) + "\n")
        
        for r in rs:
            try:
                values = []
                # QQ value
                values.append(r["QQ"])
                
                # Original metrics
                recall = float(r["METRIC.Recall"])
                precision = float(r["METRIC.Precision"])
                frac_na = float(r.get("METRIC.Frac_NA", "0"))
                values.extend([str(recall), str(precision), str(frac_na)])
                
                # Calculate confidence intervals
                try:
                    truth_tp = float(r["TRUTH.TP"])
                    truth_fn = float(r["TRUTH.FN"])
                    query_fp = float(r["QUERY.FP"])
                    query_tp = float(r["QUERY.TP"])
                    query_unk = float(r.get("QUERY.UNK", "0"))
                    
                    # Recall CI
                    recall_ci = ci.jeffreys_interval(truth_tp, truth_tp + truth_fn, ci_alpha)
                    values.extend([str(recall_ci[0]), str(recall_ci[1])])
                    
                    # Precision CI
                    precision_ci = ci.jeffreys_interval(query_tp, query_tp + query_fp, ci_alpha)
                    values.extend([str(precision_ci[0]), str(precision_ci[1])])
                    
                    # NA fraction CI
                    if total_region_size is not None and total_region_size > 0:
                        # Adjust for selected subset
                        na_ci = ci.jeffreys_interval(query_unk, total_region_size, ci_alpha)
                    else:
                        # Use subset size from record
                        subset_size = float(r.get("Subset.Size", "0"))
                        if subset_size > 0:
                            na_ci = ci.jeffreys_interval(query_unk, subset_size, ci_alpha)
                        else:
                            na_ci = (0, 0)
                    
                    values.extend([str(na_ci[0]), str(na_ci[1])])
                except (ValueError, KeyError) as e:
                    # Use zeros if we can't calculate CIs
                    logging.warning(f"Failed to calculate confidence intervals: {e}")
                    values.extend(["0", "0", "0", "0", "0", "0"])
                
                f.write("\t".join(values) + "\n")
                
            except (ValueError, KeyError) as e:
                logging.error(f"Error processing record for CI calculation: {e}")
                continue
