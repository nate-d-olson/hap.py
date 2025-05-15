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

"""
Module for processing VarScan2 variant caller output for somatic variants.
Provides feature extraction and filtering for somatic analysis.

Date:   2/10/2015
Author: Madis Rumming <mrumming@illumina.com>, Peter Krusche <pkrusche@illumina.com>
"""

import logging
import re
from typing import Any, Dict, List, Optional, Tuple, Union, cast

import pandas
from Tools.vcfextract import extractHeadersJSON, vcfExtract


def extractVarscan2SNVFeatures(
    vcfname: str, tag: str, avg_depth: Optional[Dict[str, float]] = None
) -> pandas.DataFrame:
    """Return a data frame with features collected from the given VCF, tagged by given type.

    Args:
        vcfname: Name of the VCF file
        tag: Type of variants
        avg_depth: Dictionary mapping chromosomes to average depths

    Returns:
        DataFrame with extracted features
    """
    records: List[Dict[str, Any]] = []

    if not avg_depth:
        logging.warning(
            "No average depths available, normalized depth features cannot be calculated"
        )

    hdrs = extractHeadersJSON(vcfname)

    # TODO could figure this out automatically
    nsn = "NORMAL"
    tsn = "TUMOR"
    n_sample = "S.1."
    t_sample = "S.2."

    logging.info(
        f"Normal sample name : {nsn} (prefix {n_sample}) / "
        f"tumour sample name : {tsn} (prefix {t_sample})"
    )

    features = [
        "CHROM",
        "POS",
        "REF",
        "ALT",
        "FILTER",
        "I.SSC",
        "I.GPV",
        "I.SPV",
        n_sample + "GT",
        t_sample + "GT",  # Genotype
        n_sample + "GQ",
        t_sample + "GQ",  # Genotype quality
        n_sample + "DP",
        t_sample + "DP",  # Read depth
        n_sample + "RD",
        t_sample + "RD",  # Reference depth
        n_sample + "AD",
        t_sample + "AD",  # Alternative depth
        n_sample + "FREQ",
        t_sample + "FREQ",  # Alt. frequence (FA in MuTect)
    ]

    has_warned: Dict[str, bool] = {}

    for vr in vcfExtract(vcfname, features):
        rec: Dict[str, Any] = {}
        for i, ff in enumerate(features):
            rec[ff] = vr[i]

        for q in [n_sample + "GT", t_sample + "GT"]:
            if q not in rec or rec[q] is None:
                rec[q] = "."
                if ("feat:" + q) not in has_warned:
                    logging.warning(f"Missing feature {q}")
                    has_warned["feat:" + q] = True

        # fix missing features
        for q in [
            n_sample + "GT",
            t_sample + "GT",
            n_sample + "GQ",
            t_sample + "GQ",
            n_sample + "DP",
            t_sample + "DP",
            n_sample + "AD",
            t_sample + "AD",
            n_sample + "RD",
            t_sample + "RD",
            n_sample + "FREQ",
            t_sample + "FREQ",
        ]:
            if q not in rec or rec[q] is None:
                rec[q] = 0
                if ("feat:" + q) not in has_warned:
                    logging.warning(f"Missing feature {q}")
                    has_warned["feat:" + q] = True
            else:
                if q.endswith("FREQ"):
                    try:
                        rec[q] = float(rec[q])
                    except ValueError:
                        rec[q] = float("NaN")

                else:
                    try:
                        rec[q] = int(rec[q])
                    except ValueError:
                        rec[q] = -1

        rec["tag"] = tag

        n_DP = float(rec[n_sample + "DP"])
        t_DP = float(rec[t_sample + "DP"])

        n_DP_ratio = 0.0
        t_DP_ratio = 0.0

        if avg_depth:
            chrom = rec["CHROM"]
            if chrom in avg_depth:
                n_DP_ratio = n_DP / float(avg_depth[chrom])
                t_DP_ratio = t_DP / float(avg_depth[chrom])
            elif chrom not in has_warned:
                logging.warning(f"Cannot normalize depths on {chrom}")
                has_warned[chrom] = True
        elif "DPnorm" not in has_warned:
            logging.warning("Cannot normalize depths.")
            has_warned["DPnorm"] = True

        n_allele_ref_count = rec[n_sample + "RD"]
        alleles_alt = rec["ALT"]

        if alleles_alt == ["."]:
            n_allele_alt_count = 0
        else:
            n_allele_alt_count = rec[n_sample + "AD"]

        if n_allele_alt_count + n_allele_ref_count == 0:
            n_allele_rate = 0.0
        else:
            n_allele_rate = n_allele_alt_count / float(
                n_allele_alt_count + n_allele_ref_count
            )

        t_allele_ref_count = rec[t_sample + "RD"]

        if alleles_alt == ["."]:
            t_allele_alt_count = 0
        else:
            t_allele_alt_count = rec[t_sample + "AD"]

        if t_allele_alt_count + t_allele_ref_count == 0:
            t_allele_rate = 0.0
        else:
            t_allele_rate = t_allele_alt_count / float(
                t_allele_alt_count + t_allele_ref_count
            )

        # Gather the computed data into a dict
        qrec = {
            "CHROM": rec["CHROM"],
            "POS": int(rec["POS"]),
            "REF": rec["REF"],
            "ALT": ",".join(rec["ALT"]),
            "FILTER": ",".join(rec["FILTER"]),
            "SSC": rec["I.SSC"],
            "GPV": rec["I.GPV"],
            "SPV": rec["I.SPV"],
            "N_DP": n_DP,
            "T_DP": t_DP,
            "N_DP_RATE": n_DP_ratio,
            "T_DP_RATE": t_DP_ratio,
            "N_GT": rec[n_sample + "GT"],
            "T_GT": rec[t_sample + "GT"],
            "N_GQ": rec[n_sample + "GQ"],
            "T_GQ": rec[t_sample + "GQ"],
            "N_AD": rec[n_sample + "AD"],
            "T_AD": rec[t_sample + "AD"],
            "N_FA": rec[n_sample + "FREQ"],
            "T_FA": rec[t_sample + "FREQ"],
            "N_ALT_RATE": n_allele_rate,
            "T_ALT_RATE": t_allele_rate,
            "tag": tag,
        }

        records.append(qrec)

    cols = [
        "CHROM",
        "POS",
        "REF",
        "ALT",
        "FILTER",
        "SSC",
        "GPV",
        "SPV",
        "N_DP",
        "T_DP",
        "N_DP_RATE",
        "T_DP_RATE",
        "N_GT",
        "T_GT",
        "N_GQ",
        "T_GQ",
        "N_AD",
        "T_AD",
        "N_FA",
        "T_FA",
        "N_ALT_RATE",
        "T_ALT_RATE",
        "tag",
    ]

    if records:
        df = pandas.DataFrame(records, columns=cols)
    else:
        df = pandas.DataFrame(columns=cols)

    return df


def extractVarscan2IndelFeatures(
    vcfname: str, tag: str, avg_depth: Optional[Dict[str, float]] = None
) -> pandas.DataFrame:
    """Return a data frame with INDEL features collected from the given VCF, tagged by given type.

    Args:
        vcfname: Name of the VCF file
        tag: Type of variants
        avg_depth: Dictionary mapping chromosomes to average depths

    Returns:
        DataFrame with extracted features
    """
    # For VarScan2, we use the same extraction method for SNVs and INDELs
    return extractVarscan2SNVFeatures(vcfname, tag, avg_depth)
