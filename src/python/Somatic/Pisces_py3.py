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
Module for processing Pisces variant caller output for somatic variants.
Provides feature extraction and filtering for somatic analysis.
"""

import logging
from typing import Dict, Optional

import pandas
from Tools.vcfextract import extractHeaders, vcfExtract


def extractPiscesSNVFeatures(
    vcfname: str, tag: str, avg_depth: Optional[float] = None
) -> pandas.DataFrame:
    """Return a data frame with features collected from the given VCF, tagged by given type.

    Args:
        vcfname: Name of the VCF file
        tag: Type of variants
        avg_depth: Average chromosome depths from BAM file

    Returns:
        DataFrame with extracted features
    """
    features = [
        "CHROM",
        "POS",
        "REF",
        "ALT",
        "FILTER",
        "I.DP",
        "I.EVS",
        "S.1.GT",
        "S.1.GQ",
        "S.1.AD",
        "S.1.DP",
        "S.1.VF",
        "S.1.NL",
        "S.1.SB",
        "S.1.NC",
        "S.1.AQ",
        "S.1.GQX",
    ]

    cols = [
        "CHROM",
        "POS",
        "REF",
        "ALT",
        "FILTER",
        "GQX",
        "EVS",
        "T_DP",
        "T_DP_RATE",
        "T_AF",
        "tag",
    ]

    vcfheaders = list(extractHeaders(vcfname))

    evs_featurenames: Dict[str, int] = {}
    for l in vcfheaders:
        if "##snv_scoring_features" in l:
            try:
                xl = str(l).split("=", 1)
                xl = xl[1].split(",")
                for i, n in enumerate(xl):
                    evs_featurenames[n.strip()] = i
            except Exception:
                logging.warning(f"Failed to parse EVS feature annotation: {l}")

    result = vcfExtract(vcfname, features)

    if len(result) == 0:
        result = pandas.DataFrame(columns=cols)
        return result

    # Process Pisces features
    result["GQX"] = result["S.1.GQX"].astype(float)
    result["EVS"] = result["I.EVS"].astype(float)
    result["T_DP"] = result["S.1.DP"].astype(float)

    # Calculate depth rate
    if avg_depth is not None:
        result["T_DP_RATE"] = result["T_DP"] / avg_depth
    else:
        result["T_DP_RATE"] = 1.0

    # Extract allele frequency
    result["T_AF"] = 0.0
    try:
        # In Python 3, we need to convert to string first for string operations
        ad_strings = result["S.1.AD"].astype(str)
        for i, x in enumerate(ad_strings):
            try:
                xl = [int(a) for a in x.split(",")]
                if xl[0] + xl[1] > 0:
                    result.loc[i, "T_AF"] = float(xl[1]) / float(xl[0] + xl[1])
            except Exception:
                pass
    except Exception:
        logging.exception("Cannot extract AD.")

    result["tag"] = tag

    # Select and return only the columns we need
    return result[cols].copy()


def extractPiscesIndelFeatures(
    vcfname: str, tag: str, avg_depth: Optional[float] = None
) -> pandas.DataFrame:
    """Return a data frame with INDEL features collected from Pisces VCF.

    Args:
        vcfname: Name of the VCF file
        tag: Type of variants
        avg_depth: Average chromosome depths from BAM file

    Returns:
        DataFrame with extracted features
    """
    # Same as SNV features for Pisces
    return extractPiscesSNVFeatures(vcfname, tag, avg_depth)
