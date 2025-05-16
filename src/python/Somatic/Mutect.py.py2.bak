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
Date:   2/10/2015
Author: Peter Krusche <pkrusche@illumina.com>
"""

import logging
import re

import pandas
from Tools.vcfextract import extractHeadersJSON, vcfExtract


def extractMutectSNVFeatures(vcfname, tag, avg_depth=None):
    """Return a data frame with features collected from the given VCF, tagged by given type"""
    records = []

    if not avg_depth:
        logging.warning(
            "No average depths available, normalized depth features cannot be calculated"
        )

    hdrs = extractHeadersJSON(vcfname)

    tsn = ""
    nsn = ""

    t_sample = "S.1."
    n_sample = "S.2."

    try:
        samples = hdrs["samples"]
        for f in hdrs["fields"]:
            if f["key"] == "GATKCommandLine" and f["values"]["ID"].lower() == "mutect":
                clopts = f["values"]["CommandLineOptions"]
                # ... tumor_sample_name=HCC2218_tumour ... normal_sample_name=HCC2218_normal
                m = re.search("tumor_sample_name=([^\s]+)", clopts)
                if m:
                    tsn = m.group(1)
                    for i, x in enumerate(samples):
                        if x == tsn:
                            t_sample = "S.%i." % (i + 1)
                            break
                m = re.search("normal_sample_name=([^\s]+)", clopts)
                if m:
                    nsn = m.group(1)
                    for i, x in enumerate(samples):
                        if x == nsn:
                            n_sample = "S.%i." % (i + 1)
                            break

    except Exception:
        logging.warning("Unable to detect tumour / normal sample order from VCF header")

    logging.info(
        "Normal sample name : %s (prefix %s) / tumour sample name : %s (prefix %s)"
        % (nsn, n_sample, tsn, t_sample)
    )

    features = [
        "CHROM",
        "POS",
        "REF",
        "ALT",
        "FILTER",
        "I.DB",
        "I.TLOD",
        "I.NLOD",
        "I.ECNT",
        "I.HCNT",
        "I.MAX_ED",
        "I.MIN_ED",
        n_sample + "GT",
        t_sample + "GT",
        n_sample + "DP",
        t_sample + "DP",
        n_sample + "QSS",
        t_sample + "QSS",
        n_sample + "AD",
        t_sample + "AD",
    ]

    has_warned = {}

    for vr in vcfExtract(vcfname, features):
        rec = {}
        for i, ff in enumerate(features):
            rec[ff] = vr[i]

        for q in [n_sample + "GT", t_sample + "GT"]:
            if q not in rec or rec[q] is None:
                rec[q] = "."
                if ("feat:" + q) not in has_warned:
                    logging.warning("Missing feature %s" % q)
                    has_warned["feat:" + q] = True

        # fix missing features
        for q in [
            "I.DB",
            "I.TLOD",
            "I.NLOD",
            "I.ECNT",
            "I.HCNT",
            "I.MAX_ED",
            "I.MIN_ED",
            n_sample + "GT",
            t_sample + "GT",
            n_sample + "DP",
            t_sample + "DP",
            n_sample + "QSS",
            t_sample + "QSS",
            n_sample + "AD",
            t_sample + "AD",
        ]:
            if q not in rec or rec[q] is None:
                rec[q] = 0
                if ("feat:" + q) not in has_warned:
                    logging.warning("Missing feature %s" % q)
                    has_warned["feat:" + q] = True
            else:
                # list features
                if q in [n_sample + "AD", t_sample + "AD"]:
                    if not isinstance(rec[q], list):
                        rec[q] = [rec[q], 0]
                    elif len(rec[q]) < 2:
                        rec[q] = rec[q] + [0] * (2 - len(rec[q]))

        rec["REF_COUNT_N"] = (
            rec[n_sample + "AD"][0] if isinstance(rec[n_sample + "AD"], list) else 0
        )
        rec["ALT_COUNT_N"] = (
            rec[n_sample + "AD"][1] if isinstance(rec[n_sample + "AD"], list) else 0
        )
        rec["REF_COUNT_T"] = (
            rec[t_sample + "AD"][0] if isinstance(rec[t_sample + "AD"], list) else 0
        )
        rec["ALT_COUNT_T"] = (
            rec[t_sample + "AD"][1] if isinstance(rec[t_sample + "AD"], list) else 0
        )

        # calculate VAFs
        depth_t = rec[t_sample + "DP"]
        ref_t = rec["REF_COUNT_T"]
        alt_t = rec["ALT_COUNT_T"]

        # normalise depth
        depth_t_norm = depth_t / float(avg_depth) if avg_depth else 0.0

        if depth_t <= 0:
            vaf_t = 0.0
        else:
            vaf_t = alt_t / float(depth_t)

        depth_n = rec[n_sample + "DP"]
        ref_n = rec["REF_COUNT_N"]
        alt_n = rec["ALT_COUNT_N"]

        depth_n_norm = depth_n / float(avg_depth) if avg_depth else 0.0

        if depth_n <= 0:
            vaf_n = 0.0
        else:
            vaf_n = alt_n / float(depth_n)

        rec["VAF_T"] = vaf_t
        rec["VAF_N"] = vaf_n

        rec["DEPTH_T"] = depth_t
        rec["DEPTH_T_NORM"] = depth_t_norm

        rec["DEPTH_N"] = depth_n
        rec["DEPTH_N_NORM"] = depth_n_norm

        rec["TAG"] = tag

        for k, v in list(rec.items()):
            if isinstance(v, list):
                rec[k] = ",".join(list(map(str, v)))

        records.append(rec)

    if not records:
        return None
    return pandas.DataFrame(records)


def extractMutectSNVFeatureVector(vcf, tag, avg_depth=None):
    """Get a feature vector for every call"""

    r = extractMutectSNVFeatures(vcf, tag, avg_depth)
    if r is None:
        return None

    vectors = []
    for _, vr in r.iterrows():
        rec = []
        rec.append(
            pandas.Series(
                [vr["CHROM"], vr["POS"], vr["REF"], vr["ALT"]],
                index=["CHROM", "POS", "REF", "ALT"],
            )
        )
        rec.append(pandas.Series([vr["TAG"]], index=["TAG"]))

        # basic features
        rec.append(pandas.Series([vr["FILTER"] == "PASS"], index=["PASS"]))
        if vr["I.DB"]:
            rec.append(pandas.Series([1], index=["IN_DB"]))
        else:
            rec.append(pandas.Series([0], index=["IN_DB"]))

        if vr["I.TLOD"] == ".":
            rec.append(pandas.Series([0], index=["TLOD"]))
        else:
            rec.append(pandas.Series([float(vr["I.TLOD"])], index=["TLOD"]))

        if vr["I.NLOD"] == ".":
            rec.append(pandas.Series([0], index=["NLOD"]))
        else:
            rec.append(pandas.Series([float(vr["I.NLOD"])], index=["NLOD"]))

        if vr["I.ECNT"] == ".":
            rec.append(pandas.Series([0], index=["ECNT"]))
        else:
            rec.append(pandas.Series([float(vr["I.ECNT"])], index=["ECNT"]))

        # string features, genotpyes
        gt_n = vr["S.2.GT"]
        f_name = "GT_N_" + str(gt_n)
        rec.append(pandas.Series([1], index=[f_name]))

        gt_t = vr["S.1.GT"]
        f_name = "GT_T_" + str(gt_t)
        rec.append(pandas.Series([1], index=[f_name]))

        # quality scores etc.

        rec.append(
            pandas.Series(
                [float(vr["S.2.QSS"]) if vr["S.2.QSS"] != "." else 0.0], index=["QSS_N"]
            )
        )
        rec.append(
            pandas.Series(
                [float(vr["S.1.QSS"]) if vr["S.1.QSS"] != "." else 0.0], index=["QSS_T"]
            )
        )

        rec.append(pandas.Series([vr["VAF_N"]], index=["VAF_N"]))
        rec.append(pandas.Series([vr["VAF_T"]], index=["VAF_T"]))

        rec.append(pandas.Series([vr["DEPTH_N"]], index=["DEPTH_N"]))
        rec.append(pandas.Series([vr["DEPTH_T"]], index=["DEPTH_T"]))

        rec.append(pandas.Series([vr["DEPTH_N_NORM"]], index=["DEPTH_N_NORM"]))
        rec.append(pandas.Series([vr["DEPTH_T_NORM"]], index=["DEPTH_T_NORM"]))

        rec.append(pandas.Series([vr["REF_COUNT_N"]], index=["REF_COUNT_N"]))
        rec.append(pandas.Series([vr["ALT_COUNT_N"]], index=["ALT_COUNT_N"]))
        rec.append(pandas.Series([vr["REF_COUNT_T"]], index=["REF_COUNT_T"]))
        rec.append(pandas.Series([vr["ALT_COUNT_T"]], index=["ALT_COUNT_T"]))

        if vr["VAF_T"] > vr["VAF_N"]:
            rec.append(pandas.Series([vr["VAF_T"] - vr["VAF_N"]], index=["VAF_DIFF"]))
        else:
            rec.append(pandas.Series([vr["VAF_N"] - vr["VAF_T"]], index=["VAF_DIFF"]))

        vectors.append(pandas.concat(rec))
    return pandas.DataFrame(vectors)


def extractMutectStandardFeatures(vcf, truth_regions=None, tag="UNK", avg_depth=None):
    """Get standard features for mutect SNV calls

    :param vcf: VCF file name with mutect calls
    :param truth_regions: BED file name with truth regions (or None)
    :param tag: tag for the calls
    :param avg_depth: Average depth
    :return:
    """
    return extractMutectSNVFeatures(vcf, tag, avg_depth)


def extractMutectFeatures(vcf, tag="UNK", avg_depth=None):
    """Get features for mutect SNV calls

    :param vcf: VCF file name with mutect calls
    :param tag: tag for the calls
    :param avg_depth: Average depth
    :return:
    """
    return extractMutectSNVFeatureVector(vcf, tag, avg_depth)
