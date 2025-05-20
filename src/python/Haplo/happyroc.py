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

import logging
import re

import numpy as np
import pandas
from Tools import ci

RESULT_ALLCOLUMNS = [
    "Type",
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
    "Subset.Level",
]

RESULT_ALLDTYPES = [
    str,
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
    int,
]

for count_type in [
    "TRUTH.TOTAL",
    "TRUTH.TP",
    "TRUTH.FN",
    "QUERY.TOTAL",
    "QUERY.TP",
    "QUERY.FP",
    "QUERY.UNK",
]:
    RESULT_ALLCOLUMNS.append(count_type)
    RESULT_ALLCOLUMNS.append(count_type + ".ti")
    RESULT_ALLCOLUMNS.append(count_type + ".tv")
    RESULT_ALLCOLUMNS.append(count_type + ".het")
    RESULT_ALLCOLUMNS.append(count_type + ".homalt")
    RESULT_ALLCOLUMNS.append(count_type + ".TiTv_ratio")
    RESULT_ALLCOLUMNS.append(count_type + ".het_hom_ratio")
    RESULT_ALLDTYPES += [float] * 7


def roc(
    roc_table, output_path, filter_handling=None, ci_alpha=0.05, total_region_size=None
):
    """Calculate SNP and indel ROC.

    Return a dictionary of variant types and corresponding files.

    :param filter_handling: can be None, "PASS" or "ALL" to filter rows based on the "Filter" column.
                            this is necessary because vcfeval doesn't preserve filter information
                            in GA4GH output mode, so we need to remove the corresponding rows here
    :param ci_alpha: Jeffrey's CI confidence level for recall, precision, na
    :param total_region_size: correct Subset.Size for "*" region if a subset was selected in hap.py

    """
    result = {}
    header = None
    with open(roc_table, encoding="utf-8") as rt:
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
                    except Exception:
                        pass

                try:
                    if (
                        rec["Type"] in ["SNP", "INDEL"]
                        and rec["Filter"] == "ALL"
                        and rec["Subset"] == "*"
                        and rec["Genotype"] == "*"
                        and rec["Subtype"] == "*"
                        and rec["QQ"] != "*"
                    ):  # this is the ROC score field
                        roc = "Locations." + rec["Type"]
                        if roc not in result:
                            result[roc] = [rec]
                        else:
                            result[roc].append(rec)
                except Exception:
                    pass

                try:
                    if (
                        rec["Type"] in ["SNP", "INDEL"]
                        and rec["Filter"] == "PASS"
                        and rec["Subset"] == "*"
                        and rec["Genotype"] == "*"
                        and rec["Subtype"] == "*"
                        and rec["QQ"] != "*"
                    ):  # this is the ROC score field
                        roc = "Locations." + rec["Type"] + ".PASS"
                        if roc not in result:
                            result[roc] = [rec]
                        else:
                            result[roc].append(rec)
                except Exception:
                    pass

                try:
                    if (
                        rec["Type"] in ["SNP", "INDEL"]
                        and rec["Filter"] == "SEL"
                        and rec["Subset"] == "*"
                        and rec["Genotype"] == "*"
                        and rec["Subtype"] == "*"
                        and rec["QQ"] != "*"
                    ):  # this is the ROC score field
                        roc = "Locations." + rec["Type"] + ".SEL"
                        if roc not in result:
                            result[roc] = [rec]
                        else:
                            result[roc].append(rec)
                except Exception:
                    pass

                roc = "all"
                if roc not in result:
                    result[roc] = [rec]
                else:
                    result[roc].append(rec)

    if "all" not in result:
        # minimal empty DF
        minidata = [
            {
                "Type": "SNP",
                "Subtype": "*",
                "Filter": "ALL",
                "Genotype": "*",
                "Subset": "*",
                "QQ": "*",
            }
            for _ in range(2)
        ]
        minidata[1]["Type"] = "INDEL"
        result["all"] = pandas.DataFrame(minidata, columns=RESULT_ALLCOLUMNS)
        for i, c in enumerate(RESULT_ALLCOLUMNS):
            result["all"][c] = result["all"][c].astype(
                RESULT_ALLDTYPES[i], raise_on_error=False
            )

    for k, v in list(result.items()):
        result[k] = _postprocessRocData(pandas.DataFrame(v, columns=RESULT_ALLCOLUMNS))

        # compute ratios
        for count_type in [
            "TRUTH.TOTAL",
            "TRUTH.FN",
            "TRUTH.TP",
            "QUERY.FP",
            "QUERY.TP",
            "QUERY.TOTAL",
            "QUERY.UNK",
        ]:
            result[k][count_type + ".TiTv_ratio"] = pandas.to_numeric(
                result[k][count_type + ".ti"], errors="coerce"
            ) / pandas.to_numeric(result[k][count_type + ".tv"], errors="coerce")
            result[k][count_type + ".het_hom_ratio"] = pandas.to_numeric(
                result[k][count_type + ".het"], errors="coerce"
            ) / pandas.to_numeric(result[k][count_type + ".homalt"], errors="coerce")
            result[k][count_type + ".TiTv_ratio"].replace(
                [np.inf, -np.inf], np.nan, inplace=True
            )
            result[k][count_type + ".het_hom_ratio"].replace(
                [np.inf, -np.inf], np.nan, inplace=True
            )

        if 0 < ci_alpha < 1:
            logging.info("Computing recall CIs for %s" % k)
            rc, rc_min, rc_max = ci.binomialCI(
                result[k]["TRUTH.TP"].values,
                (result[k]["TRUTH.TP"] + result[k]["TRUTH.FN"]).values,
                ci_alpha,
            )
            result[k]["METRIC.Recall.Lower"] = rc_min
            result[k]["METRIC.Recall.Upper"] = rc_max

            logging.info("Computing precision CIs for %s" % k)
            pc, pc_min, pc_max = ci.binomialCI(
                result[k]["QUERY.TP"].values,
                (result[k]["QUERY.TP"] + result[k]["QUERY.FP"]).values,
                ci_alpha,
            )
            result[k]["METRIC.Precision.Lower"] = pc_min
            result[k]["METRIC.Precision.Upper"] = pc_max

            logging.info("Computing Frac_NA CIs for %s" % k)
            fna, fna_min, fna_max = ci.binomialCI(
                result[k]["QUERY.UNK"].values, result[k]["QUERY.TOTAL"].values, ci_alpha
            )
            result[k]["METRIC.Frac_NA.Lower"] = fna_min
            result[k]["METRIC.Frac_NA.Upper"] = fna_max

        # write correct subset.size
        if total_region_size is not None:
            result[k].loc[result[k]["Subset"] == "*", "Subset.Size"] = total_region_size

        vt = re.sub("[^A-Za-z0-9\\.\\-_]", "_", k, flags=re.IGNORECASE)
        if output_path:
            result[k].to_csv(
                output_path + "." + vt + ".csv.gz", index=False, compression="gzip"
            )

    return result


def _postprocessRocData(roctable):
    """post-process ROC data by correcting the types,
    sorting and changing the ordering of the columns
    """
    if roctable.empty:
        return roctable

    def safe_int(f):
        try:
            return int(f)
        except Exception:
            return 0

    typeconv = [
        ("METRIC.Recall", None),
        ("METRIC.Precision", None),
        ("METRIC.Frac_NA", None),
        ("TRUTH.TP", safe_int),
        ("TRUTH.FN", safe_int),
        ("QUERY.TP", safe_int),
        ("QUERY.FP", safe_int),
        ("QUERY.UNK", safe_int),
        ("QUERY.TOTAL", safe_int),
        ("TRUTH.TOTAL", safe_int),
        ("FP.al", safe_int),
        ("FP.gt", safe_int),
        ("TRUTH.TOTAL.TiTv_ratio", None),
        ("TRUTH.TOTAL.het_hom_ratio", None),
        ("TRUTH.FN.TiTv_ratio", None),
        ("TRUTH.FN.het_hom_ratio", None),
        ("TRUTH.TP.TiTv_ratio", None),
        ("TRUTH.TP.het_hom_ratio", None),
        ("METRIC.F1_Score", None),
        ("QUERY.FP.TiTv_ratio", None),
        ("QUERY.FP.het_hom_ratio", None),
        ("QUERY.TP.TiTv_ratio", None),
        ("QUERY.TOTAL.TiTv_ratio", None),
        ("QUERY.TOTAL.het_hom_ratio", None),
        ("QUERY.TP.het_hom_ratio", None),
        ("QUERY.UNK.TiTv_ratio", None),
        ("QUERY.UNK.het_hom_ratio", None),
    ]

    for col, c in typeconv:
        roctable[col] = pandas.to_numeric(roctable[col], errors="coerce")
        if c:
            roctable[col] = roctable[col].apply(c)

    roctable.sort_values(
        ["Type", "Subtype", "Subset", "Filter", "Genotype", "QQ.Field", "QQ"],
        inplace=True,
    )

    return roctable[RESULT_ALLCOLUMNS]


# ROC curve support: pure-Python implementation (no C++ dependencies)
class PythonRocCurve:
    """Pure Python implementation of ROC curve generation."""

    def __init__(self):
        self.points = []

    def add_point(self, fdr, tpr, score):
        """Add a point to the ROC curve."""
        self.points.append({"fdr": fdr, "tpr": tpr, "score": score})

    def get_points(self):
        """Get all points in the ROC curve."""
        return sorted(self.points, key=lambda p: p["score"])

    def interpolate(self, num_points):
        """Interpolate the ROC curve to have a specific number of points."""
        if len(self.points) <= 1:
            return self.points

        sorted_points = sorted(self.points, key=lambda p: p["score"])
        min_score = sorted_points[0]["score"]
        max_score = sorted_points[-1]["score"]

        if min_score >= max_score:
            return sorted_points

        result = []
        for i in range(num_points):
            score = min_score + (max_score - min_score) * i / (num_points - 1)
            lower_idx = 0
            for j, p in enumerate(sorted_points):
                if p["score"] <= score:
                    lower_idx = j

            upper_idx = min(lower_idx + 1, len(sorted_points) - 1)
            if lower_idx == upper_idx:
                result.append(sorted_points[lower_idx])
            else:
                lower = sorted_points[lower_idx]
                upper = sorted_points[upper_idx]
                diff = upper["score"] - lower["score"]
                if diff <= 0:
                    result.append(lower)
                else:
                    t = (score - lower["score"]) / diff
                    fdr = lower["fdr"] + t * (upper["fdr"] - lower["fdr"])
                    tpr = lower["tpr"] + t * (upper["tpr"] - lower["tpr"])
                    result.append({"fdr": fdr, "tpr": tpr, "score": score})

        return result

    def auc(self):
        """Calculate the area under the ROC curve."""
        if len(self.points) <= 1:
            return 0.0

        sorted_points = sorted(self.points, key=lambda p: p["fdr"])
        auc = 0.0
        for i in range(1, len(sorted_points)):
            prev = sorted_points[i - 1]
            curr = sorted_points[i]
            auc += (curr["fdr"] - prev["fdr"]) * (curr["tpr"] + prev["tpr"]) / 2
        return auc


# Factory function to create the ROC curve implementation
def create_roc_curve():
    """Create a ROC curve implementation (pure Python)."""
    return PythonRocCurve()


# Test function to verify module import
def test_module():
    """Test if the ROC module is loaded correctly."""
    return {"status": "ROC curve module loaded successfully", "language_level": "3"}
