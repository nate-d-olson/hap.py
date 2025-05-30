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

"""
Haplo - Haplotype comparison module for hap.py

This package contains modules for comparing haplotypes and calculating
statistics for variant calling evaluation.
"""

# Haplo package initialization for Python 3
# Import commonly used functions directly from cython_compat

# Import Python standard libraries and compatibility code

# Import Cython modules with Python fallbacks
from .sequence_utils import complement_sequence, process_sequence, reverse_complement
from .variant_processor import VariantProcessor


# Compute ROC points in pure Python (recall, precision, F1 score, specificity)
def compute_roc_points(tp, fp, fn, tn=None, resolution=100):  # resolution unused
    """
    Compute ROC curve points from counts.

    Args:
        tp: iterable of true positives counts
        fp: iterable of false positives counts
        fn: iterable of false negatives counts
        tn: iterable of true negatives counts (optional)
        resolution: ignored

    Returns:
        tuple of lists: (recalls, precisions, f1_scores, specificities)
    """
    tps = list(tp)
    fps = list(fp)
    fns = list(fn)
    recalls = []
    precisions = []
    for t, f_p, f_n in zip(tps, fps, fns):
        denom_r = t + f_n
        recalls.append(t / denom_r if denom_r else 0.0)
        denom_p = t + f_p
        precisions.append(t / denom_p if denom_p else 0.0)
    f1_scores = []
    for r, p in zip(recalls, precisions):
        denom = p + r
        f1_scores.append(2 * p * r / denom if denom else 0.0)
    if tn is not None:
        tns = list(tn)
        specificities = []
        for t_n, f_p in zip(tns, fps):
            denom_s = t_n + f_p
            specificities.append(t_n / denom_s if denom_s else 0.0)
    else:
        specificities = [0.0] * len(recalls)
    return recalls, precisions, f1_scores, specificities


# Chromosome comparison and sorting utilities
def cmp_chromosomes(a, b):
    """Compare two chromosomes for sorting."""
    a_chr = a[3:] if a.startswith("chr") else a
    b_chr = b[3:] if b.startswith("chr") else b
    if a_chr.isdigit() and b_chr.isdigit():
        return int(a_chr) - int(b_chr)
    elif a_chr.isdigit():
        return -1
    elif b_chr.isdigit():
        return 1
    else:
        special_order = {"X": 1, "Y": 2, "M": 3, "MT": 3}
        a_val = special_order.get(a_chr, 99)
        b_val = special_order.get(b_chr, 99)
        if a_val != 99 and b_val != 99:
            return a_val - b_val
        elif a_val != 99:
            return -1
        elif b_val != 99:
            return 1
        else:
            return -1 if a < b else (1 if a > b else 0)


def sort_chromosomes(chrom_list):
    """Sort a list of chromosomes in standard order."""
    from functools import cmp_to_key

    return sorted(chrom_list, key=cmp_to_key(cmp_chromosomes))


# Version information
__version__ = "0.4.0"


# Info about package
def get_module_info():
    """Get information about the module implementation."""
    modules = {}
    for name in [
        "complement_sequence",
        "reverse_complement",
        "process_sequence",
        "VariantProcessor",
        "compute_roc_points",
        "cmp_chromosomes",
        "sort_chromosomes",
    ]:
        obj = globals().get(name)
        modules[name] = {"module": obj.__module__, "is_cython": False}
    return {"version": __version__, "modules": modules}


# Package exports
__all__ = [
    "complement_sequence",
    "reverse_complement",
    "process_sequence",
    "VariantProcessor",
    "compute_roc_points",
    "cmp_chromosomes",
    "sort_chromosomes",
    "get_module_info",
]
