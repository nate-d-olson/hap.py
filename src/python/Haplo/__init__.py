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
from .cython_compat import (
    VariantProcessor,
    cmp_chromosomes,
    complement_sequence,
    compute_roc_points,
    is_using_cython,
    reverse_complement,
    sort_chromosomes,
)

# Version information
__version__ = "0.4.0"

# Check if we're using Cython or Python implementations
_module_info = is_using_cython()


# Info about package
def get_module_info():
    """Get information about the module implementation."""
    result = {"version": __version__, "modules": _module_info}
    return result


# Package exports
__all__ = [
    "complement_sequence",
    "reverse_complement",
    "VariantProcessor",
    "compute_roc_points",
    "cmp_chromosomes",
    "sort_chromosomes",
    "get_module_info",
]
