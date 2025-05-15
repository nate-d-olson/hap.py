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
Haplo - Haplotype comparison module for hap.py

This package contains modules for comparing haplotypes and calculating
statistics for variant calling evaluation.
"""

__version__ = "0.3.15"
has_sge = False

import logging
import os

# Check if we should use the mock implementation (useful for testing/development)
use_mock = os.environ.get("HAPLO_USE_MOCK", "0").lower() in ("1", "true", "yes", "on")

if use_mock:
    # Explicitly use the mock implementation
    try:
        from Haplo.cython.mock_internal import *  # noqa

        logging.info("Using mock C++ implementation for Haplo module")
    except ImportError:
        logging.warning(
            "Could not import mock Haplo implementation. Functionality will be limited."
        )
else:
    # Try to use the real C++ implementation, fall back to mock if needed
    try:
        # Try to import the C++ implementation
        from Haplo.cython._internal import *  # noqa
    except ImportError:
        try:
            # Fall back to mock implementation
            from Haplo.cython.mock_internal import *  # noqa

            logging.warning(
                "Could not import Haplo C++ extension module. Using mock implementation instead."
            )
        except ImportError:
            logging.warning(
                "Could not import Haplo C++ extension module or mock implementation. Functionality will be limited."
            )
