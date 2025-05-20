#!/usr/bin/env python3
"""
This file exists only for backward compatibility with build systems that don't
support PEP 517/518 yet. All relevant configuration is in pyproject.toml.

This file will be removed in version 1.0.0. Please use pip install . directly.
"""

import warnings

from setuptools import setup

warnings.warn(
    "setup.py is deprecated. Please use 'pip install .' directly."
    "This file will be removed in version 1.0.0.",
    DeprecationWarning,
    stacklevel=2,
)

setup()
