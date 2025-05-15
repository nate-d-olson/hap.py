"""
Haplo package initialization - Python 3 compatible version

This module provides initialization for the Haplo package, ensuring
proper imports of all submodules and handling Python 3 compatibility.
"""

import logging
import os
import sys

# Import version information
__version__ = "0.3.12-py3"

# Set up path for submodules
current_dir = os.path.dirname(os.path.abspath(__file__))
if current_dir not in sys.path:
    sys.path.insert(0, current_dir)

# Configure fallback for Cython modules
try:
    from .cython import USING_MOCK

    if USING_MOCK:
        logging.warning("Using mock implementations for Cython modules")
except ImportError:
    logging.warning("Failed to import Cython modules")
