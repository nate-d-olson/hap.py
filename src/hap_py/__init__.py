"""hap.py - Variant calling comparison and benchmarking."""

import logging
import os
import sys
from pathlib import Path

try:
    from ._version import version as __version__
except ImportError:
    # Fallback for development installs
    __version__ = "0.4.0"


def init():
    """Initialize the hap.py package, checking for required tools."""
    # Import here to avoid circular imports
    from .external.rtg_manager import rtg_manager

    # Check for RTG tools - this will download if needed
    try:
        rtg_path = rtg_manager.get_rtg_path()
        logging.info(f"Using RTG Tools at: {rtg_path}")
    except Exception as e:
        logging.warning(f"RTG Tools setup issue: {str(e)}")
        logging.warning(
            "Some functionality may be limited. Set RTG_PATH manually if needed."
        )


# Run initialization when importing the package
init()
