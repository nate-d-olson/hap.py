"""
Python 3 compatible Cython module package with fallback mechanism

This package tries to import the Cython modules, but falls back to
mock implementations if they are not available.
"""

import logging
import os
import warnings

# Flag to control whether to use mock implementation
USE_MOCK = os.environ.get("HAPLO_USE_MOCK", "0").lower() in ("1", "true", "yes")

try:
    if USE_MOCK:
        # Force using mock implementation
        raise ImportError("Using mock implementation as requested")

    # Try to import the Cython modules
    from ._internal import *
    from .cpp_internal import *

    # Add any other Cython modules here
    # Signal that we're using the real implementation
    USING_MOCK = False

except ImportError as e:
    warnings.warn(f"Failed to import Cython modules: {e}. Using mock implementation.", stacklevel=2)
    logging.warning(f"Failed to import Cython modules: {e}. Using mock implementation.")

    # Fall back to mock implementation
    from .mock_internal import *

    # Add other mock imports here
    # Signal that we're using the mock implementation
    USING_MOCK = True

# Version information
__version__ = "0.1.0"
