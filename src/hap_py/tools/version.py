#!/usr/bin/env python3

"""Version information for hap.py"""

# Import version from the main package
try:
    from .._version import version
except ImportError:
    # Fallback version if _version.py is not available
    version = "0.4.0"

# Feature flags
has_vcfeval = True  # RTG tools are now available
