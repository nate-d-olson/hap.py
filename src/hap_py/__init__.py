"""Initialization for the :mod:`hap_py` package.

This module exposes the package version as ``__version__`` and performs
basic checks for optional dependencies when hap.py is imported. The
docstring is used by the Sphinx documentation build.
"""

# Import version information
try:
    from ._version import version as __version__
except ImportError:
    # Fallback for development installs without setuptools_scm
    __version__ = "0.4.0"
