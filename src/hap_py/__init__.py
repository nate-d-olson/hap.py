# Initialize hap_py package

# Import version information
try:
    from ._version import version as __version__
except ImportError:
    # Fallback for development installs without setuptools_scm
    __version__ = "0.4.0"
