# filepath: /Users/nolson/hap.py-update-take2/hap.py/src/python/Haplo/cython/_internal.pyx
# Internal C++ wrapper module for hap.py
# For Python 3 integration of the C++ components

# distutils: language=c++
# cython: language_level=3
# cython: boundscheck=False
# cython: wraparound=False

"""
This module provides the interface between Python and the C++ components of hap.py.
It wraps the C++ library functionality for use in the Python code.
"""

from libc.stdint cimport int32_t, int64_t
from libcpp cimport bool
from libcpp.string cimport string
from libcpp.vector cimport vector

import numpy as np

cimport numpy as np

import logging
import os

# Import Python builtins
from builtins import bytes, str

# Initialize NumPy (required for using numpy C API)
np.import_array()

# Forward declarations of C++ types
cdef extern from "Version.hh" namespace "haplotypes":
    string version_string() nogil
    string build_timestamp() nogil

# Attempt to import the C++ functionality - this is defined in a separate file
# that's generated during the build process
try:
    # Define Python-accessible version information
    def get_version():
        """Get the hap.py version string"""
        cdef string v = version_string()
        return v.decode('utf-8')

    def get_build_time():
        """Get the hap.py build timestamp"""
        cdef string t = build_timestamp()
        return t.decode('utf-8')

except Exception as e:
    logging.warning("Could not initialize C++ components: {}".format(str(e)))

    # Define fallback versions of the functions
    def get_version():
        """Fallback version function when C++ is not available"""
        return "unknown (C++ module not loaded)"

    def get_build_time():
        """Fallback build time function when C++ is not available"""
        return "unknown (C++ module not loaded)"

# Add module initialization function
def is_available():
    """Check if the C++ components are available"""
    try:
        v = get_version()
        return v != "unknown (C++ module not loaded)"
    except:
        return False

    # Add placeholder functions for any other functionality needed

# Create a test function to verify the module is working
def test_module():
    """Test if the module is working properly"""
    return {"version": get_version(), "build_time": get_build_time()}
