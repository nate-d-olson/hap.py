# filepath: /Users/nolson/hap.py-update-take2/hap.py/src/python/Haplo/cython/_internal.pyx
# Internal C++ wrapper module for hap.py
# For Python 3 integration of the C++ components

# cython: language_level=3
# distutils: language=c++
# cython: boundscheck=False
# cython: wraparound=False

"""
This module provides the interface between Python and the C++ components of hap.py.
It wraps the C++ library functionality for use in the Python code.
"""

from libcpp.string cimport string
from libcpp.vector cimport vector
from libc.stdint cimport int32_t, int64_t
from libcpp cimport bool
import numpy as np
cimport numpy as np
import os
import logging

# Import Python builtins
from builtins import bytes, str

# Initialize NumPy (required for using numpy C API)
np.import_array()

# Forward declarations of C++ types
cdef extern from "Version.hh" namespace "haplotypes":
    cdef string version_string()
    cdef string build_timestamp()

# Attempt to import the C++ functionality - this is defined in a separate file
# that's generated during the build process
try:
    # Define Python-accessible version information
    def get_version():
        """Get the hap.py version string"""
        return version_string.decode('utf-8') if isinstance(version_string, bytes) else version_string().decode('utf-8')
    
    def get_build_time():
        """Get the hap.py build timestamp"""
        return build_timestamp.decode('utf-8') if isinstance(build_timestamp, bytes) else build_timestamp().decode('utf-8')
        
except Exception as e:
    logging.warning(f"Could not initialize C++ components: {e}")
    
    # Define fallback versions of the functions
    def get_version():
        """Get the hap.py version string (fallback)"""
        return "0.0.0-fallback"
    
    def get_build_time():
        """Get the hap.py build timestamp (fallback)"""
        return "unknown"
    
    # Add placeholder functions for any other functionality needed

# Create a test function to verify the module is working
def test_module():
    """Test if the module is working properly"""
    return {"version": get_version(), "build_time": get_build_time()}
