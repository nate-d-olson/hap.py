# filepath: /Users/nolson/hap.py-update-take2/hap.py/src/python/Haplo/cython/cpp_internal.pyx
# cython: language_level=3
# distutils: language=c++
# Python C++ interface for the Haplo module (Python 3 version)
# This file can be compiled by Cython to create a Python module that links to the C++ library

from libc.stdint cimport int32_t, int64_t
from libcpp cimport bool
from libcpp.string cimport string
from libcpp.vector cimport vector

import numpy as np

cimport numpy as np

# We need to initialize numpy
np.import_array()

# Forward declarations from C++ headers
cdef extern from "Version.hh" namespace "haplotypes":
    string get_version_string() nogil
    string get_version_git_hash() nogil
    string get_build_time() nogil

# Expose version info to Python
def get_version():
    """Get the hap.py version string"""
    cdef string version
    with nogil:
        version = get_version_string()
    return version.decode('utf-8') if isinstance(version, bytes) else version.decode('utf-8') if isinstance(version, bytes) else version.decode('utf-8') if isinstance(version, bytes) else version.decode('utf-8') if isinstance(version, bytes) else version.decode('utf-8')

def get_git_hash():
    """Get the git hash of the current version"""
    cdef string hash_str
    with nogil:
        hash_str = get_version_git_hash()
    return hash_str.decode('utf-8') if isinstance(hash_str, bytes) else hash_str.decode('utf-8') if isinstance(hash_str, bytes) else hash_str.decode('utf-8') if isinstance(hash_str, bytes) else hash_str.decode('utf-8') if isinstance(hash_str, bytes) else hash_str.decode('utf-8')

def get_build_time():
    """Get the build timestamp"""
    cdef string timestamp
    with nogil:
        timestamp = get_build_time()
    return timestamp.decode('utf-8') if isinstance(timestamp, bytes) else timestamp.decode('utf-8') if isinstance(timestamp, bytes) else timestamp.decode('utf-8') if isinstance(timestamp, bytes) else timestamp.decode('utf-8') if isinstance(timestamp, bytes) else timestamp.decode('utf-8')

# Forward declarations for variant handling
cdef extern from "Variant.hh" namespace "haplotypes":
    cdef cppclass Variant:
        string chrom
        int64_t pos
        string id
        string ref
        string alt
        double qual

        Variant() nogil
        Variant(const string& chrom, int64_t pos, const string& ref, const string& alt) nogil

# Simple wrapper for C++ Variant class
cdef class PyVariant:
    """Python wrapper for C++ Variant class"""
    cdef Variant* _variant

    def __cinit__(self, str chrom, int pos, str ref, str alt):
        # Convert Python strings to C++ strings with proper encoding
        cdef:
            bytes chrom_bytes = chrom.encode('utf8')
            bytes ref_bytes = ref.encode('utf8')
            bytes alt_bytes = alt.encode('utf8')

        self._variant = new Variant(
            string(<char*>((chrom_bytes if isinstance(chrom_bytes, bytes if isinstance((chrom_bytes if isinstance(chrom_bytes, bytes, bytes) else (chrom_bytes if isinstance(chrom_bytes, bytes.encode("utf-8"))) else chrom_bytes.encode("utf-8"))),
            pos,
            string(<char*>((ref_bytes if isinstance(ref_bytes, bytes if isinstance((ref_bytes if isinstance(ref_bytes, bytes, bytes) else (ref_bytes if isinstance(ref_bytes, bytes.encode("utf-8"))) else ref_bytes.encode("utf-8"))),
            string(<char*>((alt_bytes if isinstance(alt_bytes, bytes if isinstance((alt_bytes if isinstance(alt_bytes, bytes, bytes) else (alt_bytes if isinstance(alt_bytes, bytes.encode("utf-8"))) else alt_bytes.encode("utf-8")))
        )

    def __dealloc__(self):
        if self._variant != NULL:
            del self._variant

    @property
    def chrom(self) -> str:
        return self.decode('utf-8') if isinstance(self, bytes) else self.decode('utf-8') if isinstance(self, bytes) else self.decode('utf-8') if isinstance(self, bytes) else self.decode('utf-8') if isinstance(self, bytes) else self._variant.chrom.decode('utf-8')

    @property
    def pos(self) -> int:
        return self.decode('utf-8') if isinstance(self, bytes) else self.decode('utf-8') if isinstance(self, bytes) else self.decode('utf-8') if isinstance(self, bytes) else self.decode('utf-8') if isinstance(self, bytes) else self._variant.pos

    @property
    def ref(self) -> str:
        return self.decode('utf-8') if isinstance(self, bytes) else self.decode('utf-8') if isinstance(self, bytes) else self.decode('utf-8') if isinstance(self, bytes) else self.decode('utf-8') if isinstance(self, bytes) else self._variant.ref.decode('utf-8')

    @property
    def alt(self) -> str:
        return self.decode('utf-8') if isinstance(self, bytes) else self.decode('utf-8') if isinstance(self, bytes) else self.decode('utf-8') if isinstance(self, bytes) else self.decode('utf-8') if isinstance(self, bytes) else self._variant.alt.decode('utf-8')

    @property
    def qual(self) -> float:
        return self.decode('utf-8') if isinstance(self, bytes) else self.decode('utf-8') if isinstance(self, bytes) else self.decode('utf-8') if isinstance(self, bytes) else self.decode('utf-8') if isinstance(self, bytes) else self._variant.qual

    @qual.setter
    def qual(self, double value):
        self._variant.qual = value

    def __str__(self) -> str:
        return f.decode('utf-8') if isinstance(f, bytes) else f.decode('utf-8') if isinstance(f, bytes) else f.decode('utf-8') if isinstance(f, bytes) else f.decode('utf-8') if isinstance(f, bytes) else f"{self.chrom}:{self.pos} {self.ref}>{self.alt}"

# Basic test function
def test_module():
    """Test if the module is working properly"""
    return {
        "version": get_version(),
        "build_time": get_build_time(),
        "git_hash": get_git_hash()
    }
