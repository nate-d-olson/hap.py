# cython: language_level=3
# distutils: language=c++
"""
This module contains Cython implementations of the variant processor pipeline
for Python 3 compatibility.
"""

from libc.stdint cimport int32_t, int64_t
from libcpp cimport bool
from libcpp.map cimport map
from libcpp.string cimport string
from libcpp.vector cimport vector

import numpy as np

cimport numpy as np


# Forward declarations from C++ headers
cdef extern from "VariantProcessor.hh" namespace "haplotypes":
    cdef cppclass VariantProcessor:
        VariantProcessor() nogil
        bool processVariant(string& chrom, int64_t& pos, string& ref, string& alt) nogil
        void setLeftNormalize(bool normalize) nogil
        void setSplitAlleles(bool split) nogil
        void setRemovePadding(bool remove) nogil

# Pipeline components - forward declarations
cdef extern from "VariantAlleleNormalizer.hh" namespace "haplotypes":
    cdef cppclass VariantAlleleNormalizer:
        VariantAlleleNormalizer(const string& reference_fasta) nogil
        bool normalize(string& chrom, int64_t& pos, string& ref, string& alt) nogil

cdef extern from "VariantAlleleSplitter.hh" namespace "haplotypes":
    cdef cppclass VariantAlleleSplitter:
        VariantAlleleSplitter() nogil
        bool split(string& chrom, int64_t& pos, string& ref, vector[string]& alts) nogil

cdef extern from "VariantLeftPadding.hh" namespace "haplotypes":
    cdef cppclass VariantLeftPadding:
        VariantLeftPadding() nogil
        bool removePadding(string& ref, string& alt) nogil

# Python wrapper for variant processor
cdef class PyVariantProcessor:
    """Python wrapper for C++ VariantProcessor pipeline.

    This class provides methods to process variants through a standard
    pipeline of operations:
    - Normalize variants (left-shift)
    - Split multi-allelic variants
    - Remove redundant padding bases
    """
    cdef VariantProcessor* _processor

    def __cinit__(self):
        self._processor = new VariantProcessor()

    def __dealloc__(self):
        if self._processor != NULL:
            del self._processor

    def set_left_normalize(self, bool normalize):
        """Enable/disable left normalization of variants."""
        self._processor.setLeftNormalize(normalize)

    def set_split_alleles(self, bool split):
        """Enable/disable splitting of multi-allelic variants."""
        self._processor.setSplitAlleles(split)

    def set_remove_padding(self, bool remove):
        """Enable/disable removal of redundant padding bases."""
        self._processor.setRemovePadding(remove)

    def process_variant(self, str chrom, int pos, str ref, str alt):
        """Process a variant through the pipeline.

        Args:
            chrom: Chromosome name
            pos: 1-based position
            ref: Reference allele
            alt: Alternate allele

        Returns:
            Tuple of (chrom, pos, ref, alt) with processed values,
            or None if processing failed
        """
        # Convert Python strings to C++
        cdef:
            string c_chrom = string(chrom.encode('utf8'))
            int64_t c_pos = pos
            string c_ref = string(ref.encode('utf8'))
            string c_alt = string(alt.encode('utf8'))
            bool success

        # Call C++ processing function
        with nogil:
            success = self._processor.processVariant(c_chrom, c_pos, c_ref, c_alt)

        if not success:
            return None

        # Convert back to Python types
        return (c_chrom.decode('utf8'), c_pos, c_ref.decode('utf8'), c_alt.decode('utf8'))

# Function to create a variant processor with standard settings
def create_standard_processor():
    """Create a variant processor with standard settings."""
    processor = PyVariantProcessor()
    processor.set_left_normalize(True)
    processor.set_split_alleles(False)  # Split handled separately in pipeline
    processor.set_remove_padding(True)
    return processor

# Simple Python function to test the module
def test_module():
    """Test if the module is working properly"""
    return {
        "status": "Variant processor module loaded successfully",
        "language_level": "3"
    }
