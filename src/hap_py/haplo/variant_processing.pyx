# cython: language_level=3
# distutils: language=c++

from __future__ import division, print_function

from libc.stdlib cimport free, malloc
from libc.string cimport memcpy
from libcpp.map cimport map
from libcpp.string cimport string
from libcpp.vector cimport vector

import numpy as np

cimport numpy as np

np.import_array()  # Initialize NumPy C API

# Helper functions for Python 3 string handling
cdef bytes py_str_to_bytes(s):
    """Convert Python string to bytes for C++ interop."""
    if isinstance(s, str):
        return s.encode('utf8')
    return s  # Assume already bytes

cdef string py_str_to_cpp_string(s):
    """Convert Python string to C++ string."""
    cdef bytes b = py_str_to_bytes(s)
    return string(<char*>b)

cdef str cpp_string_to_py_str(string& s):
    """Convert C++ string to Python string."""
    return s.c_str().decode('utf8')

# Safe VCF record processing with proper memory management
cdef class VariantRecord:
    cdef:
        string chrom
        long position
        string ref
        string alt
        vector[double] qual_values
        map[string, string] info_fields
        object _py_data  # Keep Python references

    def __cinit__(self):
        self._py_data = {}

    def __init__(self, chrom, position, ref, alt, qual=None, info=None):
        # Convert Python strings to C++ strings
        self.chrom = py_str_to_cpp_string(chrom)
        self.position = position
        self.ref = py_str_to_cpp_string(ref)
        self.alt = py_str_to_cpp_string(alt)

        # Store Python objects to prevent GC
        self._py_data['chrom'] = chrom
        self._py_data['ref'] = ref
        self._py_data['alt'] = alt

        # Handle quality values
        if qual is not None:
            if isinstance(qual, (list, tuple, np.ndarray)):
                for q in qual:
                    self.qual_values.push_back(float(q))
            else:
                self.qual_values.push_back(float(qual))

        # Handle info fields
        if info is not None:
            for key, value in info.items():
                self.info_fields[py_str_to_cpp_string(key)] = py_str_to_cpp_string(str(value))

    @property
    def chromosome(self):
        """Get chromosome name as Python string."""
        return cpp_string_to_py_str(self.chrom)

    @property
    def pos(self):
        """Get position."""
        return self.position

    @property
    def reference(self):
        """Get reference allele as Python string."""
        return cpp_string_to_py_str(self.ref)

    @property
    def alternate(self):
        """Get alternate allele as Python string."""
        return cpp_string_to_py_str(self.alt)

    def get_info(self, key):
        """Get info field value with proper Python 3 string handling."""
        cdef string cpp_key = py_str_to_cpp_string(key)
        if self.info_fields.find(cpp_key) != self.info_fields.end():
            return cpp_string_to_py_str(self.info_fields[cpp_key])
        return None

    def to_dict(self):
        """Convert to Python dictionary with proper string conversions."""
        result = {
            'chrom': self.chromosome,
            'pos': self.pos,
            'ref': self.reference,
            'alt': self.alternate,
            'qual': list(self.qual_values),
            'info': {cpp_string_to_py_str(k): cpp_string_to_py_str(v)
                     for k, v in self.info_fields}
        }
        return result

# Thread-safe processing with proper GIL management
def parallel_process_variants(variants, int threads=1):
    """Process variants in parallel with proper GIL handling."""
    cdef:
        vector[VariantRecord*] c_variants
        list py_variants = list(variants)  # Make sure it's a list
        VariantRecord* current_variant
        vector[VariantRecord*] results
        int i

    # Convert Python variants to C++ variants
    for var in py_variants:
        current_variant = new VariantRecord(
            var.chromosome, var.pos, var.reference, var.alternate
        )
        c_variants.push_back(current_variant)

    # Process variants in parallel without the GIL
    with nogil:
        results = _process_variants_parallel(c_variants, threads)

    # Convert results back to Python
    py_results = []
    for i in range(results.size()):
        py_results.append(results[i].to_dict())
        del results[i]  # Clean up

    # Clean up input variants
    for i in range(c_variants.size()):
        del c_variants[i]

    return py_results

cdef vector[VariantRecord*] _process_variants_parallel(vector[VariantRecord*] variants, int threads) nogil:
    """Process variants in parallel without the GIL."""
    # This would be your actual processing logic
    cdef:
        vector[VariantRecord*] results
        int i

    # For demonstration, just return the input variants
    for i in range(variants.size()):
        results.push_back(variants[i])

    return results
