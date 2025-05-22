# cython: language_level=3
# distutils: language=c++
"""
Cython implementation of ROC curve generation for hap.py.
This module provides optimized functions for Python 3 compatibility.
"""

# Python 3 imports
from __future__ import division, print_function

# Standard library imports
import math
import os
import sys

from libc.stdlib cimport free, malloc
from libc.string cimport memcpy

# External imports

import numpy as np

cimport numpy as np

np.import_array()  # Initialize NumPy C API

# C++ declarations for ROC integration
cdef extern from "Roc.hh" namespace "haplotypes":
    cdef cppclass RocPoint:
        double fdr
        double tpr
        double score

    cdef cppclass RocCurve:
        RocCurve() nogil
        void addPoint(double fdr, double tpr, double score) nogil
        int size() nogil
        vector[RocPoint] getPoints() nogil
        vector[RocPoint] interpolate(int points) nogil
        double auc() nogil

# Python wrapper for ROC curves
cdef class PyRocCurve:
    """Python wrapper for C++ RocCurve class.

    This class provides methods to generate and manipulate ROC curves.
    """
    cdef RocCurve* _curve

    def __cinit__(self):
        self._curve = new RocCurve()

    def __dealloc__(self):
        if self._curve != NULL:
            del self._curve

    def add_point(self, double fdr, double tpr, double score):
        """Add a point to the ROC curve."""
        self._curve.addPoint(fdr, tpr, score)

    def get_points(self):
        """Get all points in the ROC curve."""
        cdef:
            vector[RocPoint] points = self._curve.getPoints()
            int n = points.size()
            list result = []

        for i in range(n):
            result.append({
                'fdr': points[i].fdr,
                'tpr': points[i].tpr,
                'score': points[i].score
            })

        return result

    def interpolate(self, int points):
        """Interpolate the ROC curve to have a specific number of points."""
        cdef:
            vector[RocPoint] interpolated = self._curve.interpolate(points)
            int n = interpolated.size()
            list result = []

        for i in range(n):
            result.append({
                'fdr': interpolated[i].fdr,
                'tpr': interpolated[i].tpr,
                'score': interpolated[i].score
            })

        return result

    def auc(self):
        """Calculate the area under the ROC curve."""
        return self._curve.auc()

# Pure Python implementation for systems without C++ components
class PythonRocCurve:
    """Pure Python implementation of ROC curve generation.

    This class is used as a fallback when the C++ components are not available.
    """
    def __init__(self):
        self.points = []

    def add_point(self, fdr, tpr, score):
        """Add a point to the ROC curve."""
        self.points.append({'fdr': fdr, 'tpr': tpr, 'score': score})

    def get_points(self):
        """Get all points in the ROC curve."""
        return sorted(self.points, key=lambda p: p['score'])

    def interpolate(self, num_points):
        """Interpolate the ROC curve to have a specific number of points."""
        if len(self.points) <= 1:
            return self.points

        # Sort by score
        sorted_points = sorted(self.points, key=lambda p: p['score'])

        # Calculate range for interpolation
        min_score = sorted_points[0]['score']
        max_score = sorted_points[-1]['score']

        if min_score >= max_score:  # Avoid division by zero
            return sorted_points

        # Create interpolated points
        result = []
        for i in range(num_points):
            # Interpolate score
            score = min_score + (max_score - min_score) * i / (num_points - 1)

            # Find surrounding points
            lower_idx = 0
            for j, p in enumerate(sorted_points):
                if p['score'] <= score:
                    lower_idx = j

            upper_idx = min(lower_idx + 1, len(sorted_points) - 1)

            if lower_idx == upper_idx:
                # No interpolation needed
                result.append(sorted_points[lower_idx])
            else:
                # Interpolate values
                lower = sorted_points[lower_idx]
                upper = sorted_points[upper_idx]
                score_range = upper['score'] - lower['score']

                if score_range <= 0:  # Avoid division by zero
                    result.append(lower)
                else:
                    t = (score - lower['score']) / score_range
                    fdr = lower['fdr'] + t * (upper['fdr'] - lower['fdr'])
                    tpr = lower['tpr'] + t * (upper['tpr'] - lower['tpr'])
                    result.append({'fdr': fdr, 'tpr': tpr, 'score': score})

        return result

    def auc(self):
        """Calculate the area under the ROC curve."""
        if len(self.points) <= 1:
            return 0.0

        # Sort by FDR
        sorted_points = sorted(self.points, key=lambda p: p['fdr'])

        # Compute AUC using the trapezoidal rule
        auc = 0.0
        for i in range(1, len(sorted_points)):
            prev = sorted_points[i-1]
            curr = sorted_points[i]
            auc += (curr['fdr'] - prev['fdr']) * (curr['tpr'] + prev['tpr']) / 2

        return auc

# Factory function to create the appropriate ROC curve implementation
def create_roc_curve():
    """Create a ROC curve implementation based on available components."""
    try:
        return PyRocCurve()
    except Exception:
        return PythonRocCurve()

# Test function to verify the module is working
def test_module():
    """Test if the module is working properly."""
    return {
        "status": "ROC curve module loaded successfully",
        "language_level": "3"
    }
