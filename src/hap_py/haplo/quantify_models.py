"""
Data models for quantify operations.

This module defines the core data structures used by the quantify package,
ensuring type safety and compatibility with the original hap.py implementation.
"""

from dataclasses import dataclass
from enum import Enum
from typing import Dict, Optional, Union

import numpy as np


class VariantType(Enum):
    """Variant classification types."""

    SNP = "SNP"
    MNP = "MNP"
    INDEL = "INDEL"
    COMPLEX = "COMPLEX"
    UNKNOWN = "UNKNOWN"


class BenchmarkDecision(Enum):
    """Benchmark decision categories."""

    TP = "TP"  # True Positive
    FP = "FP"  # False Positive
    FN = "FN"  # False Negative
    UNK = "UNK"  # Unknown/Unclassified


class VariantMatchType(Enum):
    """Types of variant matches found by sophisticated matching algorithms."""

    EXACT = "exact"  # Exact coordinate and allele match
    OVERLAP = "overlap"  # Overlapping variants
    MULTIALLELIC = "multiallelic"  # Multi-allelic decomposition match
    SUPERLOCUS = "superlocus"  # Superlocus-level match
    NONE = "none"  # No match found


@dataclass
class QuantifyMetrics:
    """Core metrics calculated by quantify."""

    tp_count: int
    fp_count: int
    fn_count: int
    precision: float
    recall: float
    f1_score: float
    total_truth: int
    total_query: int

    def to_dict(self) -> Dict[str, Union[int, float]]:
        """Convert metrics to dictionary for output."""
        return {
            "TP": self.tp_count,
            "FP": self.fp_count,
            "FN": self.fn_count,
            "Precision": self.precision,
            "Recall": self.recall,
            "F1": self.f1_score,
            "Truth.total": self.total_truth,
            "Query.total": self.total_query,
        }


@dataclass
class ROCThresholds:
    """ROC curve threshold configuration."""

    field: str = "QUAL"
    min_threshold: float = 0.0
    max_threshold: float = 100.0
    delta: float = 0.5
    log_scale: bool = False

    def generate_thresholds(self) -> np.ndarray:
        """Generate array of quality thresholds for ROC analysis."""
        if self.log_scale:
            return np.logspace(
                np.log10(max(self.min_threshold, 0.001)),
                np.log10(self.max_threshold),
                int((self.max_threshold - self.min_threshold) / self.delta),
            )
        else:
            return np.arange(self.min_threshold, self.max_threshold, self.delta)


@dataclass
class StratificationRegion:
    """Stratification region definition."""

    name: str
    bed_file: Optional[str] = None
    filter_expression: Optional[str] = None

    def __post_init__(self):
        """Validate that either bed_file or filter_expression is provided."""
        if not self.bed_file and not self.filter_expression:
            raise ValueError("Either bed_file or filter_expression must be provided")
            raise ValueError("Either bed_file or filter_expression must be provided")
            raise ValueError("Either bed_file or filter_expression must be provided")
            raise ValueError("Either bed_file or filter_expression must be provided")
