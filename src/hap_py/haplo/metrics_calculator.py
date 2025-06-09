"""
Metrics calculation for variant benchmarking.

This module provides robust calculation of standard benchmarking metrics
with proper handling of edge cases and numerical stability.
"""

import logging
from typing import Dict, Optional

import numpy as np
import pandas as pd

from .quantify_models import QuantifyMetrics

logger = logging.getLogger(__name__)


class MetricsCalculator:
    """
    Calculator for variant benchmarking metrics.

    Provides methods to calculate precision, recall, F1-score, and other
    standard metrics used in variant calling benchmarking.
    """

    @staticmethod
    def calculate_basic_metrics(
        tp_count: int,
        fp_count: int,
        fn_count: int,
        total_truth: Optional[int] = None,
        total_query: Optional[int] = None,
    ) -> QuantifyMetrics:
        """
        Calculate basic benchmarking metrics.

        Args:
            tp_count: True positive count
            fp_count: False positive count
            fn_count: False negative count
            total_truth: Total variants in truth set (if known)
            total_query: Total variants in query set (if known)

        Returns:
            QuantifyMetrics object with calculated values
        """
        # Calculate precision (positive predictive value)
        if tp_count + fp_count > 0:
            precision = tp_count / (tp_count + fp_count)
        else:
            precision = 0.0 if tp_count == 0 else 1.0

        # Calculate recall (sensitivity)
        if tp_count + fn_count > 0:
            recall = tp_count / (tp_count + fn_count)
        else:
            recall = 0.0 if tp_count == 0 else 1.0

        # Calculate F1 score (harmonic mean of precision and recall)
        if precision + recall > 0:
            f1_score = 2 * (precision * recall) / (precision + recall)
        else:
            f1_score = 0.0

        # Use provided totals or derive from counts
        if total_truth is None:
            total_truth = tp_count + fn_count
        if total_query is None:
            total_query = tp_count + fp_count

        return QuantifyMetrics(
            tp_count=tp_count,
            fp_count=fp_count,
            fn_count=fn_count,
            precision=precision,
            recall=recall,
            f1_score=f1_score,
            total_truth=total_truth,
            total_query=total_query,
        )

    @staticmethod
    def calculate_stratified_metrics(
        variant_df: pd.DataFrame,
        stratification_column: str = "stratification",
        variant_type_column: str = "variant_type",
    ) -> Dict[str, Dict[str, QuantifyMetrics]]:
        """
        Calculate metrics stratified by variant type and/or region.

        Args:
            variant_df: DataFrame with variant classifications
            stratification_column: Column name for stratification categories
            variant_type_column: Column name for variant types

        Returns:
            Nested dictionary: {stratification: {variant_type: metrics}}
        """
        results = {}

        # Get unique stratification categories
        strat_categories = variant_df[stratification_column].unique()

        for strat_cat in strat_categories:
            strat_data = variant_df[variant_df[stratification_column] == strat_cat]
            results[strat_cat] = {}

            # Calculate overall metrics for this stratification
            overall_metrics = MetricsCalculator._calculate_metrics_for_subset(
                strat_data
            )
            results[strat_cat]["ALL"] = overall_metrics

            # Calculate metrics by variant type
            variant_types = strat_data[variant_type_column].unique()
            for var_type in variant_types:
                type_data = strat_data[strat_data[variant_type_column] == var_type]
                type_metrics = MetricsCalculator._calculate_metrics_for_subset(
                    type_data
                )
                results[strat_cat][var_type] = type_metrics

        return results

    @staticmethod
    def _calculate_metrics_for_subset(subset_df: pd.DataFrame) -> QuantifyMetrics:
        """Calculate metrics for a subset of variants."""
        decision_col = "benchmark_decision"
        if decision_col not in subset_df.columns:
            for alt_col in ["BD", "benchmark_decision", "decision"]:
                if alt_col in subset_df.columns:
                    decision_col = alt_col
                    break

        decision_counts = subset_df.get(decision_col, pd.Series(dtype=object)).value_counts()

        tp_count = decision_counts.get("TP", 0)
        fp_count = decision_counts.get("FP", 0)
        fn_count = decision_counts.get("FN", 0)

        return MetricsCalculator.calculate_basic_metrics(
            tp_count=tp_count, fp_count=fp_count, fn_count=fn_count
        )

    @staticmethod
    def calculate_threshold_metrics(
        variant_df: pd.DataFrame,
        quality_thresholds: np.ndarray,
        quality_column: str = "quality_score",
    ) -> pd.DataFrame:
        """
        Calculate metrics across multiple quality thresholds for ROC analysis.

        Args:
            variant_df: DataFrame with variant data
            quality_thresholds: Array of quality thresholds to evaluate
            quality_column: Column name containing quality scores

        Returns:
            DataFrame with threshold, TP, FP, FN, precision, recall, f1_score
        """
        results = []

        for threshold in quality_thresholds:
            # Filter variants above threshold
            filtered_df = variant_df[variant_df[quality_column] >= threshold].copy()

            if filtered_df.empty:
                # No variants pass threshold
                results.append(
                    {
                        "threshold": threshold,
                        "TP": 0,
                        "FP": 0,
                        "FN": 0,
                        "precision": 0.0,
                        "recall": 0.0,
                        "f1_score": 0.0,
                        "total_variants": 0,
                    }
                )
                continue

            # Calculate metrics for this threshold
            metrics = MetricsCalculator._calculate_metrics_for_subset(filtered_df)

            results.append(
                {
                    "threshold": threshold,
                    "TP": metrics.tp_count,
                    "FP": metrics.fp_count,
                    "FN": metrics.fn_count,
                    "precision": metrics.precision,
                    "recall": metrics.recall,
                    "f1_score": metrics.f1_score,
                    "total_variants": len(filtered_df),
                }
            )

        return pd.DataFrame(results)

    @staticmethod
    def calculate_summary_metrics(variant_df: pd.DataFrame) -> Dict[str, float]:
        """Return summary metrics for a set of variants."""
        metrics = MetricsCalculator._calculate_metrics_for_subset(variant_df)

        return {
            "truth_total": metrics.total_truth,
            "truth_tp": metrics.tp_count,
            "truth_fn": metrics.fn_count,
            "query_total": metrics.total_query,
            "query_tp": metrics.tp_count,
            "query_fp": metrics.fp_count,
            "fp_gt": 0,
            "fp_al": 0,
            "recall": metrics.recall,
            "precision": metrics.precision,
            "frac_na": 0.0,
            "f1_score": metrics.f1_score,
            "truth_titv": 0.0,
            "query_titv": 0.0,
            "truth_het_hom": 0.0,
            "query_het_hom": 0.0,
        }
