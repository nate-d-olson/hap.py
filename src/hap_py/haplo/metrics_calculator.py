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
        decision_counts = subset_df["benchmark_decision"].value_counts()

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
    def _titv_ratio(variants: pd.DataFrame) -> float:
        """Calculate transition/transversion ratio for SNPs."""
        if variants is None or variants.empty:
            return float("nan")

        snps = variants[
            (variants["REF"].str.len() == 1) & (variants["ALT"].str.len() == 1)
        ]
        if snps.empty:
            return float("nan")

        transitions = 0
        transversions = 0
        for ref, alt in zip(snps["REF"].str.upper(), snps["ALT"].str.upper()):
            pair = (ref, alt)
            if pair in [("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")]:
                transitions += 1
            else:
                transversions += 1

        if transversions == 0:
            return float("inf") if transitions > 0 else float("nan")

        return transitions / transversions

    @staticmethod
    def _het_hom_ratio(gt_series: pd.Series) -> float:
        """Calculate heterozygous/homozygous ratio from genotype strings."""
        if gt_series is None or gt_series.empty:
            return float("nan")

        het = 0
        hom = 0
        for gt in gt_series.astype(str):
            if gt in {"./.", ".", ""}:
                continue
            gts = gt.replace("|", "/").split("/")
            if len(gts) < 2:
                continue
            if gts[0] == gts[1]:
                hom += 1
            else:
                het += 1

        if hom == 0:
            return float("inf") if het > 0 else float("nan")

        return het / hom

    @staticmethod
    def calculate_summary_metrics(df: pd.DataFrame) -> Dict[str, float]:
        """Calculate summary metrics from a GA4GH annotated DataFrame."""

        tp = int(df.get("TP", 0).sum()) if "TP" in df.columns else 0
        fp = int(df.get("FP", 0).sum()) if "FP" in df.columns else 0
        fn = int(df.get("FN", 0).sum()) if "FN" in df.columns else 0

        truth_total = tp + fn
        query_total = tp + fp

        metrics: Dict[str, float] = {
            "truth_total": float(truth_total),
            "truth_tp": float(tp),
            "truth_fn": float(fn),
            "query_total": float(query_total),
            "query_tp": float(tp),
            "query_fp": float(fp),
        }

        if "BK" in df.columns:
            metrics["fp_gt"] = float(
                (
                    (df.get("FP") == True) & (df["BK"].astype(str).str.contains("gm"))
                ).sum()
            )
            metrics["fp_al"] = float(
                (
                    (df.get("FP") == True) & (df["BK"].astype(str).str.contains("am"))
                ).sum()
            )
        else:
            metrics["fp_gt"] = 0.0
            metrics["fp_al"] = 0.0

        precision = tp / query_total if query_total > 0 else 0.0
        recall = tp / truth_total if truth_total > 0 else 0.0
        f1 = (
            2 * precision * recall / (precision + recall)
            if precision + recall > 0
            else 0.0
        )

        metrics.update(
            {
                "precision": precision,
                "recall": recall,
                "f1_score": f1,
            }
        )

        if "BD" in df.columns and len(df) > 0:
            na_count = df["BD"].isna().sum() + (df["BD"] == ".").sum()
            metrics["frac_na"] = na_count / len(df)
        else:
            metrics["frac_na"] = 0.0

        truth_variants = df[(df.get("TP", False)) | (df.get("FN", False))]
        query_variants = df[(df.get("TP", False)) | (df.get("FP", False))]

        metrics["truth_titv"] = MetricsCalculator._titv_ratio(truth_variants)
        metrics["query_titv"] = MetricsCalculator._titv_ratio(query_variants)
        metrics["truth_het_hom"] = MetricsCalculator._het_hom_ratio(
            truth_variants.get("gt")
        )
        metrics["query_het_hom"] = MetricsCalculator._het_hom_ratio(
            query_variants.get("gt")
        )

        return metrics
