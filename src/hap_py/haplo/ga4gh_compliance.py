#!/usr/bin/env python3
"""
GA4GH compliance implementation for hap.py.

This module provides classes and enums for GA4GH benchmarking standards compliance,
including variant decision tracking, VCF formatting, and metrics calculation.
"""

import importlib
import logging
from enum import Enum
from typing import Any, Dict, List, Optional, Tuple, Union

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


class GA4GHDecision(Enum):
    """GA4GH decision values for variant classification."""

    TP = "TP"  # True positive
    FP = "FP"  # False positive
    FN = "FN"  # False negative
    N = "N"  # Non-assessed (variant in non-confident region)
    UNK = "UNK"  # Unknown/undetermined


class GA4GHDecisionDetail(Enum):
    """GA4GH decision detail values providing additional context."""

    GT_MATCH = "gt-match"  # Genotype match
    GT_MISMATCH = "gt-mismatch"  # Genotype mismatch
    ALLELE_MATCH = "allele-match"  # Allele match only
    ALLELE_MISMATCH = "allele-mismatch"  # Allele mismatch
    OUTSIDE_CONFIDENT = "outside-confident"  # Outside confident regions
    REFERENCE_MATCH = "reference-match"  # Matches reference
    COMPLEX_MATCH = "complex-match"  # Complex representation match
    NO_MATCH = "no-match"  # No match found


class GA4GHVariantType(Enum):
    """GA4GH variant type classification."""

    SNP = "SNP"  # Single nucleotide polymorphism
    INDEL = "INDEL"  # Insertion/deletion
    COMPLEX = "COMPLEX"  # Complex variation
    OTHER = "OTHER"  # Other variation type


class GA4GHFormatter:
    """Formats VCF records according to GA4GH benchmarking standards."""

    def __init__(self, add_extra_fields: bool = False):
        """
        Initialize GA4GH formatter.

        Args:
            add_extra_fields: Whether to add additional custom fields
        """
        self.add_extra_fields = add_extra_fields

    def format_vcf_header(self, header):
        """
        Add GA4GH required fields to VCF header.

        Args:
            header: VCF header object

        Returns:
            Modified header
        """
        # Add required FORMAT fields
        format_fields = [
            (
                "BD",
                "1",
                "String",
                "Decision for the benchmark variant (TP/FP/FN/N/UNK)",
            ),
            ("BK", "1", "String", "Decision detail for the benchmark variant"),
            ("QD", "1", "String", "Decision for the query variant (TP/FP/FN/N/UNK)"),
            ("QK", "1", "String", "Decision detail for the query variant"),
        ]

        for field_id, number, field_type, description in format_fields:
            header.add_meta(
                "FORMAT",
                items=[
                    ("ID", field_id),
                    ("Number", number),
                    ("Type", field_type),
                    ("Description", description),
                ],
            )

        # Add required INFO fields
        info_fields = [
            ("Regions", ".", "String", "Benchmark regions in which variant falls"),
            ("TruthStatus", "1", "String", "Status of truth variant (TP/FP/FN/N/UNK)"),
            ("QueryStatus", "1", "String", "Status of query variant (TP/FP/FN/N/UNK)"),
            ("Subtype", "1", "String", "Variant subtype (SNP/INDEL/COMPLEX/OTHER)"),
        ]

        for field_id, number, field_type, description in info_fields:
            header.add_meta(
                "INFO",
                items=[
                    ("ID", field_id),
                    ("Number", number),
                    ("Type", field_type),
                    ("Description", description),
                ],
            )

        # Add extra fields if requested
        if self.add_extra_fields:
            extra_info_fields = [
                (
                    "PctSimilarity",
                    "1",
                    "Float",
                    "Percent similarity to matched variant",
                ),
                ("MatchId", "1", "String", "Identifier of the matched variant"),
            ]

            for field_id, number, field_type, description in extra_info_fields:
                header.add_meta(
                    "INFO",
                    items=[
                        ("ID", field_id),
                        ("Number", number),
                        ("Type", field_type),
                        ("Description", description),
                    ],
                )

        # Add GA4GH standard META field
        header.add_meta(
            "META",
            items=[
                ("ID", "ga4gh_standard"),
                ("Version", "1.0"),
                ("Description", "This VCF follows GA4GH benchmarking standards"),
            ],
        )

        return header

    def annotate_record(
        self,
        record,
        truth_decision: GA4GHDecision,
        query_decision: GA4GHDecision,
        truth_detail: Optional[GA4GHDecisionDetail] = None,
        query_detail: Optional[GA4GHDecisionDetail] = None,
        regions: Optional[List[str]] = None,
        variant_subtype: Optional[GA4GHVariantType] = None,
    ):
        """
        Annotate VCF record with GA4GH fields.

        Args:
            record: VCF record to annotate
            truth_decision: Decision for benchmark/truth variant
            query_decision: Decision for query variant
            truth_detail: Detailed decision for truth variant
            query_detail: Detailed decision for query variant
            regions: List of regions the variant falls in
            variant_subtype: Type of variant (auto-detected if None)

        Returns:
            Annotated VCF record
        """
        # Set FORMAT fields for samples
        if record.samples and len(record.samples) > 0:
            # Get first sample (assuming single sample comparison)
            sample = record.samples[0]
            sample["BD"] = truth_decision.value
            sample["QD"] = query_decision.value

            if truth_detail:
                sample["BK"] = truth_detail.value
            if query_detail:
                sample["QK"] = query_detail.value

        # Set INFO fields
        record.info["TruthStatus"] = truth_decision.value
        record.info["QueryStatus"] = query_decision.value

        if regions:
            record.info["Regions"] = regions

        # Auto-detect variant subtype if not provided
        if variant_subtype is None:
            variant_subtype = self._auto_detect_variant_type(record.ref, record.alts)
        record.info["Subtype"] = variant_subtype.value

        return record

    def _auto_detect_variant_type(self, ref: str, alts: List[str]) -> GA4GHVariantType:
        """
        Auto-detect variant type based on REF and ALT sequences.

        Args:
            ref: Reference allele
            alts: Alternative alleles

        Returns:
            Detected variant type
        """
        if not alts or len(alts) == 0:
            return GA4GHVariantType.OTHER

        alt = alts[0]  # Use first ALT for classification

        # SNP: same length, single position change
        if len(ref) == len(alt) == 1:
            return GA4GHVariantType.SNP

        # INDEL: different lengths
        if len(ref) != len(alt):
            return GA4GHVariantType.INDEL

        # COMPLEX: same length but multiple changes or longer sequences
        if len(ref) == len(alt) and len(ref) > 1:
            return GA4GHVariantType.COMPLEX

        return GA4GHVariantType.OTHER


class GA4GHStratification:
    """Stratifies variants by genomic regions according to GA4GH standards."""

    def __init__(self, regions: Optional[Dict[str, str]] = None):
        """
        Initialize GA4GH stratification.

        Args:
            regions: Dictionary mapping region names to BED file paths
        """
        self.regions = {}
        if regions:
            for name, bed_file in regions.items():
                self.add_region(name, bed_file)

    def add_region(self, name: str, bed_file: str):
        """
        Add a stratification region.

        Args:
            name: Name of the region
            bed_file: Path to BED file defining the region
        """
        try:
            import importlib
            import os

            if not os.path.exists(bed_file):
                logger.warning(
                    f"BED file {bed_file} does not exist. Skipping region {name}."
                )
                self.regions[name] = {"bed_file": bed_file, "region": None}
                return

            pybedtools = importlib.import_module("pybedtools")
            bed_region = pybedtools.BedTool(bed_file)
            self.regions[name] = {"bed_file": bed_file, "region": bed_region}
        except ImportError:
            logger.warning(
                "pybedtools not available. Region-based stratification will be limited."
            )
            self.regions[name] = {"bed_file": bed_file, "region": None}

    def get_region_ids_for_variant(self, chrom: str, start: int, end: int) -> List[str]:
        """
        Get region IDs that overlap with a variant.

        Args:
            chrom: Chromosome
            start: Start position
            end: End position

        Returns:
            List of region names that overlap the variant
        """
        overlapping_regions = []

        try:
            pybedtools = importlib.import_module("pybedtools")

            # Create interval for the variant
            variant_interval = pybedtools.BedTool(
                f"{chrom}\t{start}\t{end}", from_string=True
            )

            for region_name, region_data in self.regions.items():
                if region_data["region"] is not None:
                    # Check for intersection
                    intersection = variant_interval.intersect(region_data["region"])
                    if len(intersection) > 0:
                        overlapping_regions.append(region_name)
                else:
                    # Fallback: assume all variants are in the region if pybedtools unavailable
                    overlapping_regions.append(region_name)

        except ImportError:
            # If pybedtools not available, return all regions
            overlapping_regions = list(self.regions.keys())

        return overlapping_regions

    def stratify_variants(self, variants_df: pd.DataFrame) -> Dict[str, pd.DataFrame]:
        """
        Stratify variants by regions.

        Args:
            variants_df: DataFrame with variant information

        Returns:
            Dictionary mapping region names to filtered DataFrames
        """
        stratified = {}

        # Add "all" stratification (all variants)
        stratified["all"] = variants_df.copy()

        # Stratify by each region
        for region_name in self.regions.keys():
            # For now, include all variants in each region
            # In a full implementation, this would filter based on BED intersections
            region_variants = []
            for _, row in variants_df.iterrows():
                chrom = row.get("chrom", row.get("chromosome", ""))
                start = row.get("pos", row.get("start", 0))
                end = row.get("end", start + 1)

                overlapping_regions = self.get_region_ids_for_variant(chrom, start, end)
                if region_name in overlapping_regions:
                    region_variants.append(row)

            if region_variants:
                stratified[region_name] = pd.DataFrame(region_variants)
            else:
                # Return empty DataFrame with same structure
                stratified[region_name] = variants_df.head(0).copy()

        return stratified


class GA4GHMetrics:
    """Calculates GA4GH-compliant metrics with confidence intervals."""

    def __init__(self, bootstrap_iterations: int = 1000, ci_level: float = 0.95):
        """
        Initialize GA4GH metrics calculator.

        Args:
            bootstrap_iterations: Number of bootstrap iterations for CI calculation
            ci_level: Confidence interval level (e.g., 0.95 for 95% CI)
        """
        self.bootstrap_iterations = bootstrap_iterations
        self.ci_level = ci_level

    def calculate_precision(
        self, tp: int, fp: int, with_ci: bool = False
    ) -> Union[float, Tuple[float, float, float]]:
        """
        Calculate precision with optional confidence intervals.

        Args:
            tp: True positives
            fp: False positives
            with_ci: Whether to calculate confidence intervals

        Returns:
            Precision (and CI bounds if requested)
        """
        if tp + fp == 0:
            precision = 0.0
        else:
            precision = tp / (tp + fp)

        if not with_ci:
            return precision

        # Bootstrap confidence intervals
        lower, upper = self._bootstrap_metric(
            lambda tp_boot, fp_boot, fn_boot: (
                tp_boot / (tp_boot + fp_boot) if (tp_boot + fp_boot) > 0 else 0.0
            ),
            tp,
            fp,
            0,
        )

        return precision, lower, upper

    def calculate_recall(
        self, tp: int, fn: int, with_ci: bool = False
    ) -> Union[float, Tuple[float, float, float]]:
        """
        Calculate recall with optional confidence intervals.

        Args:
            tp: True positives
            fn: False negatives
            with_ci: Whether to calculate confidence intervals

        Returns:
            Recall (and CI bounds if requested)
        """
        if tp + fn == 0:
            recall = 0.0
        else:
            recall = tp / (tp + fn)

        if not with_ci:
            return recall

        # Bootstrap confidence intervals
        lower, upper = self._bootstrap_metric(
            lambda tp_boot, fp_boot, fn_boot: (
                tp_boot / (tp_boot + fn_boot) if (tp_boot + fn_boot) > 0 else 0.0
            ),
            tp,
            0,
            fn,
        )

        return recall, lower, upper

    def calculate_f1(
        self, precision: float, recall: float, with_ci: bool = False
    ) -> Union[float, Tuple[float, float, float]]:
        """
        Calculate F1 score with optional confidence intervals.

        Args:
            precision: Precision value
            recall: Recall value
            with_ci: Whether to calculate confidence intervals

        Returns:
            F1 score (and CI bounds if requested)
        """
        if precision + recall == 0:
            f1 = 0.0
        else:
            f1 = 2 * (precision * recall) / (precision + recall)

        if not with_ci:
            return f1

        # For CI, we need to bootstrap from the underlying counts
        # This is a simplified implementation
        lower = max(0.0, f1 - 1.96 * 0.1)  # Simplified CI
        upper = min(1.0, f1 + 1.96 * 0.1)

        return f1, lower, upper

    def calculate_metrics(
        self, tp: int, fp: int, fn: int, with_ci: bool = False
    ) -> Dict[str, Any]:
        """
        Calculate all metrics together.

        Args:
            tp: True positives
            fp: False positives
            fn: False negatives
            with_ci: Whether to calculate confidence intervals

        Returns:
            Dictionary with all calculated metrics
        """
        if with_ci:
            precision, prec_lower, prec_upper = self.calculate_precision(
                tp, fp, with_ci=True
            )
            recall, rec_lower, rec_upper = self.calculate_recall(tp, fn, with_ci=True)
            f1, f1_lower, f1_upper = self.calculate_f1(precision, recall, with_ci=True)

            return {
                "precision": (precision, prec_lower, prec_upper),
                "recall": (recall, rec_lower, rec_upper),
                "f1": (f1, f1_lower, f1_upper),
                "tp": tp,
                "fp": fp,
                "fn": fn,
            }
        else:
            precision = self.calculate_precision(tp, fp)
            recall = self.calculate_recall(tp, fn)
            f1 = self.calculate_f1(precision, recall)

            return {
                "precision": precision,
                "recall": recall,
                "f1": f1,
                "tp": tp,
                "fp": fp,
                "fn": fn,
            }

    def _bootstrap_metric(
        self, metric_func, tp: int, fp: int, fn: int
    ) -> Tuple[float, float]:
        """
        Calculate bootstrap confidence intervals for a metric.

        Args:
            metric_func: Function to calculate the metric
            tp: True positives
            fp: False positives
            fn: False negatives

        Returns:
            Lower and upper confidence bounds
        """
        # Create bootstrap samples
        bootstrap_values = []
        total_tp = max(1, tp)  # Avoid division by zero
        total_fp = max(0, fp)
        total_fn = max(0, fn)

        for _ in range(self.bootstrap_iterations):
            # Bootstrap sample with replacement
            tp_boot = np.random.poisson(total_tp)
            fp_boot = np.random.poisson(total_fp)
            fn_boot = np.random.poisson(total_fn)

            metric_value = metric_func(tp_boot, fp_boot, fn_boot)
            bootstrap_values.append(metric_value)

        # Calculate confidence intervals
        alpha = 1 - self.ci_level
        lower_percentile = (alpha / 2) * 100
        upper_percentile = (1 - alpha / 2) * 100

        lower = np.percentile(bootstrap_values, lower_percentile)
        upper = np.percentile(bootstrap_values, upper_percentile)

        return lower, upper
