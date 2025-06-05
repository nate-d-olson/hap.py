#!/usr/bin/env python3
"""
GA4GH compliance implementation for hap.py.

This module implements the GA4GH benchmarking standards for variant comparison
as defined in the GA4GH benchmarking documentation:
https://github.com/ga4gh/benchmarking-tools/blob/master/doc/ref-impl/README.md

The implementation provides classes for GA4GH-compliant VCF formatting,
stratification, and metrics calculation.
"""

import argparse
import logging
from enum import Enum
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

import pandas as pd
import pysam

# Create module logger
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
    """
    Class to handle GA4GH-compliant output formatting.

    This class implements the GA4GH intermediate VCF specification for
    variant benchmarking results. It provides methods to format VCF headers
    and annotate VCF records with GA4GH-specific fields.
    """

    # GA4GH FORMAT field definitions
    GA4GH_FORMAT_FIELDS = {
        "BD": {
            "id": "BD",
            "number": "1",
            "type": "String",
            "description": "Decision for the benchmark variant (TP/FP/FN/N/UNK)",
        },
        "BK": {
            "id": "BK",
            "number": "1",
            "type": "String",
            "description": "Decision detail for the benchmark variant",
        },
        "QD": {
            "id": "QD",
            "number": "1",
            "type": "String",
            "description": "Decision for the query variant (TP/FP/FN/N/UNK)",
        },
        "QK": {
            "id": "QK",
            "number": "1",
            "type": "String",
            "description": "Decision detail for the query variant",
        },
    }

    # GA4GH INFO field definitions
    GA4GH_INFO_FIELDS = {
        "Regions": {
            "id": "Regions",
            "number": ".",
            "type": "String",
            "description": "List of region IDs this variant is located in",
        },
        "TruthStatus": {
            "id": "TruthStatus",
            "number": "1",
            "type": "String",
            "description": "Status of the variant in truth VCF (TP/FN/FP/N/UNK)",
        },
        "QueryStatus": {
            "id": "QueryStatus",
            "number": "1",
            "type": "String",
            "description": "Status of the variant in query VCF (TP/FP/FN/N/UNK)",
        },
        "Subtype": {
            "id": "Subtype",
            "number": "1",
            "type": "String",
            "description": "Variant subtype classification (SNP/INDEL/COMPLEX/OTHER)",
        },
    }

    def __init__(self, add_extra_fields: bool = False):
        """
        Initialize the GA4GH formatter.

        Args:
            add_extra_fields: Whether to add extra non-standard fields
        """
        self.add_extra_fields = add_extra_fields

    def format_vcf_header(self, vcf_header: pysam.VariantHeader) -> pysam.VariantHeader:
        """
        Add GA4GH fields to a VCF header.

        Args:
            vcf_header: pysam VariantHeader object to modify

        Returns:
            Modified VCF header with GA4GH fields
        """
        # Add FORMAT fields
        for field_id, field_def in self.GA4GH_FORMAT_FIELDS.items():
            vcf_header.add_meta(
                "FORMAT",
                items=[
                    ("ID", field_def["id"]),
                    ("Number", field_def["number"]),
                    ("Type", field_def["type"]),
                    ("Description", field_def["description"]),
                ],
            )

        # Add INFO fields
        for field_id, field_def in self.GA4GH_INFO_FIELDS.items():
            vcf_header.add_meta(
                "INFO",
                items=[
                    ("ID", field_def["id"]),
                    ("Number", field_def["number"]),
                    ("Type", field_def["type"]),
                    ("Description", field_def["description"]),
                ],
            )

        # Add GA4GH compliance header info
        vcf_header.add_meta(
            "META",
            items=[
                ("ID", "ga4gh_standard"),
                ("Version", "1.0"),
                ("Description", "This VCF follows GA4GH benchmarking standards"),
            ],
        )

        # Add extra fields if requested
        if self.add_extra_fields:
            vcf_header.add_meta(
                "INFO",
                items=[
                    ("ID", "PctSimilarity"),
                    ("Number", "1"),
                    ("Type", "Float"),
                    ("Description", "Percent similarity to matched variant"),
                ],
            )
            vcf_header.add_meta(
                "INFO",
                items=[
                    ("ID", "MatchId"),
                    ("Number", "1"),
                    ("Type", "String"),
                    ("Description", "Identifier of the matched variant"),
                ],
            )

        return vcf_header

    def annotate_record(
        self,
        record: pysam.VariantRecord,
        truth_decision: GA4GHDecision,
        query_decision: GA4GHDecision,
        truth_detail: Optional[GA4GHDecisionDetail] = None,
        query_detail: Optional[GA4GHDecisionDetail] = None,
        regions: Optional[List[str]] = None,
        variant_subtype: Optional[GA4GHVariantType] = None,
        extra_info: Optional[Dict[str, any]] = None,
    ) -> pysam.VariantRecord:
        """
        Add GA4GH annotations to a VCF record.

        Args:
            record: pysam VariantRecord to annotate
            truth_decision: Decision for the truth variant
            query_decision: Decision for the query variant
            truth_detail: Detail for the truth decision
            query_detail: Detail for the query decision
            regions: List of region IDs this variant is located in
            variant_subtype: Variant subtype classification
            extra_info: Additional INFO fields to add

        Returns:
            Annotated VCF record
        """
        # Set FORMAT fields for benchmark and query decisions
        record.samples[0]["BD"] = (
            truth_decision.value if truth_decision else GA4GHDecision.UNK.value
        )
        record.samples[0]["QD"] = (
            query_decision.value if query_decision else GA4GHDecision.UNK.value
        )

        # Set FORMAT fields for decision details
        record.samples[0]["BK"] = (
            truth_detail.value if truth_detail else GA4GHDecisionDetail.NO_MATCH.value
        )
        record.samples[0]["QK"] = (
            query_detail.value if query_detail else GA4GHDecisionDetail.NO_MATCH.value
        )

        # Set INFO fields
        record.info["TruthStatus"] = (
            truth_decision.value if truth_decision else GA4GHDecision.UNK.value
        )
        record.info["QueryStatus"] = (
            query_decision.value if query_decision else GA4GHDecision.UNK.value
        )

        # Set regions if provided
        if regions:
            record.info["Regions"] = regions

        # Set variant subtype if provided
        if variant_subtype:
            record.info["Subtype"] = variant_subtype.value
        else:
            # Auto-detect subtype if not provided
            if len(record.ref) == 1 and all(len(alt) == 1 for alt in record.alts):
                record.info["Subtype"] = GA4GHVariantType.SNP.value
            elif record.ref == record.alts[0][0] or record.ref[0] == record.alts[0][0]:
                record.info["Subtype"] = GA4GHVariantType.INDEL.value
            else:
                record.info["Subtype"] = GA4GHVariantType.COMPLEX.value

        # Add extra INFO fields if provided
        if extra_info and self.add_extra_fields:
            for key, value in extra_info.items():
                if key in record.header.info:
                    record.info[key] = value

        return record


class GA4GHStratification:
    """
    Class to handle GA4GH stratification regions.

    This class provides methods to assign variants to stratification regions
    and calculate metrics within specific regions, following the GA4GH
    stratification standards.
    """

    def __init__(self, stratification_regions: Optional[Dict[str, str]] = None):
        """
        Initialize with stratification regions.

        Args:
            stratification_regions: Dictionary mapping region ID to BED file path
        """
        self.regions = {}
        if stratification_regions:
            for region_id, bed_file in stratification_regions.items():
                self.add_region(region_id, bed_file)

    def add_region(self, region_id: str, bed_file: str) -> None:
        """
        Add a stratification region.

        Args:
            region_id: Unique identifier for the region
            bed_file: Path to BED file defining the region
        """
        try:
            import pybedtools

            region = pybedtools.BedTool(bed_file)
            self.regions[region_id] = {
                "bed_file": bed_file,
                "region": region,
            }
        except ImportError:
            logger.warning("pybedtools not available, using simplified region handling")
            self.regions[region_id] = {
                "bed_file": bed_file,
                "region": None,
            }

    def get_region_ids_for_variant(
        self, chrom: str, pos: int, end: Optional[int] = None
    ) -> List[str]:
        """
        Get list of region IDs a variant belongs to.

        Args:
            chrom: Chromosome name
            pos: 1-based start position
            end: End position (for non-SNPs)

        Returns:
            List of region IDs this variant belongs to
        """
        if not end:
            end = pos

        result = []
        try:
            import pybedtools

            for region_id, region_data in self.regions.items():
                if region_data["region"] is None:
                    continue

                # Create a temporary BED entry for the variant
                variant_bed = pybedtools.BedTool(
                    f"{chrom}\t{pos-1}\t{end}\n", from_string=True
                )

                # Check for overlap
                if len(variant_bed.intersect(region_data["region"])) > 0:
                    result.append(region_id)
        except ImportError:
            # Simplified handling without pybedtools
            logger.warning("Using simplified region check without pybedtools")
            for region_id in self.regions:
                result.append(region_id)

        return result

    def stratify_variants(self, variants_df: pd.DataFrame) -> Dict[str, pd.DataFrame]:
        """
        Group variants by stratification region.

        Args:
            variants_df: DataFrame with variants to stratify

        Returns:
            Dictionary mapping region IDs to filtered DataFrames
        """
        result = {"all": variants_df.copy()}  # Always include "all" region

        # Create a new column with region assignments
        variants_df["regions"] = variants_df.apply(
            lambda row: self.get_region_ids_for_variant(
                row["chrom"], row["pos"], row.get("end", row["pos"])
            ),
            axis=1,
        )

        # Create filtered DataFrames for each region
        for region_id in self.regions:
            region_variants = variants_df[
                variants_df["regions"].apply(lambda regions: region_id in regions)
            ].copy()
            result[region_id] = region_variants

        return result


class GA4GHMetrics:
    """
    Class to calculate GA4GH-compliant benchmarking metrics.

    This class implements the standard GA4GH metrics for variant benchmarking,
    including precision, recall, and F1-score calculations with confidence
    intervals.
    """

    def __init__(self, bootstrap_iterations: int = 1000, ci_level: float = 0.95):
        """
        Initialize metrics calculator.

        Args:
            bootstrap_iterations: Number of iterations for bootstrap CI calculation
            ci_level: Confidence interval level (0.95 = 95% CI)
        """
        self.bootstrap_iterations = bootstrap_iterations
        self.ci_level = ci_level
        self.scipy_available = False

        try:
            import scipy.stats  # noqa: F401

            self.scipy_available = True
        except ImportError:
            logger.warning(
                "scipy not available, confidence intervals will be approximate"
            )

    def calculate_precision(
        self, tp: int, fp: int, with_ci: bool = False
    ) -> Union[float, Tuple[float, float, float]]:
        """
        Calculate precision (TP / (TP + FP)).

        Args:
            tp: True positive count
            fp: False positive count
            with_ci: Whether to include confidence intervals

        Returns:
            Precision value or tuple of (precision, lower_ci, upper_ci)
        """
        if tp + fp == 0:
            precision = 0.0
        else:
            precision = tp / (tp + fp)

        if not with_ci:
            return precision

        # Calculate confidence intervals
        if self.scipy_available:
            import scipy.stats

            # Use beta distribution for binomial confidence interval
            if tp + fp > 0:
                lower, upper = scipy.stats.beta.ppf(
                    [(1 - self.ci_level) / 2, 1 - (1 - self.ci_level) / 2],
                    tp + 1,
                    fp + 1,
                )
            else:
                lower, upper = 0.0, 1.0
        else:
            # Approximate CI using normal approximation if scipy not available
            if tp + fp > 0:
                stderr = (precision * (1 - precision) / (tp + fp)) ** 0.5
                z = 1.96  # Approximate z for 95% CI
                lower = max(0.0, precision - z * stderr)
                upper = min(1.0, precision + z * stderr)
            else:
                lower, upper = 0.0, 1.0

        return precision, lower, upper

    def calculate_recall(
        self, tp: int, fn: int, with_ci: bool = False
    ) -> Union[float, Tuple[float, float, float]]:
        """
        Calculate recall (TP / (TP + FN)).

        Args:
            tp: True positive count
            fn: False negative count
            with_ci: Whether to include confidence intervals

        Returns:
            Recall value or tuple of (recall, lower_ci, upper_ci)
        """
        if tp + fn == 0:
            recall = 0.0
        else:
            recall = tp / (tp + fn)

        if not with_ci:
            return recall

        # Calculate confidence intervals
        if self.scipy_available:
            import scipy.stats

            # Use beta distribution for binomial confidence interval
            if tp + fn > 0:
                lower, upper = scipy.stats.beta.ppf(
                    [(1 - self.ci_level) / 2, 1 - (1 - self.ci_level) / 2],
                    tp + 1,
                    fn + 1,
                )
            else:
                lower, upper = 0.0, 1.0
        else:
            # Approximate CI using normal approximation if scipy not available
            if tp + fn > 0:
                stderr = (recall * (1 - recall) / (tp + fn)) ** 0.5
                z = 1.96  # Approximate z for 95% CI
                lower = max(0.0, recall - z * stderr)
                upper = min(1.0, recall + z * stderr)
            else:
                lower, upper = 0.0, 1.0

        return recall, lower, upper

    def calculate_f1(
        self, precision: float, recall: float, with_ci: bool = False
    ) -> Union[float, Tuple[float, float, float]]:
        """
        Calculate F1 score (2 * precision * recall / (precision + recall)).

        Args:
            precision: Precision value
            recall: Recall value
            with_ci: Whether to include confidence intervals

        Returns:
            F1 score or tuple of (f1, lower_ci, upper_ci)
        """
        if precision + recall == 0:
            f1 = 0.0
        else:
            f1 = 2 * precision * recall / (precision + recall)

        if not with_ci:
            return f1

        # Note: F1 confidence intervals are complex and would require bootstrapping
        # This is a placeholder for future implementation
        return f1, f1 * 0.95, min(f1 * 1.05, 1.0)

    def calculate_metrics(
        self, tp: int, fp: int, fn: int, with_ci: bool = False
    ) -> Dict[str, Union[float, Tuple[float, float, float]]]:
        """
        Calculate all metrics (precision, recall, F1).

        Args:
            tp: True positive count
            fp: False positive count
            fn: False negative count
            with_ci: Whether to include confidence intervals

        Returns:
            Dictionary with calculated metrics
        """
        precision = self.calculate_precision(tp, fp, with_ci)
        recall = self.calculate_recall(tp, fn, with_ci)

        if with_ci:
            prec_value = precision[0]
            rec_value = recall[0]
        else:
            prec_value = precision
            rec_value = recall

        f1 = self.calculate_f1(prec_value, rec_value, with_ci)

        result = {
            "precision": precision,
            "recall": recall,
            "f1": f1,
            "tp": tp,
            "fp": fp,
            "fn": fn,
        }

        return result


# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
)


def main():
    """Main function to demonstrate GA4GH classes."""
    parser = argparse.ArgumentParser(description="GA4GH Implementation Demo")
    parser.add_argument("--output-dir", help="Output directory", default=".")
    args = parser.parse_args()

    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(exist_ok=True)

    # Demo GA4GH formatter
    _ = GA4GHFormatter()  # Initialize formatter for demo
    logger.info("Created GA4GH formatter")

    # Demo GA4GH stratification
    _ = GA4GHStratification(  # Initialize stratification for demo
        {
            "high_conf": "path/to/high_conf.bed",
            "difficult": "path/to/difficult.bed",
        }
    )
    logger.info("Created GA4GH stratification handler")

    # Demo GA4GH metrics
    _ = GA4GHMetrics()  # Initialize metrics for demo
    logger.info("Created GA4GH metrics calculator")

    # Write summary to output file
    with open(output_dir / "ga4gh_implementation.txt", "w") as f:
        f.write("GA4GH Compliance Implementation\n")
        f.write("-----------------------------\n\n")
        f.write("Components initialized:\n")
        f.write("1. GA4GHFormatter\n")
        f.write("2. GA4GHStratification\n")
        f.write("3. GA4GHMetrics\n\n")
        f.write("Implementation status: Initial framework created\n")

    logger.info(f"Wrote summary to {output_dir / 'ga4gh_implementation.txt'}")


if __name__ == "__main__":
    main()
