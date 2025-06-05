#!/usr/bin/env python3
"""
GA4GH integration for the QuantifyEngine.

This module provides integration between the GA4GH compliance classes and
the QuantifyEngine to enable GA4GH-compliant variant benchmarking.
"""

import logging
from pathlib import Path
from typing import Dict, Optional, Tuple, Union

import pandas as pd
import pysam

from hap_py.haplo.ga4gh_compliance import (
    GA4GHDecision,
    GA4GHDecisionDetail,
    GA4GHFormatter,
    GA4GHMetrics,
    GA4GHStratification,
    GA4GHVariantType,
)

logger = logging.getLogger(__name__)


class GA4GHIntegration:
    """
    Class to integrate GA4GH compliance with QuantifyEngine.

    This class provides methods to integrate GA4GH standards with the
    existing QuantifyEngine class, enabling GA4GH-compliant variant
    benchmarking and output formatting.
    """

    def __init__(
        self,
        stratification_beds: Optional[Dict[str, str]] = None,
        confidence_regions: Optional[str] = None,
        bootstrap_iterations: int = 1000,
        ci_level: float = 0.95,
    ):
        """
        Initialize GA4GH integration.

        Args:
            stratification_beds: Dictionary mapping region IDs to BED file paths
            confidence_regions: Path to confident regions BED file
            bootstrap_iterations: Number of iterations for bootstrap CI calculation
            ci_level: Confidence interval level (0.95 = 95% CI)
        """
        # Initialize formatter
        self.formatter = GA4GHFormatter()

        # Initialize stratification
        self.stratification = GA4GHStratification(stratification_beds)
        if confidence_regions:
            # Add high-confidence regions if provided
            self.stratification.add_region("CONF", confidence_regions)

        # Initialize metrics
        self.metrics = GA4GHMetrics(
            bootstrap_iterations=bootstrap_iterations,
            ci_level=ci_level,
        )

    def prepare_vcf_header(
        self, header_template: pysam.VariantHeader
    ) -> pysam.VariantHeader:
        """
        Prepare a GA4GH-compliant VCF header.

        Args:
            header_template: Base VCF header to modify

        Returns:
            GA4GH-compliant VCF header
        """
        return self.formatter.format_vcf_header(header_template)

    def transform_match_to_ga4gh(
        self, match: tuple, truth_df: pd.DataFrame, query_df: pd.DataFrame
    ) -> Tuple[GA4GHDecision, GA4GHDecision, GA4GHDecisionDetail, GA4GHDecisionDetail]:
        """
        Transform a match tuple to GA4GH decisions and details.

        Args:
            match: Match tuple (truth_idx, query_idx, match_type)
            truth_df: DataFrame with truth variants
            query_df: DataFrame with query variants

        Returns:
            Tuple of (truth_decision, query_decision, truth_detail, query_detail)
        """
        truth_idx, query_idx, match_type = match

        # Default values
        truth_decision = GA4GHDecision.UNK
        query_decision = GA4GHDecision.UNK
        truth_detail = GA4GHDecisionDetail.NO_MATCH
        query_detail = GA4GHDecisionDetail.NO_MATCH

        # True positive match
        if truth_idx >= 0 and query_idx >= 0:
            truth_decision = GA4GHDecision.TP
            query_decision = GA4GHDecision.TP

            # Determine match details
            if "gt" in match_type or "genotype" in match_type:
                truth_detail = GA4GHDecisionDetail.GT_MATCH
                query_detail = GA4GHDecisionDetail.GT_MATCH
            elif "allele" in match_type:
                truth_detail = GA4GHDecisionDetail.ALLELE_MATCH
                query_detail = GA4GHDecisionDetail.ALLELE_MATCH
            elif "complex" in match_type:
                truth_detail = GA4GHDecisionDetail.COMPLEX_MATCH
                query_detail = GA4GHDecisionDetail.COMPLEX_MATCH

        # False negative (truth variant with no match)
        elif truth_idx >= 0 and query_idx < 0:
            truth_decision = GA4GHDecision.FN
            query_decision = GA4GHDecision.UNK
            truth_detail = GA4GHDecisionDetail.NO_MATCH

        # False positive (query variant with no match)
        elif truth_idx < 0 and query_idx >= 0:
            truth_decision = GA4GHDecision.UNK
            query_decision = GA4GHDecision.FP
            query_detail = GA4GHDecisionDetail.NO_MATCH

        return (truth_decision, query_decision, truth_detail, query_detail)

    def annotate_vcf_record(
        self,
        record: pysam.VariantRecord,
        truth_decision: GA4GHDecision,
        query_decision: GA4GHDecision,
        truth_detail: GA4GHDecisionDetail,
        query_detail: GA4GHDecisionDetail,
        chrom: str,
        pos: int,
        end: Optional[int] = None,
    ) -> pysam.VariantRecord:
        """
        Annotate a VCF record with GA4GH-compliant fields.

        Args:
            record: VCF record to annotate
            truth_decision: Decision for truth variant
            query_decision: Decision for query variant
            truth_detail: Detail for truth decision
            query_detail: Detail for query decision
            chrom: Chromosome name
            pos: Variant position
            end: Variant end position

        Returns:
            Annotated VCF record
        """
        # Get regions for this variant
        regions = self.stratification.get_region_ids_for_variant(chrom, pos, end)

        # Determine variant subtype
        if len(record.ref) == 1 and all(len(alt) == 1 for alt in record.alts):
            variant_subtype = GA4GHVariantType.SNP
        elif record.ref == record.alts[0][0] or record.ref[0] == record.alts[0][0]:
            variant_subtype = GA4GHVariantType.INDEL
        else:
            variant_subtype = GA4GHVariantType.COMPLEX

        # Annotate record
        return self.formatter.annotate_record(
            record,
            truth_decision=truth_decision,
            query_decision=query_decision,
            truth_detail=truth_detail,
            query_detail=query_detail,
            regions=regions,
            variant_subtype=variant_subtype,
        )

    def create_ga4gh_metrics(
        self, tp: int, fp: int, fn: int, region: str = "all", with_ci: bool = True
    ) -> Dict[str, Union[float, Tuple[float, float, float]]]:
        """
        Calculate GA4GH-compliant metrics.

        Args:
            tp: True positive count
            fp: False positive count
            fn: False negative count
            region: Region ID for these metrics
            with_ci: Whether to include confidence intervals

        Returns:
            Dictionary with metrics
        """
        # Calculate metrics with the GA4GH metrics class
        metrics = self.metrics.calculate_metrics(tp, fp, fn, with_ci=with_ci)

        # Add region information
        metrics["region"] = region

        return metrics

    def write_ga4gh_metrics_file(
        self, metrics_by_region: Dict[str, Dict], output_path: Union[str, Path]
    ) -> None:
        """
        Write GA4GH metrics to a file.

        Args:
            metrics_by_region: Dictionary mapping region IDs to metrics
            output_path: Path to write metrics file
        """
        output_path = Path(output_path)

        # Create header row
        if any(
            "precision" in metrics and isinstance(metrics["precision"], tuple)
            for metrics in metrics_by_region.values()
        ):
            # With confidence intervals
            header = [
                "Region",
                "Type",
                "TP",
                "FP",
                "FN",
                "Precision",
                "Precision_lower",
                "Precision_upper",
                "Recall",
                "Recall_lower",
                "Recall_upper",
                "F1",
                "F1_lower",
                "F1_upper",
            ]
        else:
            # Without confidence intervals
            header = [
                "Region",
                "Type",
                "TP",
                "FP",
                "FN",
                "Precision",
                "Recall",
                "F1",
            ]

        rows = [header]

        # Add rows for each region and type
        for region, metrics in metrics_by_region.items():
            variant_types = metrics.get("types", {"ALL": metrics})

            for variant_type, type_metrics in variant_types.items():
                if isinstance(type_metrics.get("precision"), tuple):
                    # With confidence intervals
                    row = [
                        region,
                        variant_type,
                        type_metrics["tp"],
                        type_metrics["fp"],
                        type_metrics["fn"],
                        type_metrics["precision"][0],
                        type_metrics["precision"][1],
                        type_metrics["precision"][2],
                        type_metrics["recall"][0],
                        type_metrics["recall"][1],
                        type_metrics["recall"][2],
                        type_metrics["f1"][0],
                        type_metrics["f1"][1],
                        type_metrics["f1"][2],
                    ]
                else:
                    # Without confidence intervals
                    row = [
                        region,
                        variant_type,
                        type_metrics["tp"],
                        type_metrics["fp"],
                        type_metrics["fn"],
                        type_metrics["precision"],
                        type_metrics["recall"],
                        type_metrics["f1"],
                    ]

                rows.append(row)

        # Write to file
        with open(output_path, "w") as f:
            for row in rows:
                f.write("\t".join(map(str, row)) + "\n")

        logger.info(f"Wrote GA4GH metrics to {output_path}")


def enhance_quantify_engine_with_ga4gh(engine):
    """
    Enhance a QuantifyEngine instance with GA4GH compliance.

    Args:
        engine: QuantifyEngine instance to enhance

    Returns:
        Enhanced QuantifyEngine
    """
    # Add GA4GH integration
    engine.ga4gh = GA4GHIntegration(
        stratification_beds=engine.stratification_regions,
        confidence_regions=engine.confident_regions,
    )

    # Enhance output methods
    original_write_vcf = engine._write_vcf_outputs

    def enhanced_write_vcf(self, output_prefix, *args, **kwargs):
        """Enhanced version of _write_vcf_outputs with GA4GH support."""
        result = original_write_vcf(output_prefix, *args, **kwargs)

        if self.quantify_method == "ga4gh":
            # Write additional GA4GH-compliant files
            try:
                # Create GA4GH VCF output
                ga4gh_metrics_path = f"{output_prefix}.ga4gh.metrics.tsv"

                # Write GA4GH metrics
                metrics_by_region = {}
                for region, results in self.results_by_region.items():
                    metrics_by_region[region] = self.ga4gh.create_ga4gh_metrics(
                        results["TP"], results["FP"], results["FN"], region=region
                    )

                self.ga4gh.write_ga4gh_metrics_file(
                    metrics_by_region, ga4gh_metrics_path
                )

                logger.info(
                    f"GA4GH-compliant outputs written to {output_prefix}.ga4gh.*"
                )
            except Exception as e:
                logger.error(f"Failed to write GA4GH-compliant outputs: {e}")

        return result

    # Replace the method
    engine._write_vcf_outputs = lambda *args, **kwargs: enhanced_write_vcf(
        engine, *args, **kwargs
    )

    # Enhance the _track_benchmarking_decisions method
    original_track_decisions = engine._track_benchmarking_decisions

    def enhanced_track_decisions(self, truth_df, query_df):
        """Enhanced version of _track_benchmarking_decisions with GA4GH support."""
        result = original_track_decisions(truth_df, query_df)

        if self.quantify_method == "ga4gh":
            # Add GA4GH-specific decision tracking
            logger.info("Adding GA4GH-specific decision tracking")

            # Track GA4GH decisions in truth_df
            for i, row in truth_df.iterrows():
                truth_decision = GA4GHDecision.UNK
                detail = GA4GHDecisionDetail.NO_MATCH

                if row.get("BD") == "TP":
                    truth_decision = GA4GHDecision.TP
                    detail = GA4GHDecisionDetail.GT_MATCH
                elif row.get("BD") == "FN":
                    truth_decision = GA4GHDecision.FN
                    detail = GA4GHDecisionDetail.NO_MATCH
                elif row.get("BD") == "FP":
                    truth_decision = GA4GHDecision.FP
                    detail = GA4GHDecisionDetail.NO_MATCH
                elif row.get("BD") == "N":
                    truth_decision = GA4GHDecision.N
                    detail = GA4GHDecisionDetail.OUTSIDE_CONFIDENT

                truth_df.at[i, "GA4GH_Decision"] = truth_decision.value
                truth_df.at[i, "GA4GH_Detail"] = detail.value

            # Track GA4GH decisions in query_df
            for i, row in query_df.iterrows():
                query_decision = GA4GHDecision.UNK
                detail = GA4GHDecisionDetail.NO_MATCH

                if row.get("BD") == "TP":
                    query_decision = GA4GHDecision.TP
                    detail = GA4GHDecisionDetail.GT_MATCH
                elif row.get("BD") == "FP":
                    query_decision = GA4GHDecision.FP
                    detail = GA4GHDecisionDetail.NO_MATCH
                elif row.get("BD") == "FN":
                    query_decision = GA4GHDecision.FN
                    detail = GA4GHDecisionDetail.NO_MATCH
                elif row.get("BD") == "N":
                    query_decision = GA4GHDecision.N
                    detail = GA4GHDecisionDetail.OUTSIDE_CONFIDENT

                query_df.at[i, "GA4GH_Decision"] = query_decision.value
                query_df.at[i, "GA4GH_Detail"] = detail.value

        return result

    # Replace the method
    engine._track_benchmarking_decisions = (
        lambda *args, **kwargs: enhanced_track_decisions(engine, *args, **kwargs)
    )

    return engine
