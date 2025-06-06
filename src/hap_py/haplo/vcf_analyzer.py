"""
VCF analysis engine for quantify operations.

This module provides memory-efficient VCF processing for large genomic datasets,
compatible with the original hap.py quantify functionality.
"""

import logging
from pathlib import Path
from typing import Dict, Iterator, List, Optional

import pandas as pd
import pysam

from .quantify_models import (
    BenchmarkDecision,
    StratificationRegion,
    VariantType,
)
from .string_handling import ensure_str

logger = logging.getLogger(__name__)


class VCFAnalyzer:
    """
    Memory-efficient VCF analyzer for quantify operations.

    This class processes VCF files from vcfeval output to extract variant
    classifications and quality metrics for benchmarking analysis.
    """

    def __init__(
        self,
        vcf_path: str,
        quality_field: str = "QUAL",
        chunk_size: int = 10000,
        stratification_regions: Optional[List[StratificationRegion]] = None,
    ):
        """
        Initialize VCF analyzer.

        Args:
            vcf_path: Path to VCF file (typically from vcfeval output)
            quality_field: Quality field to extract (QUAL, GQ, etc.)
            chunk_size: Number of variants to process in each chunk
            stratification_regions: Optional stratification regions for analysis
        """
        self.vcf_path = Path(vcf_path)
        self.quality_field = quality_field
        self.chunk_size = chunk_size
        self.stratification_regions = stratification_regions or []

        # Validate VCF file exists and is readable
        if not self.vcf_path.exists():
            raise FileNotFoundError(f"VCF file not found: {vcf_path}")

    def analyze_variants(self) -> pd.DataFrame:
        """
        Analyze all variants in the VCF file.

        Returns:
            DataFrame with columns: chrom, pos, ref, alt, variant_type,
            benchmark_decision, quality_score, sample_name, stratification
        """
        logger.info(f"Starting VCF analysis: {self.vcf_path}")

        all_variants = []
        total_processed = 0

        try:
            with pysam.VariantFile(str(self.vcf_path)) as vcf:
                # Validate quality field exists in header
                self._validate_quality_field(vcf)

                # Process variants in chunks for memory efficiency
                for chunk in self._process_chunks(vcf):
                    all_variants.extend(chunk)
                    total_processed += len(chunk)

                    if total_processed % 50000 == 0:
                        logger.info(f"Processed {total_processed} variants")

        except Exception as e:
            logger.error(f"Error processing VCF file {self.vcf_path}: {e}")
            raise

        logger.info(f"Completed VCF analysis: {total_processed} variants processed")

        if not all_variants:
            logger.warning("No variants found in VCF file")
            return pd.DataFrame()

        return pd.DataFrame(all_variants)

    def _process_chunks(self, vcf: pysam.VariantFile) -> Iterator[List[Dict]]:
        """Process VCF variants in chunks for memory efficiency."""
        chunk = []

        for record in vcf:
            variant_data = self._extract_variant_data(record)
            if variant_data:
                chunk.append(variant_data)

            if len(chunk) >= self.chunk_size:
                yield chunk
                chunk = []

        if chunk:  # Yield remaining variants
            yield chunk

    def _extract_variant_data(self, record: pysam.VariantRecord) -> Optional[Dict]:
        """
        Extract relevant data from a VCF record.

        This method extracts the key information needed for quantify analysis,
        including variant classification from vcfeval annotations.
        """
        try:
            # Extract basic variant information
            variant_data = {
                "chrom": ensure_str(record.chrom),
                "pos": record.pos,
                "ref": ensure_str(record.ref),
                "alt": ensure_str(record.alts[0]) if record.alts else ".",
                "variant_type": self._classify_variant_type(record),
                "benchmark_decision": self._extract_benchmark_decision(record),
                "quality_score": self._extract_quality_score(record),
                "sample_name": self._extract_sample_name(record),
            }

            # Add stratification information if regions are defined
            if self.stratification_regions:
                variant_data["stratification"] = self._classify_stratification(record)
            else:
                variant_data["stratification"] = "ALL"

            return variant_data

        except Exception as e:
            logger.warning(
                f"Error processing variant at {record.chrom}:{record.pos}: {e}"
            )
            return None

    def _classify_variant_type(self, record: pysam.VariantRecord) -> str:
        """Classify variant type (SNP, INDEL, COMPLEX)."""
        ref_len = len(record.ref)
        alt_len = len(record.alts[0]) if record.alts else 0

        if ref_len == alt_len:
            return VariantType.SNP.value if ref_len == 1 else VariantType.MNP.value
        if ref_len != alt_len:
            return VariantType.INDEL.value
        return VariantType.COMPLEX.value

    def _extract_benchmark_decision(self, record: pysam.VariantRecord) -> str:
        """
        Extract benchmark decision from vcfeval annotations.

        vcfeval typically adds annotations like BD (benchmark decision)
        or uses specific FILTERs to indicate TP/FP/FN status.
        """
        # Check for vcfeval-specific annotations
        if "BD" in record.info:
            decision = ensure_str(record.info["BD"])
            if decision in ["TP", "FP", "FN"]:
                return decision

        # Check FILTER field for vcfeval classifications
        if record.filter.keys():
            filters = list(record.filter.keys())
            if "TP" in filters:
                return BenchmarkDecision.TP.value
            elif "FP" in filters:
                return BenchmarkDecision.FP.value
            elif "FN" in filters:
                return BenchmarkDecision.FN.value

        # Default to unknown if no clear classification
        return BenchmarkDecision.UNK.value

    def _extract_quality_score(self, record: pysam.VariantRecord) -> float:
        """Extract quality score from specified field."""
        if self.quality_field == "QUAL":
            return float(record.qual) if record.qual is not None else 0.0

        # Check INFO field
        if self.quality_field in record.info:
            return float(record.info[self.quality_field])

        # Check FORMAT fields in samples
        for sample in record.samples.values():
            if self.quality_field in sample:
                value = sample[self.quality_field]
                if isinstance(value, (list, tuple)):
                    return float(value[0]) if value else 0.0
                return float(value)

        logger.warning(f"Quality field '{self.quality_field}' not found in record")
        return 0.0

    def _extract_sample_name(self, record: pysam.VariantRecord) -> str:
        """
        Extract sample name to determine if variant is from TRUTH or QUERY.

        vcfeval typically uses specific sample names or INFO annotations
        to distinguish between truth and query variants.
        """
        sample_names = list(record.samples.keys())

        # vcfeval often uses these sample names
        if "TRUTH" in sample_names:
            return "TRUTH"
        elif "QUERY" in sample_names:
            return "QUERY"
        elif len(sample_names) == 1:
            return sample_names[0]

        # Check for vcfeval-specific INFO annotations
        if "SYNC" in record.info:
            sync_value = record.info["SYNC"]
            if isinstance(sync_value, (list, tuple)):
                sync_value = sync_value[0]
            return "TRUTH" if sync_value == 0 else "QUERY"

        return "UNKNOWN"

    def _classify_stratification(self, record: pysam.VariantRecord) -> str:
        """Classify variant into stratification categories."""
        # This is a placeholder - full implementation would use
        # BED file overlap detection or filter expressions
        return "ALL"

    def _validate_quality_field(self, vcf: pysam.VariantFile) -> None:
        """Validate that the specified quality field exists in the VCF."""
        header = vcf.header

        if self.quality_field == "QUAL":
            return  # QUAL is always available

        # Check INFO fields
        if self.quality_field in header.info:
            return

        # Check FORMAT fields
        if self.quality_field in header.formats:
            return

        logger.warning(
            f"Quality field '{self.quality_field}' not found in VCF header. "
            f"This may cause issues during analysis."
        )
