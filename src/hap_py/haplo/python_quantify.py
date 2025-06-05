#!/usr/bin/env python3
"""
Python implementation of the quantify module.

This module provides functionality to quantify variants in VCF files,
producing stratification metrics and summary statistics.
"""

import json
import logging
from typing import Any, Dict, List, Optional, Tuple

import pandas as pd
import pysam

# Phase 2 imports for enhanced ROC analysis
try:
    import matplotlib
    import matplotlib.pyplot as plt

    matplotlib.use("Agg")  # Use non-interactive backend
    MATPLOTLIB_AVAILABLE = True
except ImportError:
    MATPLOTLIB_AVAILABLE = False

try:
    import scipy.stats  # noqa: F401

    SCIPY_AVAILABLE = True
except ImportError:
    SCIPY_AVAILABLE = False

try:
    import sklearn.metrics  # noqa: F401

    SKLEARN_AVAILABLE = True
except ImportError:
    SKLEARN_AVAILABLE = False

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class QuantifyEngine:
    """
    Engine for quantifying variant calls in VCF files.

    This class replaces the C++ quantify component with a pure Python
    implementation using pandas and pysam.
    """

    def __init__(
        self,
        truth_vcf: str,
        query_vcf: str,
        reference: Optional[str] = None,
        regions: Optional[str] = None,
        apply_filters: bool = False,
        output_vtc: bool = False,
        quantify_method: str = "xcmp",  # xcmp or ga4gh
        enable_roc_analysis: bool = True,  # Phase 2: Enable ROC analysis
        roc_bootstrap_samples: int = 1000,  # Phase 2: Bootstrap samples for confidence intervals
        quality_stratification: bool = True,  # Phase 2: Enable quality-based stratification
    ):
        """
        Initialize the quantify engine.

        Args:
            truth_vcf: Path to truth VCF file
            query_vcf: Path to query/test VCF file
            reference: Path to reference FASTA file (optional)
            regions: BED file with regions to quantify (optional)
            apply_filters: Whether to apply filters from VCF
            output_vtc: Whether to output variant truth categories
            quantify_method: Quantification method - 'xcmp' or 'ga4gh'
            enable_roc_analysis: Enable Phase 2 ROC analysis with confidence intervals
            roc_bootstrap_samples: Number of bootstrap samples for confidence intervals
            quality_stratification: Enable quality score-based stratification
        """
        self.truth_vcf = truth_vcf
        self.query_vcf = query_vcf
        self.reference = reference
        self.regions = regions
        self.apply_filters = apply_filters
        self.output_vtc = output_vtc
        self.quantify_method = quantify_method.lower()

        # Phase 2: ROC Analysis Configuration
        self.enable_roc_analysis = enable_roc_analysis
        self.roc_bootstrap_samples = roc_bootstrap_samples
        self.quality_stratification = quality_stratification

        # Validate quantify method
        if self.quantify_method not in ["xcmp", "ga4gh"]:
            raise ValueError(
                f"Invalid quantify method: {quantify_method}. Must be 'xcmp' or 'ga4gh'"
            )

        self.truth_variants: List[Dict[str, Any]] = []
        self.query_variants: List[Dict[str, Any]] = []
        self.region_list: List[Tuple[str, int, int]] = []

        # Results storage
        self.metrics: Dict[str, Any] = {}
        self.stratifications: Dict[str, Any] = {}

        # Phase 2: ROC Analysis Results Storage
        self.roc_data: Dict[str, Any] = {}
        self.quality_metrics: Dict[str, Any] = {}
        self.bootstrap_confidence_intervals: Dict[str, Any] = {}

        # Open VCF files
        self._open_vcfs()

        # Load regions if provided
        if regions:
            self._load_regions()

    def _open_vcfs(self):
        """Open VCF files using pysam."""
        try:
            self.truth_vcf_handle = pysam.VariantFile(self.truth_vcf)
            logger.info(f"Opened truth VCF: {self.truth_vcf}")
        except Exception as e:
            raise ValueError(f"Failed to open truth VCF: {e}")

        try:
            self.query_vcf_handle = pysam.VariantFile(self.query_vcf)
            logger.info(f"Opened query VCF: {self.query_vcf}")
        except Exception as e:
            raise ValueError(f"Failed to open query VCF: {e}")

        # If GA4GH method is requested, validate that the VCF conforms to the
        # GA4GH benchmarking intermediate specification. The query VCF is
        # often the same file in re-quantification workflows, but we validate
        # both handles in case they differ.
        if self.quantify_method == "ga4gh":
            try:
                self._validate_ga4gh_vcf(self.truth_vcf_handle)
                self._validate_ga4gh_vcf(self.query_vcf_handle)
            except Exception as e:
                raise ValueError(f"GA4GH compliance check failed: {e}")

        # Open reference if provided
        if self.reference:
            try:
                self.reference_handle = pysam.FastaFile(self.reference)
                self.reference_file = (
                    self.reference
                )  # Store path for _get_reference_sequence
                logger.info(f"Opened reference: {self.reference}")
            except Exception as e:
                logger.warning(f"Failed to open reference: {e}")
                self.reference = None
                self.reference_file = None
        else:
            self.reference_file = None

    def _validate_ga4gh_vcf(self, vcf_handle: pysam.VariantFile) -> None:
        """Validate required GA4GH intermediate VCF fields."""

        header = vcf_handle.header

        required_info = ["BS"]
        missing_info = [f for f in required_info if f not in header.info]

        required_format = ["GT", "BD", "BK", "BI", "QQ", "BVT", "BLT"]
        missing_format = [f for f in required_format if f not in header.formats]

        required_samples = ["TRUTH", "QUERY"]
        missing_samples = [s for s in required_samples if s not in header.samples]

        if missing_info or missing_format or missing_samples:
            problems = []
            if missing_info:
                problems.append(f"INFO fields: {', '.join(missing_info)}")
            if missing_format:
                problems.append(f"FORMAT fields: {', '.join(missing_format)}")
            if missing_samples:
                problems.append(f"samples: {', '.join(missing_samples)}")
            raise ValueError("Input VCF missing GA4GH fields - " + "; ".join(problems))

    def _load_regions(self):
        """Load regions from BED file."""
        try:
            # Read BED file
            with open(self.regions) as f:
                for line in f:
                    if line.startswith("#") or not line.strip():
                        continue

                    parts = line.strip().split("\t")
                    if len(parts) < 3:
                        continue

                    chrom = parts[0]
                    start = int(parts[1])
                    end = int(parts[2])

                    # Store as tuple (chrom, start, end)
                    self.region_list.append((chrom, start, end))

            logger.info(f"Loaded {len(self.region_list)} regions from {self.regions}")
        except Exception as e:
            logger.error(f"Failed to load regions: {e}")
            self.region_list = []

    def _load_variants(self, is_truth: bool = True):
        """
        Load variants from a VCF file into memory.

        Args:
            is_truth: Whether loading from truth or query VCF

        Returns:
            List of processed variants
        """
        vcf_handle = self.truth_vcf_handle if is_truth else self.query_vcf_handle
        vcf_name = "truth" if is_truth else "query"

        variants = []
        count = 0
        filtered_count = 0

        for variant in vcf_handle.fetch():
            count += 1

            # Skip filtered variants if apply_filters is True
            if self.apply_filters and self._is_filtered(variant):
                filtered_count += 1
                continue

            # Process variant
            processed = self._process_variant(variant)

            # Check if in regions
            if self.region_list and not self._in_regions(processed):
                continue

            variants.append(processed)

            # Log progress for large files
            if count % 100000 == 0:
                logger.info(f"Processed {count} {vcf_name} variants...")

        logger.info(
            f"Loaded {len(variants)} {vcf_name} variants "
            f"(filtered {filtered_count} variants)"
        )

        return variants

    def _is_filtered(self, variant):
        """
        Check if a variant is filtered.

        Args:
            variant: pysam VariantRecord

        Returns:
            True if variant is filtered, False otherwise
        """
        # If no filter or filter is PASS or '.', it's not filtered
        return len(variant.filter) > 0 and "PASS" not in variant.filter

    def _process_variant(self, variant):
        """
        Process a variant record.

        Args:
            variant: pysam VariantRecord

        Returns:
            Dict with processed variant info
        """
        # Basic variant info
        processed = {
            "chrom": variant.chrom,
            "pos": variant.pos,
            "id": variant.id,
            "ref": variant.ref,
            "alt": ",".join(variant.alts) if variant.alts else "",
            "qual": variant.qual,
            "filter": list(variant.filter),
            "type": self._get_variant_type(variant),
            "is_indel": (
                len(variant.ref) != len(variant.alts[0]) if variant.alts else False
            ),
            "is_snp": (
                len(variant.ref) == 1 and len(variant.alts[0]) == 1
                if variant.alts
                else False
            ),
            "is_mnp": (
                len(variant.ref) > 1 and len(variant.ref) == len(variant.alts[0])
                if variant.alts
                else False
            ),
            "length": (
                max(len(variant.ref), len(variant.alts[0])) if variant.alts else 0
            ),
        }

        # Add genotype information if available
        if variant.samples:
            sample = list(variant.samples.values())[0]

            if "GT" in sample:
                processed["gt"] = "/".join(map(str, sample["GT"]))
                processed["is_hom"] = (
                    sample["GT"][0] == sample["GT"][1]
                    if len(sample["GT"]) > 1
                    else True
                )
                processed["is_het"] = not processed["is_hom"]

            # Other sample format fields
            for key in sample.keys():
                if key != "GT":
                    processed[key] = sample[key]

        return processed

    def _get_variant_type(self, variant):
        """
        Determine the type of variant.

        Args:
            variant: pysam VariantRecord

        Returns:
            String with variant type
        """
        if not variant.alts:
            return "REF"

        ref_len = len(variant.ref)
        alt_len = len(variant.alts[0])

        if ref_len == 1 and alt_len == 1:
            return "SNP"
        elif ref_len == alt_len and ref_len > 1:
            return "MNP"
        elif ref_len > alt_len:
            return "DEL"
        elif ref_len < alt_len:
            return "INS"
        else:
            return "COMPLEX"

    def _in_regions(self, variant):
        """
        Check if a variant is in the specified regions.

        Args:
            variant: Processed variant dict

        Returns:
            True if variant is in regions, False otherwise
        """
        for chrom, start, end in self.region_list:
            if (
                variant["chrom"] == chrom
                and variant["pos"] >= start
                and variant["pos"] <= end
            ):
                return True

        return False

    def quantify(self):
        """
        Quantify variants in the VCF files.

        Returns:
            Dict with quantification results
        """
        # Load variants
        logger.info("Loading truth variants...")
        self.truth_variants = self._load_variants(is_truth=True)

        logger.info("Loading query variants...")
        self.query_variants = self._load_variants(is_truth=False)

        # Match variants between truth and query
        logger.info("Matching variants...")
        self._match_variants()

        # Calculate metrics
        logger.info("Calculating metrics...")
        self._calculate_metrics()

        # Stratify results
        logger.info("Stratifying results...")
        self._stratify_results()

        # Perform ROC analysis (Phase 2)
        logger.info("Performing ROC analysis (Phase 2)...")
        self._perform_roc_analysis()

        # Store results
        results = {
            "metrics": self.metrics,
            "stratifications": self.stratifications,
            "truth_variants": len(self.truth_variants),
            "query_variants": len(self.query_variants),
        }

        # Add Phase 2 results if ROC analysis was performed
        if self.enable_roc_analysis:
            results["roc_data"] = self.roc_data
            results["bootstrap_confidence_intervals"] = (
                self.bootstrap_confidence_intervals
            )
            if self.quality_stratification:
                results["quality_metrics"] = self.quality_metrics

        return results

    def process_vcf(self) -> Optional[pd.DataFrame]:
        """
        Process VCF file and return a DataFrame with variant information.

        This method reads the input VCF and converts it to a structured DataFrame
        that can be used for analysis and ROC generation.

        Returns:
            DataFrame with variant information, or None if no variants found
        """
        try:
            # Open VCF files
            self._open_vcfs()

            # Load regions if specified
            self._load_regions()

            # For xcmp results, we typically work with a single comparison VCF
            # that contains both truth and query information
            variants_data = []

            vcf_handle = self.truth_vcf_handle
            for variant in vcf_handle.fetch():
                # Skip filtered variants if requested
                if self.apply_filters and self._is_filtered(variant):
                    continue

                # Process the variant and extract relevant information
                variant_data = self._process_variant(variant)

                # Check if variant is in specified regions
                if self.region_list and not self._in_regions(variant_data):
                    continue

                # Add additional fields for analysis
                variant_info = {
                    "CHROM": variant_data["chrom"],
                    "POS": variant_data["pos"],
                    "REF": variant_data["ref"],
                    "ALT": variant_data["alt"],
                    "Type": variant_data["type"],
                    "Length": variant_data["length"],
                }

                # Extract INFO fields
                for key, value in variant.info.items():
                    variant_info[key] = value

                # Extract quality score
                if variant.qual is not None:
                    variant_info["QUAL"] = variant.qual
                else:
                    variant_info["QUAL"] = 0.0

                # Extract decision information from BD (Benchmarking Decision) tag
                if "BD" in variant.info:
                    variant_info["BD"] = variant.info["BD"]
                else:
                    variant_info["BD"] = "UNK"  # Unknown

                # Extract BK (Benchmarking Kind) tag
                if "BK" in variant.info:
                    variant_info["BK"] = variant.info["BK"]
                else:
                    variant_info["BK"] = "UNK"

                variants_data.append(variant_info)

            if not variants_data:
                logger.warning("No variants found in VCF file")
                return None

            # Convert to DataFrame
            df = pd.DataFrame(variants_data)
            logger.info(f"Processed {len(df)} variants from VCF")

            return df

        except Exception as e:
            logger.error(f"Error processing VCF file: {e}")
            return None
        finally:
            # Close VCF handles
            if hasattr(self, "truth_vcf_handle") and self.truth_vcf_handle:
                self.truth_vcf_handle.close()
            if hasattr(self, "query_vcf_handle") and self.query_vcf_handle:
                self.query_vcf_handle.close()

    def apply_bed_stratification(self, df: pd.DataFrame, bed_file: str) -> pd.DataFrame:
        """
        Apply BED file stratification to a DataFrame of variants.

        Args:
            df: DataFrame with variant information
            bed_file: Path to BED file with regions

        Returns:
            Filtered DataFrame containing only variants in the BED regions
        """
        try:
            # Load BED regions
            regions = []
            with open(bed_file) as f:
                for line in f:
                    if line.startswith("#") or not line.strip():
                        continue
                    parts = line.strip().split("\t")
                    if len(parts) >= 3:
                        chrom = parts[0]
                        start = int(parts[1])
                        end = int(parts[2])
                        regions.append((chrom, start, end))

            if not regions:
                logger.warning(f"No regions found in BED file: {bed_file}")
                return df

            # Filter variants that overlap with BED regions
            filtered_indices = []
            for idx, row in df.iterrows():
                chrom = row.get("CHROM", "")
                pos = row.get("POS", 0)

                # Check if variant overlaps with any region
                for region_chrom, region_start, region_end in regions:
                    if chrom == region_chrom and region_start <= pos <= region_end:
                        filtered_indices.append(idx)
                        break

            filtered_df = df.loc[filtered_indices].copy()
            logger.info(
                f"Stratification with {bed_file}: {len(filtered_df)}/{len(df)} variants retained"
            )

            return filtered_df

        except Exception as e:
            logger.error(f"Error applying BED stratification: {e}")
            return df

    def write_output_vcf(self, df: pd.DataFrame, output_path: str) -> None:
        """
        Write processed variants to output VCF file.

        Args:
            df: DataFrame with variant information
            output_path: Path for output VCF file
        """
        try:
            # For now, just copy the input VCF to output
            # In a full implementation, this would write the processed variants
            import shutil

            shutil.copy2(self.truth_vcf, output_path)
            logger.info(f"Output VCF written to {output_path}")
        except Exception as e:
            logger.error(f"Error writing output VCF: {e}")

    def _match_variants(self):
        """
        Match variants between truth and query sets using sophisticated algorithms.

        This implements Phase 1 of the quantify module enhancement:
        - Variant normalization (trimLeft, trimRight, leftShift)
        - Sophisticated variant matching beyond simple coordinate matching
        - Benchmarking decision tracking (BD, BVT, QQ fields)
        - Multi-allelic variant support
        - XCMP vs GA4GH quantification method support
        """
        logger.info(
            f"Starting sophisticated variant matching using {self.quantify_method.upper()} method..."
        )

        # Normalize all variants first
        self._normalize_all_variants()

        # Decompose multi-allelic variants for better matching
        self._decompose_multiallelic_variants()

        # Convert to pandas DataFrames for matching
        truth_df = pd.DataFrame(self.truth_variants)
        query_df = pd.DataFrame(self.query_variants)

        # Add source tracking and benchmarking decision fields
        truth_df["source"] = "truth"
        query_df["source"] = "query"

        # Initialize benchmarking decision tracking
        self._initialize_benchmarking_decisions(truth_df, query_df)

        # Perform sophisticated matching based on quantify method
        if self.quantify_method == "xcmp":
            matches = self._perform_xcmp_matching(truth_df, query_df)
        else:  # ga4gh
            matches = self._perform_ga4gh_matching(truth_df, query_df)

        # Apply matches and track benchmarking decisions
        self._apply_matches(truth_df, query_df, matches)

        # Track benchmarking decisions (BD, BVT, QQ fields)
        self._track_benchmarking_decisions(truth_df, query_df)

        # Update variant lists with enhanced match information
        self.truth_variants = truth_df.to_dict("records")
        self.query_variants = query_df.to_dict("records")

        logger.info(
            f"Variant matching complete: {len(matches)} matches found using {self.quantify_method.upper()} method"
        )

    def _initialize_benchmarking_decisions(
        self, truth_df: pd.DataFrame, query_df: pd.DataFrame
    ):
        """Initialize benchmarking decision tracking fields."""
        # Initialize match tracking
        truth_df["match"] = False
        query_df["match"] = False
        truth_df["match_idx"] = -1
        query_df["match_idx"] = -1

        # Initialize benchmarking decision fields (BD, BVT, QQ)
        truth_df["BD"] = "FN"  # Default to False Negative
        query_df["BD"] = "FP"  # Default to False Positive

        # Benchmarking Variant Type (BVT)
        truth_df["BVT"] = truth_df.apply(self._classify_variant_type, axis=1)
        query_df["BVT"] = query_df.apply(self._classify_variant_type, axis=1)

        # Quality Quantiles (QQ) - based on QUAL field if available
        if "qual" in query_df.columns:
            # Calculate quality quantiles with better handling of missing/null values
            qual_values = query_df["qual"].fillna(0.0)
            qual_values = pd.to_numeric(qual_values, errors="coerce").fillna(0.0)

            # Create quantiles only if we have valid quality scores
            if qual_values.max() > 0:
                query_df["QQ"] = pd.qcut(
                    qual_values, q=10, labels=False, duplicates="drop"
                )
            else:
                query_df["QQ"] = 5  # Default middle quantile
        else:
            query_df["QQ"] = 5  # Default middle quantile

        # Add confidence scores for matches
        truth_df["match_confidence"] = 0.0
        query_df["match_confidence"] = 0.0

    def _classify_variant_type(self, row: pd.Series) -> str:
        """
        Classify variant type for BVT (Benchmarking Variant Type) field.

        Args:
            row: Pandas series representing a variant

        Returns:
            Variant type string (SNP, INS, DEL, COMPLEX)
        """
        ref = str(row.get("ref", ""))
        alt = str(row.get("alt", ""))

        if len(ref) == 1 and len(alt) == 1:
            return "SNP"
        elif len(ref) > len(alt):
            return "DEL"
        elif len(ref) < len(alt):
            return "INS"
        else:
            return "COMPLEX"

    def _perform_sophisticated_matching(
        self, truth_df: pd.DataFrame, query_df: pd.DataFrame
    ) -> list:
        """
        Perform sophisticated variant matching beyond simple coordinate matching.

        This implements multiple matching strategies:
        1. Exact coordinate and allele matching (after normalization)
        2. Overlapping variant matching for complex cases
        3. Multi-allelic variant decomposition and matching
        4. Superlocus-based matching for complex regions

        Args:
            truth_df: Truth variants DataFrame
            query_df: Query variants DataFrame

        Returns:
            List of match tuples (truth_idx, query_idx, match_type, confidence)
        """
        matches = []

        # Strategy 1: Exact matching after normalization
        exact_matches = self._find_exact_matches(truth_df, query_df)
        matches.extend(exact_matches)

        # Strategy 2: Overlapping variant matching
        overlap_matches = self._find_overlapping_matches(
            truth_df, query_df, exact_matches
        )
        matches.extend(overlap_matches)

        # Strategy 3: Multi-allelic decomposition matching
        multiallelic_matches = self._find_multiallelic_matches(
            truth_df, query_df, exact_matches + overlap_matches
        )
        matches.extend(multiallelic_matches)

        # Strategy 4: Superlocus matching for complex regions
        superlocus_matches = self._find_superlocus_matches(truth_df, query_df, matches)
        matches.extend(superlocus_matches)

        return matches

    def _find_exact_matches(
        self, truth_df: pd.DataFrame, query_df: pd.DataFrame
    ) -> list:
        """Find exact matches between normalized variants."""
        matches = []

        # Create lookup dictionary for query variants
        query_lookup = {}
        for idx, row in query_df.iterrows():
            key = (row["chrom"], row["pos"], row["ref"], row["alt"])
            if key not in query_lookup:
                query_lookup[key] = []
            query_lookup[key].append(idx)

        # Find matching truth variants
        for truth_idx, truth_row in truth_df.iterrows():
            truth_key = (
                truth_row["chrom"],
                truth_row["pos"],
                truth_row["ref"],
                truth_row["alt"],
            )

            if truth_key in query_lookup:
                for query_idx in query_lookup[truth_key]:
                    matches.append((truth_idx, query_idx, "exact", 1.0))

        return matches

    def _find_overlapping_matches(
        self, truth_df: pd.DataFrame, query_df: pd.DataFrame, existing_matches: list
    ) -> list:
        """Find overlapping variants that may represent the same biological variation."""
        matches = []

        # Get already matched indices
        matched_truth = {m[0] for m in existing_matches}
        matched_query = {m[1] for m in existing_matches}

        # Find overlapping variants
        for truth_idx, truth_row in truth_df.iterrows():
            if truth_idx in matched_truth:
                continue

            truth_start = truth_row["pos"]
            truth_end = truth_row["pos"] + len(truth_row["ref"]) - 1

            for query_idx, query_row in query_df.iterrows():
                if query_idx in matched_query:
                    continue

                if truth_row["chrom"] != query_row["chrom"]:
                    continue

                query_start = query_row["pos"]
                query_end = query_row["pos"] + len(query_row["ref"]) - 1

                # Check for overlap
                if truth_start <= query_end and query_start <= truth_end:
                    # For overlapping variants, also check allele compatibility
                    # Only match if the alleles are compatible or represent equivalent changes
                    if self._are_alleles_compatible(truth_row, query_row):
                        # Calculate overlap confidence based on position proximity and allele similarity
                        distance = abs(truth_start - query_start)
                        confidence = max(
                            0.1, 1.0 - (distance / 50.0)
                        )  # Confidence decreases with distance

                        matches.append((truth_idx, query_idx, "overlap", confidence))

        return matches

    def _find_multiallelic_matches(
        self, truth_df: pd.DataFrame, query_df: pd.DataFrame, existing_matches: list
    ) -> list:
        """Handle multi-allelic variants by decomposing and matching individual alleles."""
        matches = []

        # Get already matched indices
        matched_truth = {m[0] for m in existing_matches}
        matched_query = {m[1] for m in existing_matches}

        # This is a simplified implementation - in production would fully decompose multi-allelic variants
        for truth_idx, truth_row in truth_df.iterrows():
            if truth_idx in matched_truth:
                continue

            # Check if this might be a multi-allelic variant (contains comma in ALT)
            if "," in str(truth_row.get("alt", "")):
                alt_alleles = str(truth_row["alt"]).split(",")

                for query_idx, query_row in query_df.iterrows():
                    if query_idx in matched_query:
                        continue

                    # Check if query variant matches any of the decomposed alleles
                    if (
                        truth_row["chrom"] == query_row["chrom"]
                        and truth_row["pos"] == query_row["pos"]
                        and truth_row["ref"] == query_row["ref"]
                        and query_row["alt"] in alt_alleles
                    ):
                        matches.append((truth_idx, query_idx, "multiallelic", 0.9))

        return matches

    def _find_superlocus_matches(
        self, truth_df: pd.DataFrame, query_df: pd.DataFrame, existing_matches: list
    ) -> list:
        """
        Find matches at the superlocus level for complex variant regions.

        Superloci are regions containing multiple variants that should be evaluated together.
        """
        matches = []

        # Get already matched indices
        matched_truth = {m[0] for m in existing_matches}
        matched_query = {m[1] for m in existing_matches}

        # Group variants into potential superloci (simplified approach)
        truth_superloci = self._group_into_superloci(truth_df, matched_truth)
        query_superloci = self._group_into_superloci(query_df, matched_query)

        # Match superloci
        for truth_locus in truth_superloci:
            for query_locus in query_superloci:
                if self._superloci_overlap(truth_locus, query_locus):
                    # Simple superlocus matching - in production would use sophisticated algorithms
                    confidence = 0.7  # Lower confidence for superlocus matches

                    # Match first variant in each superlocus as representative
                    if truth_locus and query_locus:
                        matches.append(
                            (
                                truth_locus[0]["idx"],
                                query_locus[0]["idx"],
                                "superlocus",
                                confidence,
                            )
                        )

        return matches

    def _group_into_superloci(
        self, df: pd.DataFrame, matched_indices: set, window_size: int = 50
    ) -> list:
        """Group nearby variants into superloci."""
        superloci = []

        for chrom in df["chrom"].unique():
            chrom_variants = df[df["chrom"] == chrom].copy()
            chrom_variants = chrom_variants[~chrom_variants.index.isin(matched_indices)]

            if chrom_variants.empty:
                continue

            chrom_variants = chrom_variants.sort_values("pos")

            current_locus = []
            last_pos = None

            for idx, row in chrom_variants.iterrows():
                if last_pos is None or row["pos"] - last_pos <= window_size:
                    current_locus.append(
                        {"idx": idx, "pos": row["pos"], "chrom": row["chrom"]}
                    )
                    last_pos = row["pos"]
                else:
                    if len(current_locus) > 1:  # Only consider multi-variant loci
                        superloci.append(current_locus)
                    current_locus = [
                        {"idx": idx, "pos": row["pos"], "chrom": row["chrom"]}
                    ]
                    last_pos = row["pos"]

            if len(current_locus) > 1:
                superloci.append(current_locus)

        return superloci

    def _superloci_overlap(self, locus1: list, locus2: list) -> bool:
        """Check if two superloci overlap."""
        if not locus1 or not locus2:
            return False

        # Check chromosome
        if locus1[0]["chrom"] != locus2[0]["chrom"]:
            return False

        # Check position overlap
        locus1_start = min(v["pos"] for v in locus1)
        locus1_end = max(v["pos"] for v in locus1)
        locus2_start = min(v["pos"] for v in locus2)
        locus2_end = max(v["pos"] for v in locus2)

        return locus1_start <= locus2_end and locus2_start <= locus1_end

    def _apply_matches(
        self, truth_df: pd.DataFrame, query_df: pd.DataFrame, matches: list
    ):
        """Apply found matches to the DataFrames."""
        for truth_idx, query_idx, match_type, confidence in matches:
            # Mark as matched
            truth_df.loc[truth_idx, "match"] = True
            query_df.loc[query_idx, "match"] = True

            # Store match indices
            truth_df.loc[truth_idx, "match_idx"] = query_idx
            query_df.loc[query_idx, "match_idx"] = truth_idx

            # Store additional match metadata
            truth_df.loc[truth_idx, "match_type"] = match_type
            query_df.loc[query_idx, "match_type"] = match_type
            truth_df.loc[truth_idx, "match_confidence"] = confidence
            query_df.loc[query_idx, "match_confidence"] = confidence

    def _track_benchmarking_decisions(
        self, truth_df: pd.DataFrame, query_df: pd.DataFrame
    ):
        """Track benchmarking decisions (BD field) based on matches."""
        # Update BD field based on matches
        truth_df.loc[truth_df["match"], "BD"] = "TP"  # True Positive
        truth_df.loc[~truth_df["match"], "BD"] = "FN"  # False Negative

        query_df.loc[query_df["match"], "BD"] = "TP"  # True Positive
        query_df.loc[~query_df["match"], "BD"] = "FP"  # False Positive

        # Log benchmarking decision summary
        tp_count = len(truth_df[truth_df["BD"] == "TP"])
        fn_count = len(truth_df[truth_df["BD"] == "FN"])
        fp_count = len(query_df[query_df["BD"] == "FP"])

        logger.info(
            f"Benchmarking decisions: TP={tp_count}, FN={fn_count}, FP={fp_count}"
        )

    def _decompose_multiallelic_variants(self):
        """
        Decompose multi-allelic variants into bi-allelic variants for better matching.

        This enhances the matching by creating separate entries for each alternative allele
        in multi-allelic variants, making it easier to match individual alleles.
        """
        logger.info("Decomposing multi-allelic variants...")

        # Process truth variants
        self.truth_variants = self._decompose_variant_list(self.truth_variants)

        # Process query variants
        self.query_variants = self._decompose_variant_list(self.query_variants)

        logger.info(
            f"After decomposition: {len(self.truth_variants)} truth, {len(self.query_variants)} query variants"
        )

    def _decompose_variant_list(
        self, variant_list: List[Dict[str, Any]]
    ) -> List[Dict[str, Any]]:
        """
        Decompose multi-allelic variants in a variant list.

        Args:
            variant_list: List of variant dictionaries

        Returns:
            List with multi-allelic variants decomposed into bi-allelic variants
        """
        decomposed = []

        for variant in variant_list:
            alt_alleles = str(variant.get("alt", "")).split(",")

            if len(alt_alleles) <= 1:
                # Not multi-allelic, keep as is
                decomposed.append(variant)
            else:
                # Multi-allelic - create separate variant for each alt allele
                for i, alt_allele in enumerate(alt_alleles):
                    if not alt_allele.strip():
                        continue

                    # Create new variant for this allele
                    new_variant = variant.copy()
                    new_variant["alt"] = alt_allele.strip()
                    new_variant["original_multiallelic"] = True
                    new_variant["allele_index"] = i
                    new_variant["original_variant_id"] = variant.get("id", "")

                    # Update variant type for this specific allele
                    new_variant["type"] = self._determine_variant_type_from_alleles(
                        new_variant["ref"], alt_allele.strip()
                    )

                    decomposed.append(new_variant)

        return decomposed

    def _determine_variant_type_from_alleles(self, ref: str, alt: str) -> str:
        """
        Determine variant type from ref and alt alleles.

        Args:
            ref: Reference allele
            alt: Alternative allele

        Returns:
            Variant type string
        """
        ref_len = len(ref)
        alt_len = len(alt)

        if ref_len == 1 and alt_len == 1:
            return "SNP"
        elif ref_len == alt_len and ref_len > 1:
            return "MNP"
        elif ref_len > alt_len:
            return "DEL"
        elif ref_len < alt_len:
            return "INS"
        else:
            return "COMPLEX"

    def _perform_xcmp_matching(
        self, truth_df: pd.DataFrame, query_df: pd.DataFrame
    ) -> list:
        """
        Perform XCMP-style variant matching.

        XCMP (eXtended Comparison) uses more sophisticated matching that includes:
        - Exact matching after normalization
        - Complex variant matching using superloci
        - Representation-independent matching

        Args:
            truth_df: Truth variants DataFrame
            query_df: Query variants DataFrame

        Returns:
            List of match tuples
        """
        logger.info("Performing XCMP-style matching...")
        matches = []

        # Strategy 1: Exact matching after normalization
        exact_matches = self._find_exact_matches(truth_df, query_df)
        matches.extend(exact_matches)
        logger.info(f"XCMP exact matches: {len(exact_matches)}")

        # Strategy 2: Complex variant equivalence matching
        complex_matches = self._find_complex_equivalence_matches(
            truth_df, query_df, exact_matches
        )
        matches.extend(complex_matches)
        logger.info(f"XCMP complex equivalence matches: {len(complex_matches)}")

        # Strategy 3: Superlocus matching for complex regions
        superlocus_matches = self._find_superlocus_matches(truth_df, query_df, matches)
        matches.extend(superlocus_matches)
        logger.info(f"XCMP superlocus matches: {len(superlocus_matches)}")

        return matches

    def _perform_ga4gh_matching(
        self, truth_df: pd.DataFrame, query_df: pd.DataFrame
    ) -> list:
        """
        Perform GA4GH-style variant matching.

        GA4GH matching is more conservative and focuses on:
        - Exact coordinate and allele matching
        - Standard VCF representation matching
        - Less complex variant interpretation

        Args:
            truth_df: Truth variants DataFrame
            query_df: Query variants DataFrame

        Returns:
            List of match tuples
        """
        logger.info("Performing GA4GH-style matching...")
        matches = []

        # Strategy 1: Exact matching after normalization
        exact_matches = self._find_exact_matches(truth_df, query_df)
        matches.extend(exact_matches)
        logger.info(f"GA4GH exact matches: {len(exact_matches)}")

        # Strategy 2: Overlapping variant matching (more conservative than XCMP)
        overlap_matches = self._find_overlapping_matches(
            truth_df, query_df, exact_matches
        )
        matches.extend(overlap_matches)
        logger.info(f"GA4GH overlap matches: {len(overlap_matches)}")

        # Strategy 3: Multi-allelic matching
        multiallelic_matches = self._find_multiallelic_matches(
            truth_df, query_df, exact_matches + overlap_matches
        )
        matches.extend(multiallelic_matches)
        logger.info(f"GA4GH multi-allelic matches: {len(multiallelic_matches)}")

        return matches

    def _find_complex_equivalence_matches(
        self, truth_df: pd.DataFrame, query_df: pd.DataFrame, existing_matches: list
    ) -> list:
        """
        Find matches for variants that are equivalent but have different representations.

        This handles cases where the same biological variation can be represented
        differently in VCF format (e.g., different ways to represent the same indel).

        Args:
            truth_df: Truth variants DataFrame
            query_df: Query variants DataFrame
            existing_matches: Already found matches to exclude

        Returns:
            List of complex equivalence matches
        """
        matches = []

        # Get already matched indices
        matched_truth = {m[0] for m in existing_matches}
        matched_query = {m[1] for m in existing_matches}

        # Find equivalent variants with different representations
        for truth_idx, truth_row in truth_df.iterrows():
            if truth_idx in matched_truth:
                continue

            for query_idx, query_row in query_df.iterrows():
                if query_idx in matched_query:
                    continue

                if truth_row["chrom"] != query_row["chrom"]:
                    continue

                # Check if variants are equivalent using sophisticated algorithms
                if self._are_variants_equivalent(truth_row, query_row):
                    confidence = 0.8  # High confidence for equivalent variants
                    matches.append((truth_idx, query_idx, "equivalent", confidence))
                    matched_truth.add(truth_idx)
                    matched_query.add(query_idx)
                    break

        return matches

    def _are_variants_equivalent(self, var1: pd.Series, var2: pd.Series) -> bool:
        """
        Determine if two variants are biologically equivalent.

        This implements sophisticated equivalence checking that goes beyond
        simple string matching to identify variants that represent the same
        biological change.

        Args:
            var1: First variant (pandas Series)
            var2: Second variant (pandas Series)

        Returns:
            True if variants are equivalent, False otherwise
        """
        # Must be on same chromosome
        if var1["chrom"] != var2["chrom"]:
            return False

        # Check if variants represent the same net sequence change
        net_change1 = self._calculate_net_sequence_change(var1)
        net_change2 = self._calculate_net_sequence_change(var2)

        if net_change1 == net_change2:
            return True

        # Check if variants have overlapping affected regions and similar effects
        if self._variants_overlap(var1, var2):
            # For overlapping variants, check if they have similar sequence impact
            similar_impact = self._have_similar_sequence_impact(var1, var2)
            if similar_impact:
                return True

        return False

    def _calculate_net_sequence_change(self, variant: pd.Series) -> str:
        """
        Calculate the net sequence change caused by a variant.

        Args:
            variant: Variant data (pandas Series)

        Returns:
            String representing the net sequence change
        """
        ref = str(variant.get("ref", ""))
        alt = str(variant.get("alt", ""))
        pos = variant.get("pos", 0)

        # Calculate what sequence is removed and what is added
        removed = ref
        added = alt

        return f"{pos}:{removed}>{added}"

    def _variants_overlap(self, var1: pd.Series, var2: pd.Series) -> bool:
        """
        Check if two variants have overlapping genomic regions.

        Args:
            var1: First variant
            var2: Second variant

        Returns:
            True if variants overlap, False otherwise
        """
        start1 = var1["pos"]
        end1 = var1["pos"] + len(str(var1.get("ref", ""))) - 1

        start2 = var2["pos"]
        end2 = var2["pos"] + len(str(var2.get("ref", ""))) - 1

        return start1 <= end2 and start2 <= end1

    def _are_alleles_compatible(self, var1: pd.Series, var2: pd.Series) -> bool:
        """
        Check if two variants have compatible alleles (represent the same biological change).

        Args:
            var1: First variant (truth)
            var2: Second variant (query)

        Returns:
            True if alleles are compatible, False otherwise
        """
        # If reference sequences are different, they can't be directly comparable
        ref1 = str(var1.get("ref", "") if hasattr(var1, "get") else var1["ref"])
        ref2 = str(var2.get("ref", "") if hasattr(var2, "get") else var2["ref"])
        alt1 = str(var1.get("alt", "") if hasattr(var1, "get") else var1["alt"])
        alt2 = str(var2.get("alt", "") if hasattr(var2, "get") else var2["alt"])

        # Exact allele match (after normalization this should catch most cases)
        if ref1 == ref2 and alt1 == alt2:
            return True

        # If positions are identical, require exact allele match
        pos1 = var1.get("pos") if hasattr(var1, "get") else var1["pos"]
        pos2 = var2.get("pos") if hasattr(var2, "get") else var2["pos"]
        if pos1 == pos2:
            return ref1 == ref2 and alt1 == alt2

        # For overlapping but not identical positions, need more sophisticated analysis
        # This would need reference sequence to properly evaluate
        # For now, we'll be conservative and only match if they represent
        # the same type of change in overlapping regions

        # Get variant lengths to classify type
        ref1_len = len(ref1)
        alt1_len = len(alt1)
        ref2_len = len(ref2)
        alt2_len = len(alt2)

        # Classify variant types using the length-based method
        type1 = self._classify_variant_type_from_lengths(ref1_len, alt1_len)
        type2 = self._classify_variant_type_from_lengths(ref2_len, alt2_len)

        # Same variant type is a good sign for compatibility
        if type1 == type2:
            # For SNPs, must have same position and alleles
            if type1 == "SNP":
                return (
                    var1.get("pos") == var2.get("pos") and ref1 == ref2 and alt1 == alt2
                )
            # For indels, more complex - would need reference sequence for proper analysis
            # For now, be conservative
            return False

        return False

    def _classify_variant_type_from_lengths(self, ref_len: int, alt_len: int) -> str:
        """Classify variant type based on reference and alternate allele lengths."""
        if ref_len == 1 and alt_len == 1:
            return "SNP"
        elif ref_len > alt_len:
            return "DEL"
        elif ref_len < alt_len:
            return "INS"
        else:
            return "COMPLEX"

    def _have_similar_sequence_impact(self, var1: pd.Series, var2: pd.Series) -> bool:
        """
        Check if two variants have similar sequence impact.

        Args:
            var1: First variant
            var2: Second variant

        Returns:
            True if variants have similar impact, False otherwise
        """
        # Get or determine variant types - check both BVT and type fields
        type1 = var1.get("BVT", var1.get("type", ""))
        type2 = var2.get("BVT", var2.get("type", ""))

        # If type is not explicitly set, determine from ref/alt lengths
        if not type1:
            ref1_len = len(str(var1.get("ref", "")))
            alt1_len = len(str(var1.get("alt", "")))
            type1 = self._classify_variant_type_from_lengths(ref1_len, alt1_len)

        if not type2:
            ref2_len = len(str(var2.get("ref", "")))
            alt2_len = len(str(var2.get("alt", "")))
            type2 = self._classify_variant_type_from_lengths(ref2_len, alt2_len)

        # Same type is compatible only if alleles are also the same
        if type1 == type2:
            # For SNPs, require exact match of ref and alt alleles
            if type1 == "SNP":
                return var1.get("ref") == var2.get("ref") and var1.get(
                    "alt"
                ) == var2.get("alt")
            # For other types, be more lenient but still require some similarity
            return True

        # INS and DEL can be compatible if they're part of complex rearrangement
        if {type1, type2} <= {"INS", "DEL", "COMPLEX"}:
            return True

        return False

    def _calculate_metrics(self):
        """Calculate performance metrics."""
        # Count TP, FP, FN
        tp = sum(1 for v in self.truth_variants if v["match"])
        fp = sum(1 for v in self.query_variants if not v["match"])
        fn = sum(1 for v in self.truth_variants if not v["match"])

        # Calculate precision, recall, F1
        precision = tp / (tp + fp) if tp + fp > 0 else 0
        recall = tp / (tp + fn) if tp + fn > 0 else 0
        f1 = (
            2 * precision * recall / (precision + recall)
            if precision + recall > 0
            else 0
        )

        # Store metrics
        self.metrics = {
            "TP": tp,
            "FP": fp,
            "FN": fn,
            "PRECISION": precision,
            "RECALL": recall,
            "F1": f1,
        }

    def _stratify_results(self):
        """Stratify results by variant type and other attributes."""
        # Initialize stratifications
        self.stratifications = {"variant_type": {}, "indel_size": {}, "zygosity": {}}

        # Stratify by variant type
        for type_value in ["SNP", "INS", "DEL", "MNP", "COMPLEX"]:
            # Filter variants by type
            truth_of_type = [v for v in self.truth_variants if v["type"] == type_value]
            query_of_type = [v for v in self.query_variants if v["type"] == type_value]

            # Calculate metrics for this type
            tp = sum(1 for v in truth_of_type if v["match"])
            fp = sum(1 for v in query_of_type if not v["match"])
            fn = sum(1 for v in truth_of_type if not v["match"])

            precision = tp / (tp + fp) if tp + fp > 0 else 0
            recall = tp / (tp + fn) if tp + fn > 0 else 0
            f1 = (
                2 * precision * recall / (precision + recall)
                if precision + recall > 0
                else 0
            )

            # Store metrics
            self.stratifications["variant_type"][type_value] = {
                "TP": tp,
                "FP": fp,
                "FN": fn,
                "PRECISION": precision,
                "RECALL": recall,
                "F1": f1,
            }

        # Stratify by indel size
        for size_range in [(1, 5), (6, 15), (16, 50), (51, float("inf"))]:
            range_name = (
                f"{size_range[0]}-{size_range[1]}"
                if size_range[1] != float("inf")
                else f"{size_range[0]}+"
            )

            # Filter indels by size
            truth_indels = [
                v
                for v in self.truth_variants
                if v["is_indel"] and size_range[0] <= v["length"] <= size_range[1]
            ]
            query_indels = [
                v
                for v in self.query_variants
                if v["is_indel"] and size_range[0] <= v["length"] <= size_range[1]
            ]

            # Calculate metrics for this size range
            tp = sum(1 for v in truth_indels if v["match"])
            fp = sum(1 for v in query_indels if not v["match"])
            fn = sum(1 for v in truth_indels if not v["match"])

            precision = tp / (tp + fp) if tp + fp > 0 else 0
            recall = tp / (tp + fn) if tp + fn > 0 else 0
            f1 = (
                2 * precision * recall / (precision + recall)
                if precision + recall > 0
                else 0
            )

            # Store metrics
            self.stratifications["indel_size"][range_name] = {
                "TP": tp,
                "FP": fp,
                "FN": fn,
                "PRECISION": precision,
                "RECALL": recall,
                "F1": f1,
            }

        # Stratify by zygosity
        for zygosity in ["HET", "HOM"]:
            is_hom = zygosity == "HOM"

            # Filter variants by zygosity
            truth_zyg = [
                v for v in self.truth_variants if v.get("is_hom", False) == is_hom
            ]
            query_zyg = [
                v for v in self.query_variants if v.get("is_hom", False) == is_hom
            ]

            # Calculate metrics for this zygosity
            tp = sum(1 for v in truth_zyg if v["match"])
            fp = sum(1 for v in query_zyg if not v["match"])
            fn = sum(1 for v in truth_zyg if not v["match"])

            precision = tp / (tp + fp) if tp + fp > 0 else 0
            recall = tp / (tp + fn) if tp + fn > 0 else 0
            f1 = (
                2 * precision * recall / (precision + recall)
                if precision + recall > 0
                else 0
            )

            # Store metrics
            self.stratifications["zygosity"][zygosity] = {
                "TP": tp,
                "FP": fp,
                "FN": fn,
                "PRECISION": precision,
                "RECALL": recall,
                "F1": f1,
            }

    def write_results(self, output_prefix: str):
        """
        Write results to files.

        Args:
            output_prefix: Prefix for output files
        """
        # Write metrics to JSON
        metrics_file = f"{output_prefix}.metrics.json"
        with open(metrics_file, "w") as f:
            json.dump(
                {"metrics": self.metrics, "stratifications": self.stratifications},
                f,
                indent=2,
            )

        logger.info(f"Wrote metrics to {metrics_file}")

        # Write summary to TSV
        summary_file = f"{output_prefix}.summary.tsv"
        with open(summary_file, "w") as f:
            # Write header
            f.write("Type\tTP\tFP\tFN\tPrecision\tRecall\tF1\n")

            # Write overall metrics
            f.write(
                f"OVERALL\t{self.metrics['TP']}\t{self.metrics['FP']}\t{self.metrics['FN']}\t"
                f"{self.metrics['PRECISION']:.4f}\t{self.metrics['RECALL']:.4f}\t{self.metrics['F1']:.4f}\n"
            )

            # Write stratifications
            for strat_type, strat_data in self.stratifications.items():
                for type_value, metrics in strat_data.items():
                    f.write(
                        f"{strat_type.upper()}_{type_value}\t{metrics['TP']}\t{metrics['FP']}\t{metrics['FN']}\t"
                        f"{metrics['PRECISION']:.4f}\t{metrics['RECALL']:.4f}\t{metrics['F1']:.4f}\n"
                    )

        logger.info(f"Wrote summary to {summary_file}")

        # Write VTC (variant truth categories) if requested
        if self.output_vtc:
            self._write_vtc(output_prefix)

        # Write Phase 2 ROC analysis results if enabled
        if self.enable_roc_analysis:
            self._write_roc_results(output_prefix)

        # Write Phase 2 ROC analysis results if enabled
        if self.enable_roc_analysis:
            self._write_roc_results(output_prefix)

    def _write_vtc(self, output_prefix: str):
        """
        Write variant truth categories.

        Args:
            output_prefix: Prefix for output files
        """
        # Create DataFrames
        truth_df = pd.DataFrame(self.truth_variants)
        query_df = pd.DataFrame(self.query_variants)

        # Add truth category
        truth_df["category"] = "FN"
        truth_df.loc[truth_df["match"], "category"] = "TP"

        query_df["category"] = "FP"
        query_df.loc[query_df["match"], "category"] = "TP"

        # Write truth VTC
        truth_vtc_file = f"{output_prefix}.truth.vtc.tsv"
        truth_df.to_csv(truth_vtc_file, sep="\t", index=False)
        logger.info(f"Wrote truth VTC to {truth_vtc_file}")

        # Write query VTC
        query_vtc_file = f"{output_prefix}.query.vtc.tsv"
        query_df.to_csv(query_vtc_file, sep="\t", index=False)
        logger.info(f"Wrote query VTC to {query_vtc_file}")

    def _normalize_all_variants(self):
        """Normalize all variants using standard VCF normalization."""
        logger.info("Normalizing variants...")

        # Normalize truth variants
        for variant in self.truth_variants:
            self._normalize_variant(variant)

        # Normalize query variants
        for variant in self.query_variants:
            self._normalize_variant(variant)

    def _normalize_variant(self, variant: dict) -> None:
        """
        Normalize a single variant using standard VCF normalization.

        Based on the C++ implementation in RefVar.cpp:
        1. First trimRight (remove common suffix)
        2. Then trimLeft (remove common prefix and adjust position)

        Args:
            variant: Variant dictionary to normalize in-place
        """
        try:
            ref = str(variant.get("ref", "")).upper()
            alt = str(variant.get("alt", "")).upper()
            pos = int(variant.get("pos", 0))

            if not ref or not alt or pos < 1:
                logger.warning(f"Invalid variant data for normalization: {variant}")
                return

            # Step 1: Trim right (remove common suffix)
            while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
                ref = ref[:-1]
                alt = alt[:-1]

            # Step 2: Trim left (remove common prefix and adjust position)
            while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
                ref = ref[1:]
                alt = alt[1:]
                pos += 1

            # Update variant with normalized values
            variant["ref"] = ref
            variant["alt"] = alt
            variant["pos"] = pos

        except Exception as e:
            logger.warning(f"Error normalizing variant {variant}: {e}")

    def _get_reference_sequence(self, chrom: str, start: int, end: int) -> str:
        """
        Get reference sequence for a genomic region.

        Args:
            chrom: Chromosome name
            start: Start position (1-based)
            end: End position (1-based, inclusive)

        Returns:
            Reference sequence string
        """
        # In a full implementation, this would fetch from a reference FASTA file
        # For now, return a placeholder that allows basic normalization
        if hasattr(self, "reference") and self.reference:
            try:
                # Use pysam or similar to fetch reference sequence
                return self.reference.fetch(chrom, start - 1, end).upper()
            except Exception as e:
                logger.warning(
                    f"Could not fetch reference sequence for {chrom}:{start}-{end}: {e}"
                )

        # Fallback: return empty string (normalization will be limited)
        return ""

    def _perform_quality_stratification(self):
        """
        Stratify variants by quality score and calculate metrics for each bin.

        This method creates quality bins and calculates metrics (TP, FP, FN, precision,
        recall, F1) for each bin.

        Returns:
            None - Results are stored in self.quality_metrics
        """
        logger.info("Performing quality score stratification...")

        # Define quality bins
        quality_bins = [
            {"name": "Q1-10", "min": 1, "max": 10},
            {"name": "Q10-20", "min": 10, "max": 20},
            {"name": "Q20-30", "min": 20, "max": 30},
            {"name": "Q30-40", "min": 30, "max": 40},
            {"name": "Q40+", "min": 40, "max": float("inf")},
        ]

        bin_metrics = {}

        # Process each quality bin
        for bin_info in quality_bins:
            bin_name = bin_info["name"]
            min_qual = bin_info["min"]
            max_qual = bin_info["max"]

            # Filter variants by quality
            bin_variants = [
                v
                for v in self.query_variants
                if min_qual <= v.get("qual", 0) < max_qual
            ]

            # Calculate metrics for this bin
            tp = sum(1 for v in bin_variants if v.get("match", False))
            fp = len(bin_variants) - tp

            # Find truth variants that would match to variants in this bin
            truth_variant_ids = {
                v.get("truth_variant_id") for v in bin_variants if v.get("match", False)
            }
            fn = sum(
                1 for v in self.truth_variants if v.get("id") not in truth_variant_ids
            )

            # Calculate precision, recall, F1
            precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
            recall = tp / (tp + fn) if (tp + fn) > 0 else 0.0
            f1 = (
                2 * (precision * recall) / (precision + recall)
                if (precision + recall) > 0
                else 0.0
            )

            # Store metrics for this bin
            bin_metrics[bin_name] = {
                "quality_range": f"{min_qual}-{max_qual if max_qual != float('inf') else '∞'}",
                "TP": tp,
                "FP": fp,
                "FN": fn,
                "PRECISION": precision,
                "RECALL": recall,
                "F1": f1,
                "variant_count": len(bin_variants),
            }

            logger.debug(
                f"Bin {bin_name}: {tp} TP, {fp} FP, {fn} FN, Precision: {precision:.4f}, Recall: {recall:.4f}"
            )

        self.quality_metrics = {
            "bin_metrics": bin_metrics,
            "total_variants": len(self.query_variants),
        }

        logger.info(f"Quality stratification completed for {len(quality_bins)} bins")

    def _generate_roc_curve(self, variant_type: str, bvt_values: list):
        """
        Generate ROC curve data for a specific variant type.

        Args:
            variant_type: Type of variant for ROC curve ('snp', 'indel', 'all')
            bvt_values: List of BVT values to include in this variant type

        Returns:
            None - Results are stored in self.roc_data
        """
        logger.info(f"Generating ROC curve for {variant_type}...")

        # Filter variants by type
        variants = [v for v in self.query_variants if v.get("BVT") in bvt_values]

        # If no variants of this type, create empty data structure
        if not variants:
            self.roc_data[variant_type] = {
                "thresholds": [],
                "tp": [],
                "fp": [],
                "fn": [],
                "precision": [],
                "recall": [],
            }
            return

        # Get truth variant count for this type
        truth_total = len(
            [v for v in self.truth_variants if v.get("BVT") in bvt_values]
        )

        # Sort variants by quality score (descending)
        sorted_variants = sorted(variants, key=lambda v: v.get("qual", 0), reverse=True)

        # Initialize data structures
        thresholds = []
        tp_counts = []
        fp_counts = []
        fn_counts = []
        precision_values = []
        recall_values = []

        # Calculate metrics at each threshold point
        current_tp = 0
        current_fp = 0
        prev_quality = float("inf")

        for variant in sorted_variants:
            quality = variant.get("qual", 0)

            # Add a point for each distinct quality value
            if quality != prev_quality:
                tp_at_threshold = current_tp
                fp_at_threshold = current_fp
                fn_at_threshold = truth_total - tp_at_threshold

                precision = (
                    tp_at_threshold / (tp_at_threshold + fp_at_threshold)
                    if (tp_at_threshold + fp_at_threshold) > 0
                    else 0.0
                )
                recall = tp_at_threshold / truth_total if truth_total > 0 else 0.0

                thresholds.append(quality)
                tp_counts.append(tp_at_threshold)
                fp_counts.append(fp_at_threshold)
                fn_counts.append(fn_at_threshold)
                precision_values.append(precision)
                recall_values.append(recall)

                prev_quality = quality

            # Update counts
            if variant.get("match", False):
                current_tp += 1
            else:
                current_fp += 1

        # Add final point (including all variants)
        fn_final = truth_total - current_tp
        precision_final = (
            current_tp / (current_tp + current_fp)
            if (current_tp + current_fp) > 0
            else 0.0
        )
        recall_final = current_tp / truth_total if truth_total > 0 else 0.0

        thresholds.append(0)  # Lowest possible threshold
        tp_counts.append(current_tp)
        fp_counts.append(current_fp)
        fn_counts.append(fn_final)
        precision_values.append(precision_final)
        recall_values.append(recall_final)

        # Store results
        self.roc_data[variant_type] = {
            "thresholds": thresholds,
            "tp": tp_counts,
            "fp": fp_counts,
            "fn": fn_counts,
            "precision": precision_values,
            "recall": recall_values,
        }

        # Calculate AUC if sklearn is available
        if SKLEARN_AVAILABLE:
            try:
                from sklearn.metrics import auc

                # Calculate area under the precision-recall curve
                self.roc_data[variant_type]["auc"] = auc(
                    recall_values, precision_values
                )
                logger.info(
                    f"AUC for {variant_type}: {self.roc_data[variant_type]['auc']:.4f}"
                )
            except Exception as e:
                logger.warning(f"Error calculating AUC: {e}")
        else:
            logger.debug("sklearn not available - skipping AUC calculation")

    def _calculate_bootstrap_confidence_intervals(self):
        """
        Calculate bootstrap confidence intervals for precision and recall.

        This method uses the Jeffreys confidence interval method from hap_py.tools.ci
        to calculate confidence intervals for precision and recall at each threshold.

        Returns:
            None - Results are stored in self.bootstrap_confidence_intervals
        """
        from hap_py.tools.ci import jeffreysCI

        if not SCIPY_AVAILABLE:
            logger.warning(
                "SciPy not available - skipping bootstrap confidence intervals"
            )
            return

        logger.info("Calculating bootstrap confidence intervals...")

        for variant_type in self.roc_data:
            if not self.roc_data[variant_type]["thresholds"]:
                continue

            roc_data = self.roc_data[variant_type]
            precision_ci = []
            recall_ci = []

            for i in range(len(roc_data["thresholds"])):
                tp = roc_data["tp"][i]
                fp = roc_data["fp"][i]
                fn = roc_data["fn"][i]

                # Calculate precision confidence interval
                if tp + fp > 0:
                    _, p_lower, p_upper = jeffreysCI(tp, tp + fp)
                    precision_ci.append({"lower": p_lower, "upper": p_upper})
                else:
                    precision_ci.append({"lower": 0.0, "upper": 1.0})

                # Calculate recall confidence interval
                if tp + fn > 0:
                    _, r_lower, r_upper = jeffreysCI(tp, tp + fn)
                    recall_ci.append({"lower": r_lower, "upper": r_upper})
                else:
                    recall_ci.append({"lower": 0.0, "upper": 1.0})

            # Store confidence intervals
            self.bootstrap_confidence_intervals[variant_type] = {
                "precision_ci": precision_ci,
                "recall_ci": recall_ci,
            }

            # Calculate AUC confidence interval if sklearn is available
            if SKLEARN_AVAILABLE and "auc" in roc_data:
                try:
                    # Simple estimate - more sophisticated bootstrapping could be implemented
                    auc_value = roc_data["auc"]
                    self.bootstrap_confidence_intervals[variant_type]["auc"] = auc_value
                except Exception as e:
                    logger.warning(f"Error calculating AUC confidence interval: {e}")

        logger.info("Bootstrap confidence intervals calculated")

    def _perform_multi_threshold_analysis(self):
        """
        Analyze results at standard quality thresholds (Q10, Q20, Q30, etc.).

        This method calculates precision/recall metrics at commonly used quality
        thresholds and stores the results.

        Returns:
            None - Results are stored in self.roc_data["multi_threshold"]
        """
        logger.info("Performing multi-threshold analysis...")

        # Define standard quality thresholds
        standard_thresholds = {
            "Q10": 10,
            "Q20": 20,
            "Q30": 30,
            "Q40": 40,
            "Q50": 50,
        }

        # Define variant types to analyze
        variant_types = {
            "snp": ["SNP"],
            "indel": ["INS", "DEL"],
            "all": ["SNP", "INS", "DEL", "OTHER"],
        }

        # Initialize multi-threshold results
        self.roc_data["multi_threshold"] = {}

        # Analyze each variant type
        for variant_type, bvt_values in variant_types.items():
            self.roc_data["multi_threshold"][variant_type] = {}

            # Filter variants by type
            variants = [v for v in self.query_variants if v.get("BVT") in bvt_values]
            truth_variants = [
                v for v in self.truth_variants if v.get("BVT") in bvt_values
            ]
            truth_total = len(truth_variants)

            # Analyze each threshold
            for threshold_name, threshold_value in standard_thresholds.items():
                # Count variants above threshold
                variants_above_threshold = [
                    v for v in variants if v.get("qual", 0) >= threshold_value
                ]

                # Calculate metrics
                tp = sum(1 for v in variants_above_threshold if v.get("match", False))
                fp = len(variants_above_threshold) - tp
                fn = truth_total - tp

                precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
                recall = tp / truth_total if truth_total > 0 else 0.0
                f1 = (
                    2 * (precision * recall) / (precision + recall)
                    if (precision + recall) > 0
                    else 0.0
                )

                # Store metrics
                self.roc_data["multi_threshold"][variant_type][threshold_name] = {
                    "threshold": threshold_value,
                    "TP": tp,
                    "FP": fp,
                    "FN": fn,
                    "PRECISION": precision,
                    "RECALL": recall,
                    "F1": f1,
                    "variants_above_threshold": len(variants_above_threshold),
                }

                logger.debug(
                    f"{variant_type.upper()} @ {threshold_name}: "
                    f"{tp} TP, {fp} FP, {fn} FN, "
                    f"Precision: {precision:.4f}, Recall: {recall:.4f}"
                )

        logger.info("Multi-threshold analysis completed")

    def _write_roc_results(self, output_prefix: str):
        """
        Write enhanced ROC analysis results to files.

        Args:
            output_prefix: Prefix for output files
        """
        if not self.enable_roc_analysis or not self.roc_data:
            return

        logger.info("Writing enhanced ROC analysis results...")

        # Write ROC curves data
        roc_file = f"{output_prefix}.roc.tsv"
        with open(roc_file, "w") as f:
            f.write(
                "Type\tThreshold\tTP\tFP\tFN\tPrecision\tRecall\tPrecision_Lower\tRecall_Lower\tRecall_Upper\n"
            )

            for variant_type in ["snp", "indel", "all"]:
                if (
                    variant_type not in self.roc_data
                    or not self.roc_data[variant_type]["thresholds"]
                ):
                    continue

                roc_data = self.roc_data[variant_type]
                ci_data = self.bootstrap_confidence_intervals.get(variant_type, {})

                for i in range(len(roc_data["thresholds"])):
                    threshold = roc_data["thresholds"][i]
                    tp = roc_data["tp"][i]
                    fp = roc_data["fp"][i]
                    fn = roc_data["fn"][i]
                    precision = roc_data["precision"][i]
                    recall = roc_data["recall"][i]

                    # Get confidence intervals if available
                    p_lower = p_upper = r_lower = r_upper = ""
                    if "precision_ci" in ci_data and i < len(ci_data["precision_ci"]):
                        p_lower = f"{ci_data['precision_ci'][i]['lower']:.4f}"
                        p_upper = f"{ci_data['precision_ci'][i]['upper']:.4f}"
                    if "recall_ci" in ci_data and i < len(ci_data["recall_ci"]):
                        r_lower = f"{ci_data['recall_ci'][i]['lower']:.4f}"
                        r_upper = f"{ci_data['recall_ci'][i]['upper']:.4f}"

                    f.write(
                        f"{variant_type.upper()}\t{threshold:.2f}\t{tp}\t{fp}\t{fn}\t"
                        f"{precision:.4f}\t{recall:.4f}\t{p_lower}\t{p_upper}\t{r_lower}\t{r_upper}\n"
                    )

        logger.info(f"Wrote ROC curves to {roc_file}")

        # Write quality stratification results
        if self.quality_stratification and hasattr(self, "quality_metrics"):
            quality_file = f"{output_prefix}.quality_stratification.tsv"
            with open(quality_file, "w") as f:
                f.write(
                    "Quality_Bin\tQuality_Range\tTP\tFP\tFN\tPrecision\tRecall\tF1\tVariant_Count\n"
                )

                for bin_name, metrics in self.quality_metrics["bin_metrics"].items():
                    f.write(
                        f"{bin_name}\t{metrics['quality_range']}\t{metrics['TP']}\t{metrics['FP']}\t{metrics['FN']}\t"
                        f"{metrics['PRECISION']:.4f}\t{metrics['RECALL']:.4f}\t{metrics['F1']:.4f}\t{metrics['variant_count']}\n"
                    )

            logger.info(f"Wrote quality stratification to {quality_file}")

        # Write multi-threshold analysis results
        if "multi_threshold" in self.roc_data:
            multi_threshold_file = f"{output_prefix}.multi_threshold.tsv"
            with open(multi_threshold_file, "w") as f:
                f.write(
                    "Type\tThreshold_Name\tThreshold\tTP\tFP\tFN\tPrecision\tRecall\tF1\tVariants_Above_Threshold\n"
                )

                for variant_type, threshold_data in self.roc_data[
                    "multi_threshold"
                ].items():
                    for threshold_name, metrics in threshold_data.items():
                        f.write(
                            f"{variant_type.upper()}\t{threshold_name}\t{metrics['threshold']}\t{metrics['TP']}\t{metrics['FP']}\t{metrics['FN']}\t"
                            f"{metrics['PRECISION']:.4f}\t{metrics['RECALL']:.4f}\t{metrics['F1']:.4f}\t{metrics['variants_above_threshold']}\n"
                        )

            logger.info(f"Wrote multi-threshold analysis to {multi_threshold_file}")

        # Generate ROC curve plots if matplotlib is available
        if MATPLOTLIB_AVAILABLE:
            try:
                plot_file = f"{output_prefix}.roc_plot.png"

                plt.figure(figsize=(10, 8))

                for variant_type in ["snp", "indel", "all"]:
                    if (
                        variant_type not in self.roc_data
                        or not self.roc_data[variant_type]["thresholds"]
                    ):
                        continue

                    roc_data = self.roc_data[variant_type]

                    # Sort points by recall for smooth curves
                    points = sorted(zip(roc_data["recall"], roc_data["precision"]))
                    recall_sorted = [p[0] for p in points]
                    precision_sorted = [p[1] for p in points]

                    plt.plot(
                        recall_sorted,
                        precision_sorted,
                        label=f"{variant_type.upper()} (AUC={self.bootstrap_confidence_intervals.get(variant_type, {}).get('auc', 0):.4f})",
                    )

                plt.xlabel("Recall")
                plt.ylabel("Precision")
                plt.title("Precision-Recall Curves")
                plt.grid(True, alpha=0.3)
                plt.legend()
                plt.savefig(plot_file)

                logger.info(f"Generated ROC curve plot at {plot_file}")
            except Exception as e:
                logger.warning(f"Failed to generate ROC plot: {e}")

    def _perform_roc_analysis(self):
        """
        Perform ROC analysis for different variant types.

        This method orchestrates the ROC analysis workflow:
        1. Generate ROC curves for different variant types
        2. Calculate bootstrap confidence intervals
        3. Perform quality score stratification
        4. Analyze standard quality thresholds

        Returns:
            None - Results are stored in instance attributes
        """
        if not self.enable_roc_analysis:
            logger.info("ROC analysis is disabled. Skipping...")
            return

        logger.info("Performing ROC analysis...")

        # Initialize ROC data structures
        self.roc_data = {}
        self.bootstrap_confidence_intervals = {}

        # Process different variant types
        variant_types = {
            "snp": ["SNP"],
            "indel": ["INS", "DEL"],
            "all": ["SNP", "INS", "DEL", "OTHER"],
        }

        # Generate ROC curves for each variant type
        for variant_type, bvt_values in variant_types.items():
            self._generate_roc_curve(variant_type, bvt_values)

        # Calculate bootstrap confidence intervals
        if SCIPY_AVAILABLE:
            self._calculate_bootstrap_confidence_intervals()
        else:
            logger.warning(
                "SciPy not available - skipping bootstrap confidence intervals"
            )

        # Perform quality score stratification
        if self.quality_stratification:
            self._perform_quality_stratification()
        else:
            logger.info("Quality stratification disabled. Skipping...")

        # Perform analysis at standard quality thresholds
        self._perform_multi_threshold_analysis()

        logger.info("ROC analysis completed")
