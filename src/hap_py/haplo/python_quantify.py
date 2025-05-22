#!/usr/bin/env python3
"""
Python implementation of the quantify module.

This module provides functionality to quantify variants in VCF files,
producing stratification metrics and summary statistics.
"""

import json
import logging
from typing import Optional

import pandas as pd
import pysam

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
        """
        self.truth_vcf = truth_vcf
        self.query_vcf = query_vcf
        self.reference = reference
        self.regions = regions
        self.apply_filters = apply_filters
        self.output_vtc = output_vtc

        self.truth_variants = []
        self.query_variants = []
        self.region_list = []

        # Results storage
        self.metrics = {}
        self.stratifications = {}

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

        # Open reference if provided
        if self.reference:
            try:
                self.reference_handle = pysam.FastaFile(self.reference)
                logger.info(f"Opened reference: {self.reference}")
            except Exception as e:
                logger.warning(f"Failed to open reference: {e}")
                self.reference = None

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

        # Return results
        return {"metrics": self.metrics, "stratifications": self.stratifications}

    def _match_variants(self):
        """Match variants between truth and query sets."""
        # Convert to pandas DataFrames for easier matching
        truth_df = pd.DataFrame(self.truth_variants)
        query_df = pd.DataFrame(self.query_variants)

        # Add 'source' column
        truth_df["source"] = "truth"
        query_df["source"] = "query"

        # Combine dataframes
        all_variants = pd.concat([truth_df, query_df])

        # Group by chromosome, position, reference, and alternate alleles
        # This identifies matching variants
        grouped = all_variants.groupby(["chrom", "pos", "ref", "alt"])

        # Initialize match columns
        truth_df["match"] = False
        query_df["match"] = False
        truth_df["match_idx"] = -1
        query_df["match_idx"] = -1

        # Mark matches
        for (chrom, pos, ref, alt), group in grouped:
            if (
                len(group) > 1
                and "truth" in group["source"].values
                and "query" in group["source"].values
            ):
                # Get truth and query indices
                truth_idx = group[group["source"] == "truth"].index
                query_idx = group[group["source"] == "query"].index

                # Mark as matched
                truth_df.loc[truth_idx, "match"] = True
                query_df.loc[query_idx, "match"] = True

                # Store match indices
                for t_idx in truth_idx:
                    truth_df.loc[t_idx, "match_idx"] = query_idx[0]

                for q_idx in query_idx:
                    query_df.loc[q_idx, "match_idx"] = truth_idx[0]

        # Update variant lists with match information
        self.truth_variants = truth_df.to_dict("records")
        self.query_variants = query_df.to_dict("records")

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


def main():
    """Command-line entry point."""
    import argparse

    parser = argparse.ArgumentParser(description="Quantify variant calls in VCF files")

    parser.add_argument("--truth-vcf", required=True, help="Path to truth VCF file")
    parser.add_argument("--query-vcf", required=True, help="Path to query VCF file")
    parser.add_argument("--reference", help="Path to reference FASTA file")
    parser.add_argument("--regions", help="BED file with regions to quantify")
    parser.add_argument(
        "--apply-filters", action="store_true", help="Apply filters from VCF"
    )
    parser.add_argument(
        "--output-vtc", action="store_true", help="Output variant truth categories"
    )
    parser.add_argument(
        "--output-prefix", required=True, help="Prefix for output files"
    )

    args = parser.parse_args()

    # Run quantification
    engine = QuantifyEngine(
        truth_vcf=args.truth_vcf,
        query_vcf=args.query_vcf,
        reference=args.reference,
        regions=args.regions,
        apply_filters=args.apply_filters,
        output_vtc=args.output_vtc,
    )

    results = engine.quantify()
    engine.write_results(args.output_prefix)

    # Print summary
    print("\nQuantification Summary:")
    print(f"Total TP: {results['metrics']['TP']}")
    print(f"Total FP: {results['metrics']['FP']}")
    print(f"Total FN: {results['metrics']['FN']}")
    print(f"Precision: {results['metrics']['PRECISION']:.4f}")
    print(f"Recall: {results['metrics']['RECALL']:.4f}")
    print(f"F1 Score: {results['metrics']['F1']:.4f}")


if __name__ == "__main__":
    main()
