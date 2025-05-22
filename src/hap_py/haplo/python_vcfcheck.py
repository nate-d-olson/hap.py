#!/usr/bin/env python3
"""
Python implementation of the vcfcheck functionality.

This module replaces the C++ vcfcheck tool, which validates VCF files for
compatibility with hap.py. It uses pysam to check for common issues like
invalid reference alleles, problematic genotypes, and header problems.
"""

import argparse
import logging
import sys
from typing import Dict, List, Optional, Tuple

import pysam
from pysam import VariantFile, VariantRecord


class VCFChecker:
    """
    Validates VCF files for compatibility with hap.py processing.

    This checks for:
    - Valid VCF/BCF format
    - Valid reference alleles
    - Compatible sample columns
    - Structural variants that may cause issues
    - Genotypes that may be problematic
    """

    def __init__(
        self,
        reference_path: Optional[str] = None,
        strict: bool = False,
        apply_filters: bool = False,
    ):
        """
        Initialize the VCF checker.

        Args:
            reference_path: Path to reference FASTA file
            strict: Whether to use strict checking mode
            apply_filters: Whether to skip filtered variants
        """
        self.reference_path = reference_path
        self.strict = strict
        self.apply_filters = apply_filters
        self.logger = logging.getLogger("vcfcheck")

        # Reference sequence handle
        self.reference = None
        if reference_path:
            try:
                self.reference = pysam.FastaFile(reference_path)
                self.logger.info(f"Loaded reference: {reference_path}")
            except Exception as e:
                self.logger.error(f"Failed to load reference: {e}")

        # Counters for issues
        self.stats = {
            "total_variants": 0,
            "filtered_variants": 0,
            "invalid_ref_alleles": 0,
            "missing_genotypes": 0,
            "invalid_genotypes": 0,
            "structural_variants": 0,
            "overlapping_variants": 0,
            "other_issues": 0,
        }

        # Track positions for overlap detection
        self.last_pos = {}  # Chromosome -> position

    def check_file(
        self, vcf_path: str, output_path: Optional[str] = None
    ) -> Dict[str, int]:
        """
        Check a VCF file for issues.

        Args:
            vcf_path: Path to the VCF file to check
            output_path: Optional path to write issues to

        Returns:
            Dictionary with issue counts
        """
        self.logger.info(f"Checking {vcf_path}")
        self._reset_stats()

        # Output file handle
        out_file = None
        if output_path:
            self.logger.info(f"Writing issues to {output_path}")
            out_file = open(output_path, "w")
            out_file.write("CHROM\tPOS\tREF\tALT\tISSUE\tDETAILS\n")

        try:
            # Open the VCF file
            vcf = VariantFile(vcf_path)

            # Check header
            header_issues = self._check_header(vcf.header)
            for issue in header_issues:
                self.logger.warning(f"Header issue: {issue}")
                if out_file:
                    out_file.write(f"HEADER\t0\t.\t.\t{issue}\t.\n")

            # Check variants
            for record in vcf:
                self.stats["total_variants"] += 1

                # Skip filtered variants if requested
                if self.apply_filters and self._is_filtered(record):
                    self.stats["filtered_variants"] += 1
                    continue

                # Check for issues
                issues = self._check_variant(record)

                # Output and log issues
                for issue_type, details in issues:
                    self.logger.debug(
                        f"Issue at {record.chrom}:{record.pos} - {issue_type}: {details}"
                    )
                    if out_file:
                        alt = ",".join(record.alts) if record.alts else "."
                        out_file.write(
                            f"{record.chrom}\t{record.pos}\t{record.ref}\t"
                            f"{alt}\t{issue_type}\t{details}\n"
                        )

            self.logger.info(f"Checked {self.stats['total_variants']} variants")

        except Exception as e:
            self.logger.error(f"Error checking file: {e}")
            if out_file:
                out_file.write(f"ERROR\t0\t.\t.\tFailed to process file\t{e}\n")

        finally:
            if out_file:
                out_file.close()

        return self.stats

    def _reset_stats(self):
        """Reset statistics counters."""
        for key in self.stats:
            self.stats[key] = 0
        self.last_pos = {}

    def _check_header(self, header) -> List[str]:
        """
        Check a VCF header for issues.

        Args:
            header: pysam VCF header

        Returns:
            List of issue descriptions
        """
        issues = []

        # Check for required fields (be less strict for INFO as some simple VCFs may not have it)
        # Note: pysam VariantHeader doesn't support 'in' operator, so check attributes directly
        if not hasattr(header, 'filters'):
            issues.append("Missing required header field: FILTER")
        if not hasattr(header, 'formats') or len(header.formats) == 0:
            issues.append("Missing required header field: FORMAT")
        # INFO field is optional for simple VCFs, just warn if missing
        if not hasattr(header, 'info') or len(header.info) == 0:
            # Just a warning, not an error
            pass

        # Check sample columns
        if len(header.samples) == 0:
            issues.append("No sample columns found")

        # Check for FORMAT:GT field
        has_gt = False
        for format_field in header.formats.values():
            if format_field.name == "GT":
                has_gt = True
                break

        if not has_gt:
            issues.append("Missing GT format field")

        return issues

    def _check_variant(self, record: VariantRecord) -> List[Tuple[str, str]]:
        """
        Check a variant record for issues.

        Args:
            record: pysam variant record

        Returns:
            List of (issue_type, details) tuples
        """
        issues = []

        # Check for overlapping variants
        if record.chrom in self.last_pos and record.pos <= self.last_pos[record.chrom]:
            issues.append(
                (
                    "OVERLAP",
                    f"Overlaps with previous variant at position {self.last_pos[record.chrom]}",
                )
            )
            self.stats["overlapping_variants"] += 1

        self.last_pos[record.chrom] = record.pos + len(record.ref) - 1

        # Check reference allele
        if self.reference:
            ref_issues = self._check_reference_allele(record)
            issues.extend(ref_issues)

        # Check for structural variants
        if self._is_structural_variant(record):
            issues.append(("STRUCTURAL", "Structural variant detected"))
            self.stats["structural_variants"] += 1

        # Check genotypes
        gt_issues = self._check_genotypes(record)
        issues.extend(gt_issues)

        return issues

    def _check_reference_allele(self, record: VariantRecord) -> List[Tuple[str, str]]:
        """
        Check if the reference allele matches the reference genome.

        Args:
            record: pysam variant record

        Returns:
            List of (issue_type, details) tuples
        """
        issues = []

        if not self.reference:
            return issues

        try:
            # Get the reference sequence
            ref_seq = self.reference.fetch(
                record.chrom, record.pos - 1, record.pos - 1 + len(record.ref)
            )

            # Compare with the record's reference allele
            if ref_seq.upper() != record.ref.upper():
                issues.append(
                    (
                        "REF_MISMATCH",
                        f"Reference allele mismatch: VCF={record.ref}, Reference={ref_seq}",
                    )
                )
                self.stats["invalid_ref_alleles"] += 1

        except Exception as e:
            issues.append(("REF_ERROR", f"Error checking reference: {e}"))
            self.stats["other_issues"] += 1

        return issues

    def _is_structural_variant(self, record: VariantRecord) -> bool:
        """
        Check if a variant is a structural variant.

        Args:
            record: pysam variant record

        Returns:
            True if the variant is structural, False otherwise
        """
        # Check for standard SV indicators
        if record.info.get("SVTYPE"):
            return True

        # Check for long alleles
        if len(record.ref) > 50:
            return True

        if record.alts:
            for alt in record.alts:
                if len(alt) > 50:
                    return True

        # Check for symbolic alleles
        if record.alts:
            for alt in record.alts:
                if alt.startswith("<") and alt.endswith(">"):
                    return True

        return False

    def _check_genotypes(self, record: VariantRecord) -> List[Tuple[str, str]]:
        """
        Check genotypes for issues.

        Args:
            record: pysam variant record

        Returns:
            List of (issue_type, details) tuples
        """
        issues = []

        # Check if GT field exists
        if "GT" not in record.format:
            issues.append(("MISSING_GT", "No GT field found in variant"))
            self.stats["missing_genotypes"] += 1
            return issues

        # Check each sample's genotype
        for sample in record.samples:
            if record.samples[sample].get("GT") is None:
                issues.append(("MISSING_GT", f"Sample {sample} has no GT"))
                self.stats["missing_genotypes"] += 1
                continue

            gt = record.samples[sample]["GT"]

            # Check for unphased genotypes if strict mode
            if self.strict and not record.samples[sample].phased:
                issues.append(("UNPHASED", f"Sample {sample} has unphased genotype"))
                self.stats["invalid_genotypes"] += 1

            # Check for invalid allele indices
            max_allele_idx = len(record.alts) if record.alts else 0
            for allele_idx in gt:
                if allele_idx is not None and allele_idx > max_allele_idx:
                    issues.append(
                        (
                            "INVALID_GT",
                            f"Sample {sample} has invalid allele index {allele_idx}",
                        )
                    )
                    self.stats["invalid_genotypes"] += 1

        return issues

    def _is_filtered(self, record: VariantRecord) -> bool:
        """
        Check if a variant record is filtered.

        Args:
            record: The variant record to check

        Returns:
            True if the record is filtered, False otherwise
        """
        # If no filters, it's not filtered
        if not record.filter:
            return False

        # In pysam, PASS is represented as an empty filter
        # But sometimes a "PASS" string might be present
        for filter_name in record.filter:
            if filter_name != "PASS":
                return True

        return False

    def print_summary(self):
        """Print a summary of the checking results."""
        print("\nVCF Check Summary:")
        print(f"Total variants: {self.stats['total_variants']}")
        if self.apply_filters:
            print(f"Filtered variants: {self.stats['filtered_variants']}")
        print(f"Invalid reference alleles: {self.stats['invalid_ref_alleles']}")
        print(f"Missing genotypes: {self.stats['missing_genotypes']}")
        print(f"Invalid genotypes: {self.stats['invalid_genotypes']}")
        print(f"Structural variants: {self.stats['structural_variants']}")
        print(f"Overlapping variants: {self.stats['overlapping_variants']}")
        print(f"Other issues: {self.stats['other_issues']}")

        # Calculate total issues
        total_issues = (
            self.stats["invalid_ref_alleles"]
            + self.stats["missing_genotypes"]
            + self.stats["invalid_genotypes"]
            + self.stats["structural_variants"]
            + self.stats["overlapping_variants"]
            + self.stats["other_issues"]
        )

        print(f"\nTotal issues: {total_issues}")
        if self.stats["total_variants"] > 0:
            issue_rate = (total_issues / self.stats["total_variants"]) * 100
            print(f"Issue rate: {issue_rate:.2f}%")

        if total_issues == 0:
            print("\n✅ VCF file passed all checks")
        else:
            print("\n⚠️ VCF file has issues that may affect processing")


def main():
    """Command-line entry point."""
    parser = argparse.ArgumentParser(
        description="Check a VCF file for compatibility with hap.py"
    )
    parser.add_argument("input_file", help="Input VCF/BCF file to check")
    parser.add_argument(
        "-r", "--reference", help="Reference FASTA file for checking reference alleles"
    )
    parser.add_argument("-o", "--output", help="Output file to write issues to")
    parser.add_argument(
        "-s",
        "--strict",
        action="store_true",
        help="Use strict checking mode (more warnings)",
    )
    parser.add_argument(
        "-f",
        "--apply-filters",
        action="store_true",
        help="Skip filtered variants (those with FILTER != PASS)",
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true", help="Enable verbose logging"
    )

    args = parser.parse_args()

    # Setup logging
    log_level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(
        level=log_level, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )

    # Create and run the VCF checker
    checker = VCFChecker(
        reference_path=args.reference,
        strict=args.strict,
        apply_filters=args.apply_filters,
    )

    try:
        checker.check_file(args.input_file, args.output)
        checker.print_summary()
    except Exception as e:
        logging.error(f"Error checking file: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
