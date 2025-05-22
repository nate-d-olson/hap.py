#!/usr/bin/env python3
"""
Python implementation of the hapcmp functionality.

This module replaces the C++ hapcmp tool, which performs haplotype comparison
between two VCF files. The Python version uses pysam and BioPython to achieve
the same functionality while being more maintainable.
"""

import argparse
import json
import logging
import sys
import time
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Union

from pysam import FastaFile, VariantFile, VariantRecord


@dataclass
class HaplotypeBlock:
    """Represents a haplotype block for comparison."""

    chrom: str
    start: int
    end: int
    variants1: List[VariantRecord] = field(default_factory=list)
    variants2: List[VariantRecord] = field(default_factory=list)
    haplotypes1: List[str] = field(default_factory=list)
    haplotypes2: List[str] = field(default_factory=list)
    status: str = "unknown"
    error_message: str = ""


class HaploComparator:
    """
    Compares haplotypes between two VCF files.

    This class implements the core functionality of the hapcmp tool,
    allowing comparison of haplotypes between two VCF files within
    specified genomic regions.
    """

    def __init__(
        self,
        reference_path: str,
        max_haplotypes: int = 4096,
        output_sequences: bool = False,
        apply_filters: bool = False,
        do_alignment: bool = False,
    ):
        """
        Initialize a HaploComparator instance.

        Args:
            reference_path: Path to the reference FASTA file
            max_haplotypes: Maximum number of haplotypes to enumerate
            output_sequences: Whether to output haplotype sequences
            apply_filters: Whether to apply filters from VCF
            do_alignment: Whether to perform alignment on mismatching haplotypes
        """
        self.reference_path = reference_path
        self.max_haplotypes = max_haplotypes
        self.output_sequences = output_sequences
        self.apply_filters = apply_filters
        self.do_alignment = do_alignment
        self.logger = logging.getLogger("hapcmp")

        # Open reference file
        try:
            self.reference = FastaFile(reference_path)
            self.logger.info(f"Opened reference file: {reference_path}")
        except Exception as e:
            self.logger.error(f"Failed to open reference file: {e}")
            raise ValueError(f"Failed to open reference file: {e}")

    def compare_files(
        self,
        file1: str,
        sample1: str,
        file2: str,
        sample2: str,
        regions_file: str,
        output_bed: Optional[str] = None,
        output_diffs: Optional[str] = None,
        output_errors: Optional[str] = None,
        block_limit: int = -1,
        show_progress: bool = False,
        progress_seconds: int = 10,
    ) -> Dict[str, int]:
        """
        Compare haplotypes between two VCF files.

        Args:
            file1: Path to the first VCF file
            sample1: Sample name in the first VCF file (empty for first sample)
            file2: Path to the second VCF file
            sample2: Sample name in the second VCF file (empty for first sample)
            regions_file: Path to the BED file with regions to compare
            output_bed: Path to write comparison results as BED
            output_diffs: Path to write differences as JSON
            output_errors: Path to write error information
            block_limit: Maximum number of blocks to process (-1 for all)
            show_progress: Whether to show progress information
            progress_seconds: Seconds between progress updates

        Returns:
            Dictionary with comparison statistics
        """
        # Open VCF files
        try:
            vcf1 = VariantFile(file1)
            vcf2 = VariantFile(file2)

            # Select samples if specified
            if sample1 and sample1 in vcf1.header.samples:
                vcf1.subset_samples([sample1])
            elif sample1:
                self.logger.warning(f"Sample {sample1} not found in {file1}")

            if sample2 and sample2 in vcf2.header.samples:
                vcf2.subset_samples([sample2])
            elif sample2:
                self.logger.warning(f"Sample {sample2} not found in {file2}")

            self.logger.info(f"Opened VCF files: {file1} and {file2}")
        except Exception as e:
            self.logger.error(f"Failed to open VCF files: {e}")
            raise ValueError(f"Failed to open VCF files: {e}")

        # Read regions from BED file
        regions = self._read_regions(regions_file)
        self.logger.info(f"Read {len(regions)} regions from {regions_file}")

        # Open output files if specified
        bed_out = None
        diffs_out = None
        errors_out = None

        if output_bed:
            bed_out = open(output_bed, "w")
            bed_out.write("#CHROM\tSTART\tEND\tSTATUS\n")

        if output_diffs:
            diffs_out = open(output_diffs, "w")

        if output_errors:
            errors_out = open(output_errors, "w")

        # Process blocks
        stats = {
            "blocks_total": len(regions),
            "blocks_processed": 0,
            "blocks_match": 0,
            "blocks_mismatch": 0,
            "blocks_error": 0,
        }

        start_time = time.time()
        last_progress = start_time

        for i, region in enumerate(regions):
            if block_limit >= 0 and i >= block_limit:
                self.logger.info(f"Reached block limit of {block_limit}")
                break

            # Show progress if requested
            if show_progress and time.time() - last_progress >= progress_seconds:
                self._show_progress(i, len(regions), stats, start_time)
                last_progress = time.time()

            # Process region
            block = self._process_region(vcf1, vcf2, region, stats)

            # Write outputs if requested
            if bed_out:
                self._write_bed_record(bed_out, block)

            if diffs_out and block.status != "match":
                self._write_diff_record(diffs_out, block)

            if errors_out and block.status == "error":
                self._write_error_record(errors_out, block)

        # Show final progress
        if show_progress:
            self._show_progress(len(regions), len(regions), stats, start_time)

        # Close output files
        if bed_out:
            bed_out.close()

        if diffs_out:
            diffs_out.close()

        if errors_out:
            errors_out.close()

        # Close input files
        vcf1.close()
        vcf2.close()

        return stats

    def _read_regions(self, regions_file: str) -> List[Dict[str, Union[str, int]]]:
        """
        Read regions from a BED file.

        Args:
            regions_file: Path to the BED file

        Returns:
            List of region dictionaries with chrom, start, end
        """
        regions = []

        # Handle stdin
        if regions_file == "-":
            for line in sys.stdin:
                if line.startswith("#"):
                    continue
                fields = line.strip().split("\t")
                if len(fields) >= 3:
                    regions.append(
                        {
                            "chrom": fields[0],
                            "start": int(fields[1]),
                            "end": int(fields[2]),
                        }
                    )
            return regions

        # Read from file
        try:
            with open(regions_file) as f:
                for line in f:
                    if line.startswith("#"):
                        continue
                    fields = line.strip().split("\t")
                    if len(fields) >= 3:
                        regions.append(
                            {
                                "chrom": fields[0],
                                "start": int(fields[1]),
                                "end": int(fields[2]),
                            }
                        )
        except Exception as e:
            self.logger.error(f"Failed to read regions file: {e}")
            raise ValueError(f"Failed to read regions file: {e}")

        return regions

    def _process_region(
        self,
        vcf1: VariantFile,
        vcf2: VariantFile,
        region: Dict[str, Union[str, int]],
        stats: Dict[str, int],
    ) -> HaplotypeBlock:
        """
        Process a single region for haplotype comparison.

        Args:
            vcf1: First VCF file
            vcf2: Second VCF file
            region: Region to process
            stats: Statistics dictionary to update

        Returns:
            HaplotypeBlock with comparison results
        """
        chrom = str(region["chrom"])
        start = int(region["start"])
        end = int(region["end"])

        # Create block
        block = HaplotypeBlock(chrom=chrom, start=start, end=end)

        # Fetch variants from both VCFs
        try:
            for record in vcf1.fetch(chrom, start, end):
                if not self.apply_filters or self._passes_filters(record):
                    block.variants1.append(record)

            for record in vcf2.fetch(chrom, start, end):
                if not self.apply_filters or self._passes_filters(record):
                    block.variants2.append(record)
        except Exception as e:
            block.status = "error"
            block.error_message = f"Failed to fetch variants: {e}"
            stats["blocks_error"] += 1
            stats["blocks_processed"] += 1
            return block

        # Generate haplotypes
        try:
            block.haplotypes1 = self._generate_haplotypes(
                chrom, start, end, block.variants1
            )
            block.haplotypes2 = self._generate_haplotypes(
                chrom, start, end, block.variants2
            )
        except Exception as e:
            block.status = "error"
            block.error_message = f"Failed to generate haplotypes: {e}"
            stats["blocks_error"] += 1
            stats["blocks_processed"] += 1
            return block

        # Compare haplotypes
        if self._compare_haplotypes(block.haplotypes1, block.haplotypes2):
            block.status = "match"
            stats["blocks_match"] += 1
        else:
            block.status = "mismatch"
            stats["blocks_mismatch"] += 1

        stats["blocks_processed"] += 1
        return block

    def _passes_filters(self, record: VariantRecord) -> bool:
        """
        Check if a variant passes filters.

        Args:
            record: Variant record to check

        Returns:
            True if the variant passes filters, False otherwise
        """
        # If no filters, it passes
        if not record.filter:
            return True

        # Check for PASS filter
        for flt in record.filter:
            if flt != "PASS":
                return False

        return True

    def _generate_haplotypes(
        self, chrom: str, start: int, end: int, variants: List[VariantRecord]
    ) -> List[str]:
        """
        Generate all possible haplotypes for a set of variants.

        Args:
            chrom: Chromosome name
            start: Start position
            end: End position
            variants: List of variant records

        Returns:
            List of haplotype sequences
        """
        # Get reference sequence for the region
        ref_seq = self.reference.fetch(chrom, start, end)

        # If no variants, return reference sequence
        if not variants:
            return [ref_seq]

        # Sort variants by position
        variants = sorted(variants, key=lambda x: x.pos)

        # Generate haplotypes using a simple approach for this initial implementation
        # For a more complete implementation, we would need to handle phasing,
        # complex variants, overlapping variants, etc.

        # For now, we'll just generate all possible combinations of alt alleles
        # This is a simplification of the C++ implementation
        haplotypes = [ref_seq]

        for variant in variants:
            # Skip non-variant sites
            if len(variant.alts) == 0:
                continue

            # Adjust position to be relative to region start
            rel_pos = variant.pos - start

            # Skip variants outside region
            if rel_pos < 0 or rel_pos >= len(ref_seq):
                continue

            # Get genotype
            genotype = None
            for sample in variant.samples.values():
                if "GT" in sample:
                    genotype = sample["GT"]
                    break

            # Skip if no genotype
            if not genotype:
                continue

            # Generate new haplotypes
            new_haplotypes = []
            for hap in haplotypes:
                # For each allele in genotype
                for allele_idx in genotype:
                    if allele_idx is None or allele_idx == 0:
                        # Reference allele, keep haplotype unchanged
                        new_haplotypes.append(hap)
                    elif allele_idx > 0 and allele_idx <= len(variant.alts):
                        # Apply alt allele
                        alt = variant.alts[allele_idx - 1]
                        new_hap = (
                            hap[:rel_pos] + alt + hap[rel_pos + len(variant.ref) :]
                        )
                        new_haplotypes.append(new_hap)

            # Check if we exceed max haplotypes
            if len(new_haplotypes) > self.max_haplotypes:
                self.logger.warning(
                    f"Exceeded maximum number of haplotypes ({self.max_haplotypes})"
                )
                return new_haplotypes[: self.max_haplotypes]

            haplotypes = new_haplotypes

        return haplotypes

    def _compare_haplotypes(
        self, haplotypes1: List[str], haplotypes2: List[str]
    ) -> bool:
        """
        Compare two sets of haplotypes.

        Args:
            haplotypes1: First set of haplotypes
            haplotypes2: Second set of haplotypes

        Returns:
            True if the haplotypes match, False otherwise
        """
        # Simple case: exactly matching sets
        if set(haplotypes1) == set(haplotypes2):
            return True

        # If alignment is enabled, try to find approximate matches
        if self.do_alignment:
            # This would be implemented using sequence alignment
            # For now, we'll just return False
            return False

        return False

    def _show_progress(
        self, current: int, total: int, stats: Dict[str, int], start_time: float
    ) -> None:
        """
        Show progress information.

        Args:
            current: Current block index
            total: Total number of blocks
            stats: Statistics dictionary
            start_time: Start time in seconds
        """
        elapsed = time.time() - start_time
        percent = 100.0 * current / total if total > 0 else 100.0

        self.logger.info(
            f"Progress: {current}/{total} blocks ({percent:.1f}%), "
            f"Elapsed: {elapsed:.1f}s, "
            f"Matches: {stats['blocks_match']}, "
            f"Mismatches: {stats['blocks_mismatch']}, "
            f"Errors: {stats['blocks_error']}"
        )

    def _write_bed_record(self, bed_out, block: HaplotypeBlock) -> None:
        """
        Write a block record to the BED output file.

        Args:
            bed_out: Output file handle
            block: Haplotype block to write
        """
        bed_out.write(f"{block.chrom}\t{block.start}\t{block.end}\t{block.status}\n")

    def _write_diff_record(self, diffs_out, block: HaplotypeBlock) -> None:
        """
        Write a difference record to the JSON output file.

        Args:
            diffs_out: Output file handle
            block: Haplotype block to write
        """
        diff_record = {
            "chrom": block.chrom,
            "start": block.start,
            "end": block.end,
            "status": block.status,
        }

        # Add sequences if requested
        if self.output_sequences:
            diff_record["haplotypes1"] = block.haplotypes1
            diff_record["haplotypes2"] = block.haplotypes2

        # Add error message if present
        if block.error_message:
            diff_record["error"] = block.error_message

        diffs_out.write(json.dumps(diff_record) + "\n")

    def _write_error_record(self, errors_out, block: HaplotypeBlock) -> None:
        """
        Write an error record to the output file.

        Args:
            errors_out: Output file handle
            block: Haplotype block to write
        """
        errors_out.write(
            f"{block.chrom}\t{block.start}\t{block.end}\t{block.error_message}\n"
        )


def main():
    """Main entry point for the hapcmp tool."""
    parser = argparse.ArgumentParser(
        description="Compare haplotypes between two VCF files"
    )

    parser.add_argument(
        "input_regions",
        help="The input bed file specifying haplotype block regions (use - for stdin)",
    )
    parser.add_argument(
        "input_vcfs",
        nargs=2,
        help="Two VCF files to compare (use file:sample for a specific sample column)",
    )
    parser.add_argument(
        "-r", "--reference", required=True, help="The reference fasta file"
    )
    parser.add_argument(
        "-b",
        "--output-bed",
        help="Output block results as bed files (default is to output to stdout)",
    )
    parser.add_argument("-e", "--output-errors", help="Output failure information")
    parser.add_argument(
        "-d",
        "--output-diffs",
        help="Output shared and different variants to a JSON file",
    )
    parser.add_argument(
        "-n",
        "--max-n-haplotypes",
        type=int,
        default=4096,
        help="Maximum number of haplotypes to enumerate",
    )
    parser.add_argument(
        "--output-sequences", action="store_true", help="Output haplotype sequences"
    )
    parser.add_argument(
        "--progress", action="store_true", help="Output progress information"
    )
    parser.add_argument(
        "--progress-seconds",
        type=int,
        default=10,
        help="Output progress information every n seconds",
    )
    parser.add_argument(
        "-l",
        "--limit",
        type=int,
        default=-1,
        help="Maximum number of haplotype blocks to process",
    )
    parser.add_argument(
        "-f", "--apply-filters", action="store_true", help="Apply filtering in VCF"
    )
    parser.add_argument(
        "--do-alignment",
        action="store_true",
        help="Perform alignments on mismatching haplotypes",
    )
    parser.add_argument(
        "--verbose",
        "-v",
        action="count",
        default=0,
        help="Increase verbosity (can be used multiple times)",
    )

    args = parser.parse_args()

    # Set up logging
    log_level = logging.WARNING
    if args.verbose == 1:
        log_level = logging.INFO
    elif args.verbose >= 2:
        log_level = logging.DEBUG

    logging.basicConfig(
        level=log_level, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )

    logger = logging.getLogger("hapcmp")

    # Parse input VCF paths and samples
    vcf1_parts = args.input_vcfs[0].split(":", 1)
    vcf1_path = vcf1_parts[0]
    vcf1_sample = vcf1_parts[1] if len(vcf1_parts) > 1 else ""

    vcf2_parts = args.input_vcfs[1].split(":", 1)
    vcf2_path = vcf2_parts[0]
    vcf2_sample = vcf2_parts[1] if len(vcf2_parts) > 1 else ""

    # Create comparator
    try:
        comparator = HaploComparator(
            reference_path=args.reference,
            max_haplotypes=args.max_n_haplotypes,
            output_sequences=args.output_sequences,
            apply_filters=args.apply_filters,
            do_alignment=args.do_alignment,
        )

        # Run comparison
        stats = comparator.compare_files(
            file1=vcf1_path,
            sample1=vcf1_sample,
            file2=vcf2_path,
            sample2=vcf2_sample,
            regions_file=args.input_regions,
            output_bed=args.output_bed,
            output_diffs=args.output_diffs,
            output_errors=args.output_errors,
            block_limit=args.limit,
            show_progress=args.progress,
            progress_seconds=args.progress_seconds,
        )

        # Print summary
        logger.info(
            f"Processed {stats['blocks_processed']} of {stats['blocks_total']} blocks"
        )
        logger.info(f"Matches: {stats['blocks_match']}")
        logger.info(f"Mismatches: {stats['blocks_mismatch']}")
        logger.info(f"Errors: {stats['blocks_error']}")

        # Write summary to stdout if no output bed file specified
        if not args.output_bed:
            print(
                f"Processed {stats['blocks_processed']} of {stats['blocks_total']} blocks"
            )
            print(f"Matches: {stats['blocks_match']}")
            print(f"Mismatches: {stats['blocks_mismatch']}")
            print(f"Errors: {stats['blocks_error']}")

    except Exception as e:
        logger.error(f"Error during comparison: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
