#!/usr/bin/env python3
"""
Python implementation of the preprocess module.

This module provides functionality to preprocess VCF files,
including variant decomposition, normalization, and left-shifting.
This is a pure Python implementation that replaces the C++ preprocess component.
"""

import logging
import os
import tempfile
from enum import Enum
from typing import List, Optional, Tuple, Union

import pysam
from pysam import FastaFile, TabixFile, VariantFile, VariantRecord

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class DecomposeLevel(Enum):
    """Enumeration for decompose level options."""

    NONE = 0  # No decomposition
    CONSERVATIVE = 1  # Only decompose when there's no ambiguity
    AGGRESSIVE = 2  # Always decompose to smallest possible units


class PreprocessEngine:
    """
    Engine for preprocessing variant calls in VCF files.

    This class replaces the C++ preprocess component with a pure Python
    implementation using pysam.

    Features:
    - Variant decomposition (split multi-allelic variants)
    - Left-alignment of variants
    - Normalization (trim common prefixes/suffixes)
    - Special handling for haploid regions (e.g., chrX in males)
    """

    def __init__(
        self,
        input_vcf: str,
        reference_fasta: str,
        output_vcf: Optional[str] = None,
        decompose_level: Union[DecomposeLevel, int] = DecomposeLevel.CONSERVATIVE,
        left_shift: bool = True,
        regions: Optional[str] = None,
        haploid_x: bool = False,
        output_bcf: bool = False,
        pass_only: bool = False,
    ):
        """
        Initialize the preprocess engine.

        Args:
            input_vcf: Path to input VCF file
            reference_fasta: Path to reference FASTA file
            output_vcf: Path to output VCF file (default: auto-generate)
            decompose_level: Level of variant decomposition
            left_shift: Whether to left-align variants
            regions: BED file or region string for processing specific regions
            haploid_x: Whether to treat X chromosome as haploid
            output_bcf: Whether to output BCF instead of VCF
            pass_only: Whether to only include variants with PASS filter
        """
        self.input_vcf = input_vcf
        self.reference_fasta = reference_fasta
        self.output_vcf = output_vcf

        # Convert decompose_level to enum if it's an int
        if isinstance(decompose_level, int):
            self.decompose_level = DecomposeLevel(decompose_level)
        else:
            self.decompose_level = decompose_level

        self.left_shift = left_shift
        self.regions = regions
        self.haploid_x = haploid_x
        self.output_bcf = output_bcf
        self.pass_only = pass_only

        # Tracking stats
        self.stats = {
            "total_variants": 0,
            "decomposed_variants": 0,
            "left_shifted_variants": 0,
            "normalized_variants": 0,
            "filtered_variants": 0,
            "haploid_adjusted_variants": 0,
        }

        # Open files
        self._open_files()

    def _open_files(self):
        """Open input/output files and reference."""
        try:
            self.vcf_in = VariantFile(self.input_vcf)
            logger.info(f"Opened input VCF: {self.input_vcf}")
        except Exception as e:
            raise ValueError(f"Failed to open input VCF: {e}")

        try:
            self.reference = FastaFile(self.reference_fasta)
            logger.info(f"Opened reference: {self.reference_fasta}")
        except Exception as e:
            raise ValueError(f"Failed to open reference FASTA: {e}")

        # Open regions file if provided
        self.region_list = None
        if self.regions:
            if os.path.exists(self.regions):
                try:
                    self.region_tabix = TabixFile(self.regions)
                    logger.info(f"Opened regions BED: {self.regions}")
                except Exception as e:
                    logger.warning(f"Failed to open regions as tabix file: {e}")
                    # Try to parse as a BED file
                    self.region_list = self._parse_bed_file(self.regions)
            else:
                # Assume it's a region string like "chr1:1000-2000"
                self.region_list = [self.regions]

        # Create output file if not provided
        if not self.output_vcf:
            suffix = "bcf" if self.output_bcf else "vcf.gz"
            fd, self.output_vcf = tempfile.mkstemp(suffix=f".{suffix}")
            os.close(fd)
            logger.info(f"Created temporary output file: {self.output_vcf}")

        # Prepare output header
        self._prepare_output_header()

    def _prepare_output_header(self):
        """Prepare the output VCF/BCF header."""
        header = self.vcf_in.header.copy()

        # Add INFO fields for preprocessing steps
        header.add_line(
            "##INFO=<ID=DECOMPOSED,Number=0,Type=Flag,"
            'Description="Variant was decomposed from a multi-allelic variant">'
        )
        header.add_line(
            "##INFO=<ID=NORMALIZED,Number=0,Type=Flag,"
            'Description="Variant was normalized by trimming common prefix/suffix">'
        )
        header.add_line(
            "##INFO=<ID=LEFTSHIFTED,Number=0,Type=Flag,"
            'Description="Variant was left-shifted">'
        )
        header.add_line(
            "##INFO=<ID=ORIGINAL_POS,Number=1,Type=Integer,"
            'Description="Original position before preprocessing">'
        )
        header.add_line(
            "##INFO=<ID=ORIGINAL_ALLELES,Number=.,Type=String,"
            'Description="Original alleles before preprocessing">'
        )

        # Add metadata about preprocessing
        header.add_line(f'##preprocess_version="{__version__}"')
        header.add_line(
            f'##preprocess_command="python_preprocess.py '
            f"--decompose {self.decompose_level.value} "
            f"--left-shift {str(self.left_shift).lower()} "
            f'--haploid-x {str(self.haploid_x).lower()}"'
        )

        # Create output file
        mode = "wb" if self.output_bcf else "wz"
        self.vcf_out = VariantFile(self.output_vcf, mode=mode, header=header)

    def _parse_bed_file(self, bed_file: str) -> List[str]:
        """Parse a BED file into a list of region strings."""
        regions = []
        with open(bed_file) as f:
            for line in f:
                if line.startswith("#") or not line.strip():
                    continue

                fields = line.strip().split("\t")
                if len(fields) >= 3:
                    chrom, start, end = fields[0], int(fields[1]), int(fields[2])
                    regions.append(f"{chrom}:{start}-{end}")

        return regions

    def is_in_regions(self, record: VariantRecord) -> bool:
        """Check if a variant record is in the specified regions."""
        if not self.regions:
            return True

        if self.region_tabix:
            # Use tabix for efficient region checking
            try:
                for _ in self.region_tabix.fetch(
                    record.contig, record.pos - 1, record.pos
                ):
                    return True
                return False
            except ValueError:
                # Contig not in the index
                return False
        elif self.region_list:
            # Check against parsed regions
            for region in self.region_list:
                if ":" in region and "-" in region:
                    # Parse region string
                    r_chrom, pos_range = region.split(":", 1)
                    r_start, r_end = map(int, pos_range.split("-", 1))

                    if (
                        record.contig == r_chrom
                        and record.pos >= r_start
                        and record.pos <= r_end
                    ):
                        return True
                elif record.contig == region:
                    # Whole chromosome
                    return True

            return False

        return True

    def normalize_variant(
        self, chrom: str, pos: int, ref: str, alt: str
    ) -> Tuple[int, str, str]:
        """
        Normalize a variant by trimming common prefixes/suffixes.

        Args:
            chrom: Chromosome name
            pos: 1-based position
            ref: Reference allele
            alt: Alternative allele

        Returns:
            Tuple of (normalized_pos, normalized_ref, normalized_alt)
        """
        # Skip if one of the alleles is symbolic
        if alt.startswith("<") or ref.startswith("<"):
            return pos, ref, alt

        # Skip if one of the alleles is empty
        if not ref or not alt:
            return pos, ref, alt

        # Trim common suffix
        while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
            ref = ref[:-1]
            alt = alt[:-1]

        # Trim common prefix and adjust position
        new_pos = pos
        while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
            ref = ref[1:]
            alt = alt[1:]
            new_pos += 1

        # Make sure we don't end up with empty alleles
        if not ref or not alt:
            ref = "."
            alt = "."

        return new_pos, ref, alt

    def left_shift_variant(
        self, chrom: str, pos: int, ref: str, alt: str
    ) -> Tuple[int, str, str]:
        """
        Left-shift a variant as far as possible.

        Args:
            chrom: Chromosome name
            pos: 1-based position
            ref: Reference allele
            alt: Alternative allele

        Returns:
            Tuple of (shifted_pos, shifted_ref, shifted_alt)
        """
        # Skip if one of the alleles is symbolic
        if alt.startswith("<") or ref.startswith("<"):
            return pos, ref, alt

        # Skip if one of the alleles is empty or ref=alt
        if not ref or not alt or ref == alt:
            return pos, ref, alt

        # Left-shift logic: only applies to indels
        is_indel = len(ref) != len(alt)
        if not is_indel or not self.left_shift:
            return pos, ref, alt

        # Get the sequence context
        try:
            # We need to look backwards to allow for shifting
            # Get more sequence than we'd theoretically need
            context_size = max(100, 2 * max(len(ref), len(alt)))
            seq_start = max(0, pos - context_size - 1)
            seq_end = pos + len(ref) + 10

            # Get the sequence (0-based)
            sequence = self.reference.fetch(chrom, seq_start, seq_end)

            # Convert to 0-based for calculations
            pos_0 = pos - 1 - seq_start
            ref_in_seq = sequence[pos_0 : pos_0 + len(ref)]

            # Sanity check that our sequence fetch worked
            if ref_in_seq.upper() != ref.upper():
                logger.warning(
                    f"Reference mismatch at {chrom}:{pos}: expected '{ref}', "
                    f"got '{ref_in_seq}' - skipping left-shift"
                )
                return pos, ref, alt

            # Left-shift logic
            shifted = False
            new_pos = pos
            new_ref = ref
            new_alt = alt

            # If it's a deletion
            if len(ref) > len(alt):
                # Get the deleted sequence
                deleted = ref[len(alt) :]
                prefix = alt

                # Try to shift left
                current_pos = pos_0
                while current_pos > 0:
                    # Check if the base to the left matches the last base of the deletion
                    if sequence[current_pos - 1] == deleted[-1]:
                        # Shift left
                        deleted = deleted[:-1] + sequence[current_pos - 1]
                        current_pos -= 1
                        shifted = True
                    else:
                        break

                if shifted:
                    # Recalculate ref/alt based on the new position
                    new_pos = seq_start + current_pos + 1  # Convert back to 1-based
                    new_ref = prefix + deleted
                    new_alt = prefix

            # If it's an insertion
            elif len(ref) < len(alt):
                # Get the inserted sequence
                inserted = alt[len(ref) :]
                prefix = ref

                # Try to shift left
                current_pos = pos_0
                while current_pos > 0:
                    # Check if the base to the left matches the last base of the insertion
                    if sequence[current_pos - 1] == inserted[-1]:
                        # Shift left
                        inserted = inserted[:-1] + sequence[current_pos - 1]
                        current_pos -= 1
                        shifted = True
                    else:
                        break

                if shifted:
                    # Recalculate ref/alt based on the new position
                    new_pos = seq_start + current_pos + 1  # Convert back to 1-based
                    new_ref = prefix
                    new_alt = prefix + inserted

            # Normalize after left-shifting
            if shifted:
                new_pos, new_ref, new_alt = self.normalize_variant(
                    chrom, new_pos, new_ref, new_alt
                )
                self.stats["left_shifted_variants"] += 1

            return new_pos, new_ref, new_alt

        except Exception as e:
            logger.warning(f"Error during left-shifting at {chrom}:{pos}: {e}")
            return pos, ref, alt

    def decompose_variant(self, record: VariantRecord) -> List[VariantRecord]:
        """
        Decompose a multi-allelic variant into multiple bi-allelic variants.

        Args:
            record: Input variant record

        Returns:
            List of decomposed variant records
        """
        # Skip if decomposition is disabled or there's only one alt allele
        if self.decompose_level == DecomposeLevel.NONE or len(record.alts) <= 1:
            return [record]

        result = []
        original_pos = record.pos
        original_alleles = [record.ref] + list(record.alts)

        for i, alt in enumerate(record.alts):
            # Create a new record for each alt allele
            new_record = self.vcf_out.new_record()

            # Copy basic fields
            new_record.contig = record.contig
            new_record.pos = record.pos
            new_record.id = record.id
            new_record.ref = record.ref
            new_record.alts = (alt,)

            # Handle filters
            if self.pass_only and any(f != "PASS" for f in record.filter.keys()):
                continue

            for f in record.filter.keys():
                new_record.filter.add(f)

            # Copy all INFO fields
            for key, value in record.info.items():
                if key in new_record.header.info:
                    new_record.info[key] = value

            # Add decomposition INFO
            new_record.info["DECOMPOSED"] = True
            new_record.info["ORIGINAL_POS"] = original_pos
            new_record.info["ORIGINAL_ALLELES"] = original_alleles

            # Copy format fields and sample data
            for sample in record.samples:
                # Make sure the sample exists in the new record
                if sample in new_record.samples:
                    # Copy all format fields
                    for field in record.format.keys():
                        if field in new_record.header.formats:
                            # Handle GT field specially to select the right allele
                            if field == "GT":
                                gt = record.samples[sample]["GT"]
                                # Adjust GT to be bi-allelic (0 = ref, 1 = this alt)
                                new_gt = []
                                for g in gt:
                                    if g is None:
                                        new_gt.append(None)
                                    elif g == 0:
                                        new_gt.append(0)  # Reference allele stays 0
                                    elif g == i + 1:
                                        new_gt.append(1)  # This alt becomes 1
                                    else:
                                        new_gt.append(0)  # Other alts become 0

                                new_record.samples[sample]["GT"] = tuple(new_gt)
                            else:
                                # For other fields, just copy the value
                                try:
                                    new_record.samples[sample][field] = record.samples[
                                        sample
                                    ][field]
                                except Exception as e:
                                    logger.warning(
                                        f"Failed to copy format field {field} for sample {sample}: {e}"
                                    )

            # Normalize and left-shift if needed
            if self.left_shift or self.decompose_level == DecomposeLevel.AGGRESSIVE:
                new_pos, new_ref, new_alt = self.left_shift_variant(
                    new_record.contig, new_record.pos, new_record.ref, alt
                )

                if (
                    new_pos != new_record.pos
                    or new_ref != new_record.ref
                    or new_alt != alt
                ):
                    # Position or alleles changed, update the record
                    new_record.pos = new_pos
                    new_record.ref = new_ref
                    new_record.alts = (new_alt,)
                    new_record.info["LEFTSHIFTED"] = True

                # Normalize the variant
                new_pos, new_ref, new_alt = self.normalize_variant(
                    new_record.contig, new_record.pos, new_record.ref, new_alt
                )

                if (
                    new_pos != new_record.pos
                    or new_ref != new_record.ref
                    or new_alt != new_record.alts[0]
                ):
                    # Position or alleles changed, update the record
                    new_record.pos = new_pos
                    new_record.ref = new_ref
                    new_record.alts = (new_alt,)
                    new_record.info["NORMALIZED"] = True
                    self.stats["normalized_variants"] += 1

            result.append(new_record)

        if len(result) > 1:
            self.stats["decomposed_variants"] += 1

        return result

    def handle_haploid_regions(self, record: VariantRecord) -> VariantRecord:
        """
        Handle special cases for haploid regions like chrX in males.

        Args:
            record: Input variant record

        Returns:
            Modified variant record
        """
        if not self.haploid_x:
            return record

        # Check if this is chromosome X
        chrom = record.contig.lower().replace("chr", "")
        if chrom != "x":
            return record

        # Only process if GT field exists
        if "GT" not in record.format:
            return record

        # Convert heterozygous genotypes to homozygous
        for sample in record.samples:
            gt = record.samples[sample]["GT"]
            if gt is None or len(gt) < 2:
                continue

            # Check if it's heterozygous
            if gt[0] != gt[1] and gt[0] is not None and gt[1] is not None:
                # Convert to homozygous for the alt allele if present
                if gt[0] > 0:
                    new_gt = (gt[0], gt[0])
                elif gt[1] > 0:
                    new_gt = (gt[1], gt[1])
                else:
                    new_gt = (0, 0)

                record.samples[sample]["GT"] = new_gt
                self.stats["haploid_adjusted_variants"] += 1

        return record

    def process(self) -> str:
        """
        Process the input VCF file.

        Returns:
            Path to the output VCF/BCF file
        """
        logger.info(
            f"Starting preprocessing with options: "
            f"decompose={self.decompose_level.name}, "
            f"left_shift={self.left_shift}, "
            f"haploid_x={self.haploid_x}"
        )

        # Process each variant
        for record in self.vcf_in:
            self.stats["total_variants"] += 1

            # Skip if variant doesn't pass filters
            if self.pass_only and any(f != "PASS" for f in record.filter.keys()):
                self.stats["filtered_variants"] += 1
                continue

            # Skip if not in specified regions
            if not self.is_in_regions(record):
                continue

            # Decompose multi-allelic variants
            decomposed_records = self.decompose_variant(record)

            # Process each decomposed record
            for decomp_record in decomposed_records:
                # Handle haploid regions
                if self.haploid_x:
                    decomp_record = self.handle_haploid_regions(decomp_record)

                # Write to output
                self.vcf_out.write(decomp_record)

        # Close files
        self.vcf_out.close()
        self.vcf_in.close()

        # Index the output file if it's VCF.GZ
        if not self.output_bcf and self.output_vcf.endswith(".vcf.gz"):
            try:
                pysam.tabix_index(self.output_vcf, preset="vcf")
                logger.info(f"Indexed output VCF: {self.output_vcf}")
            except Exception as e:
                logger.warning(f"Failed to index output VCF: {e}")

        # Log statistics
        logger.info("Preprocessing statistics:")
        logger.info(f"  Total variants processed: {self.stats['total_variants']}")
        logger.info(f"  Decomposed variants: {self.stats['decomposed_variants']}")
        logger.info(f"  Left-shifted variants: {self.stats['left_shifted_variants']}")
        logger.info(f"  Normalized variants: {self.stats['normalized_variants']}")
        logger.info(
            f"  Haploid adjusted variants: {self.stats['haploid_adjusted_variants']}"
        )
        logger.info(f"  Filtered variants: {self.stats['filtered_variants']}")

        return self.output_vcf


# Define version
__version__ = "1.0.0"


def main():
    """Command-line interface for preprocess."""
    import argparse

    parser = argparse.ArgumentParser(
        description="Preprocess VCF files - normalize, decompose, and left-shift variants"
    )
    parser.add_argument("input_vcf", help="Input VCF/BCF file")
    parser.add_argument(
        "-r", "--reference", required=True, help="Reference FASTA file (required)"
    )
    parser.add_argument(
        "-o", "--output", help="Output VCF/BCF file (default: auto-generate)"
    )
    parser.add_argument(
        "-V",
        "--decompose",
        type=int,
        default=1,
        choices=[0, 1, 2],
        help="Decomposition level: 0=none, 1=conservative, 2=aggressive (default: 1)",
    )
    parser.add_argument(
        "-L",
        "--left-shift",
        action="store_true",
        default=True,
        help="Left-align variants (default: True)",
    )
    parser.add_argument(
        "--no-left-shift",
        action="store_false",
        dest="left_shift",
        help="Disable left-alignment of variants",
    )
    parser.add_argument(
        "-l",
        "--regions",
        help="Regions to process (BED file or region string like 'chr1:1000-2000')",
    )
    parser.add_argument(
        "--haploid-x",
        action="store_true",
        help="Treat X chromosome as haploid (convert het to hom)",
    )
    parser.add_argument("--bcf", action="store_true", help="Output BCF instead of VCF")
    parser.add_argument(
        "--pass-only",
        action="store_true",
        help="Only include variants with PASS filter",
    )
    parser.add_argument("--verbose", action="store_true", help="Enable verbose logging")

    args = parser.parse_args()

    # Set logging level
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    # Create and run the preprocess engine
    try:
        engine = PreprocessEngine(
            input_vcf=args.input_vcf,
            reference_fasta=args.reference,
            output_vcf=args.output,
            decompose_level=args.decompose,
            left_shift=args.left_shift,
            regions=args.regions,
            haploid_x=args.haploid_x,
            output_bcf=args.bcf,
            pass_only=args.pass_only,
        )

        output_file = engine.process()
        print(f"Preprocessing complete: {output_file}")

    except Exception as e:
        logger.error(f"Error during preprocessing: {e}")
        return 1

    return 0


if __name__ == "__main__":
    import sys

    sys.exit(main())
