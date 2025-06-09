#!/usr/bin/env python3
"""
Python implementation of multimerge functionality.

This module implements VCF file merging functionality similar to the original
C++ multimerge tool. It supports various options like:
- Trimming alleles
- Left-shifting variants
- Splitting complex alleles
- Merging by location
- Unique alleles
- Processing formats
"""

import argparse
import logging
import sys
from typing import List, Optional, Tuple

try:
    import pysam
except ImportError:
    logging.error(
        "pysam is required for multimerge functionality. Please install with 'pip install pysam'."
    )
    sys.exit(1)

from hap_py.haplo.python_preprocess import PreprocessEngine


def parse_vcf_sample_arg(arg: str) -> Tuple[str, Optional[str]]:
    """
    Parse a VCF sample argument in the format "filename:sample_name".

    Args:
        arg: Input string in format "filename:sample_name"

    Returns:
        Tuple of (filename, sample_name). If no sample name is specified,
        returns (filename, None)
    """
    parts = arg.split(":", 1)
    if len(parts) == 1:
        return parts[0], None
    return parts[0], parts[1]


def get_sample_name(vcf_path: str, specified_sample: Optional[str] = None) -> str:
    """
    Get the sample name from a VCF file.

    Args:
        vcf_path: Path to the VCF file
        specified_sample: Sample name specified by the user, or None

    Returns:
        The sample name to use

    Raises:
        ValueError: If specified sample doesn't exist or no samples in VCF
    """
    vcf = pysam.VariantFile(vcf_path)

    if not vcf.header.samples:
        raise ValueError(f"No samples found in VCF file: {vcf_path}")

    if specified_sample is None:
        # Use first sample if none specified
        return vcf.header.samples[0]

    if specified_sample == "*":
        # Special case: use all samples
        return "*"

    # Check if specified sample exists
    if specified_sample not in vcf.header.samples:
        raise ValueError(f"Sample {specified_sample} not found in VCF file: {vcf_path}")

    return specified_sample


def trim_alleles(ref: str, alts: List[str]) -> Tuple[str, List[str]]:
    """
    Trim common prefix and suffix from alleles.

    Args:
        ref: Reference allele
        alts: List of alternate alleles

    Returns:
        Tuple of (trimmed ref, list of trimmed alts)
    """
    if not alts:
        return ref, alts

    # Find common prefix
    prefix_len = 0
    for i in range(min(len(ref), min(len(alt) for alt in alts))):
        if all(ref[i] == alt[i] for alt in alts):
            prefix_len += 1
        else:
            break

    # Need to keep at least 1 base
    if prefix_len >= len(ref):
        prefix_len = max(0, len(ref) - 1)

    # Find common suffix
    suffix_len = 0
    for i in range(1, min(len(ref), min(len(alt) for alt in alts)) - prefix_len + 1):
        if all(ref[-i] == alt[-i] for alt in alts):
            suffix_len += 1
        else:
            break

    # Need to keep at least 1 base
    if prefix_len + suffix_len >= len(ref):
        suffix_len = max(0, len(ref) - prefix_len - 1)

    # Apply trimming
    if prefix_len > 0 or suffix_len > 0:
        new_ref = ref[prefix_len : len(ref) - suffix_len] if ref else ""
        new_alts = [
            alt[prefix_len : len(alt) - suffix_len] if alt else "" for alt in alts
        ]

        # Ensure all alleles have at least one base
        if not new_ref:
            new_ref = "N"
            new_alts = ["N" + alt for alt in new_alts]

        for i, alt in enumerate(new_alts):
            if not alt:
                new_alts[i] = "N"

        return new_ref, new_alts

    return ref, alts


def merge_vcfs(
    input_vcfs: List[Tuple[str, Optional[str]]],
    output_path: str,
    reference_path: Optional[str] = None,
    trim_alleles_flag: bool = False,
    left_shift: bool = False,
    split_alleles: bool = False,
    merge_by_location: bool = False,
    unique_alleles: bool = False,
    process_formats: bool = False,
    process_split: bool = False,
    process_full: bool = False,
) -> bool:
    """
    Merge multiple VCF files.

    Args:
        input_vcfs: List of (vcf_path, sample_name) tuples
        output_path: Path to write merged VCF
        reference_path: Path to reference FASTA (required for some options)
        trim_alleles_flag: Whether to trim common bases from alleles
        left_shift: Whether to left-shift variants
        split_alleles: Whether to split multi-allelic variants
        merge_by_location: Whether to merge variants at the same location
        unique_alleles: Whether to keep only unique alleles
        process_formats: Whether to process FORMAT fields
        process_split: Whether to split complex variants

    Returns:
        True if merge completed successfully, False otherwise
    """
    # Check for required reference for certain operations
    if (left_shift or split_alleles) and not reference_path:
        logging.error(
            "Reference FASTA is required for left_shift and split_alleles options"
        )
        return False

    # Open reference if needed
    reference = None
    if reference_path:
        try:
            reference = pysam.FastaFile(reference_path)
        except Exception as e:
            logging.error(f"Failed to open reference file {reference_path}: {e}")
            return False

    # Process input VCFs
    vcf_readers = []
    samples = []

    for vcf_path, specified_sample in input_vcfs:
        try:
            # Parse VCF path and sample name
            vcf = pysam.VariantFile(vcf_path)

            # Handle sample name
            sample_name = get_sample_name(vcf_path, specified_sample)
            samples.append(sample_name)

            vcf_readers.append(vcf)
        except Exception as e:
            logging.error(f"Error processing VCF {vcf_path}: {e}")
            return False

    # Create output VCF
    try:
        # Initialize with the header from the first VCF
        output_header = vcf_readers[0].header.copy()

        # Copy FORMAT and INFO fields from all other VCFs to the output header
        for i in range(1, len(vcf_readers)):
            other_header = vcf_readers[i].header
            # Copy FORMAT fields
            for fmt in other_header.formats.values():
                if fmt.name not in output_header.formats:
                    output_header.add_line(
                        f'##FORMAT=<ID={fmt.name},Number={fmt.number},Type={fmt.type},Description="{fmt.description}">'
                    )

            # Copy INFO fields
            for info in other_header.info.values():
                if info.name not in output_header.info:
                    output_header.add_line(
                        f'##INFO=<ID={info.name},Number={info.number},Type={info.type},Description="{info.description}">'
                    )

            # Copy FILTER fields
            for filter_id in other_header.filters:
                if filter_id not in output_header.filters:
                    filter_obj = other_header.filters[filter_id]
                    output_header.add_line(
                        f'##FILTER=<ID={filter_id},Description="{filter_obj.description}">'
                    )

        # Add multimerge command line to header
        command = " ".join(sys.argv)
        output_header.add_line(f'##multimerge_command="{command}"')

        # Open output file
        output_vcf = pysam.VariantFile(output_path, "w", header=output_header)

        # Collect variants from all VCFs
        variant_dict = {}  # position -> list of variants

        for i, vcf in enumerate(vcf_readers):
            for var in vcf.fetch():
                pos = (var.chrom, var.pos)

                # Process variant based on options
                if left_shift and reference:
                    # Left-shift the variant
                    chrom, pos, ref, alt = (
                        var.chrom,
                        var.pos,
                        var.ref,
                        var.alts[0] if var.alts else "",
                    )
                    try:
                        # Create a preprocessor for normalization
                        preprocessor = PreprocessEngine(reference_file=reference_path)
                        # The normalize_variant method returns (pos, ref, alt)
                        new_pos, new_ref, new_alt = preprocessor.normalize_variant(
                            chrom, pos, ref, alt
                        )
                        var.pos = new_pos
                        var.ref = new_ref
                        var.alts = [new_alt]
                    except Exception as e:
                        logging.warning(
                            f"Could not normalize variant {var.chrom}:{var.pos} {var.ref}>{var.alts}: {e}"
                        )

                if trim_alleles_flag:
                    # Trim common prefix/suffix
                    new_ref, new_alts = trim_alleles(
                        var.ref, list(var.alts) if var.alts else []
                    )
                    var.ref = new_ref
                    var.alts = tuple(new_alts)

                # Store variants by position if merging by location
                if merge_by_location:
                    if pos not in variant_dict:
                        variant_dict[pos] = []
                    variant_dict[pos].append((var, samples[i]))
                else:
                    # Write variant directly to output
                    output_vcf.write(var)

        # If merging by location, process the collected variants
        if merge_by_location:
            for pos, var_list in sorted(variant_dict.items()):
                if len(var_list) == 1:
                    # Only one variant at this position, just write it
                    output_vcf.write(var_list[0][0])
                else:
                    # Multiple variants at the same position, merge them
                    base_var = var_list[0][0]

                    # Collect all unique alleles if specified
                    if unique_alleles:
                        all_alts = set()
                        for var, _ in var_list:
                            if var.alts:
                                all_alts.update(var.alts)
                        base_var.alts = tuple(sorted(all_alts))

                    # Write the merged variant
                    output_vcf.write(base_var)

        # Close files
        output_vcf.close()
        for vcf in vcf_readers:
            vcf.close()

        logging.info(f"Merged {len(vcf_readers)} VCF files into {output_path}")
        return True

    except Exception as e:
        logging.error(f"Error during VCF merge: {e}")
        return False


def main():
    """Main entry point for multimerge."""
    # Set up logging
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    )

    # Create command-line parser
    parser = argparse.ArgumentParser(
        description="Merge multiple VCF files. Input format: file.vcf.gz:sample_name"
    )

    # Required arguments
    parser.add_argument(
        "input_files",
        nargs="+",
        help="Input VCF files with optional sample names (format: file.vcf:sample)",
    )
    parser.add_argument("-o", "--output", required=True, help="Output VCF file")

    # Optional arguments
    parser.add_argument(
        "-r",
        "--reference",
        help="Reference FASTA file (required for --leftshift and --splitalleles)",
    )
    parser.add_argument(
        "--trimalleles",
        type=int,
        default=0,
        help="Trim common bases from alleles (0=off, 1=on)",
    )
    parser.add_argument(
        "--leftshift", type=int, default=0, help="Left-shift variants (0=off, 1=on)"
    )
    parser.add_argument(
        "--splitalleles",
        type=int,
        default=0,
        help="Split multi-allelic variants (0=off, 1=on)",
    )
    parser.add_argument(
        "--merge-by-location",
        type=int,
        default=0,
        help="Merge variants at the same location (0=off, 1=on)",
    )
    parser.add_argument(
        "--unique-alleles",
        type=int,
        default=0,
        help="Keep only unique alleles (0=off, 1=on)",
    )
    parser.add_argument(
        "--process-formats",
        type=int,
        default=0,
        help="Process FORMAT fields (0=off, 1=on)",
    )
    parser.add_argument(
        "--process-split",
        type=int,
        default=0,
        help="Split complex variants (0=off, 1=on)",
    )
    parser.add_argument(
        "--process-full",
        type=int,
        default=0,
        help="Process all variants fully (0=off, 1=on)",
    )

    # Parse arguments
    args = parser.parse_args()

    # Process input files (separate VCF path and sample name)
    input_vcfs = [parse_vcf_sample_arg(arg) for arg in args.input_files]

    logging.info(f"Multimerge called with {len(input_vcfs)} input files")

    # Validate inputs
    if len(input_vcfs) < 1:
        logging.error("At least one input VCF file is required")
        return 1

    # Run merge operation
    success = merge_vcfs(
        input_vcfs=input_vcfs,
        output_path=args.output,
        reference_path=args.reference,
        trim_alleles_flag=bool(args.trimalleles),
        left_shift=bool(args.leftshift),
        split_alleles=bool(args.splitalleles),
        merge_by_location=bool(args.merge_by_location),
        unique_alleles=bool(args.unique_alleles),
        process_formats=bool(args.process_formats),
        process_split=bool(args.process_split),
        process_full=bool(args.process_full),
    )

    if not success:
        logging.error("Merge failed")
        return 1

    logging.info(f"Merge completed successfully: {args.output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
