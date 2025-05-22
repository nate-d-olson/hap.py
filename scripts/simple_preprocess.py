#!/usr/bin/env python3
"""
Simplified standalone preprocess demo.

This script demonstrates a simple version of VCF preprocessing.
"""

import argparse
import logging
import os
import sys
import tempfile

try:
    import pysam
except ImportError:
    print("Error: pysam is required. Install it with 'pip install pysam'")
    sys.exit(1)

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def preprocess_vcf(input_vcf, reference_fasta, output_vcf=None, decompose=True, 
                  left_shift=True, pass_only=False):
    """
    Preprocess a VCF file.
    
    Args:
        input_vcf: Path to input VCF file
        reference_fasta: Path to reference FASTA file
        output_vcf: Path to output VCF file
        decompose: Whether to decompose multi-allelic variants
        left_shift: Whether to left-align variants
        pass_only: Whether to only include PASS variants
    
    Returns:
        Path to output VCF file
    """
    logger.info(f"Processing {input_vcf}")
    logger.info(f"Using reference {reference_fasta}")
    
    # Create output file if not provided
    if not output_vcf:
        fd, output_vcf = tempfile.mkstemp(suffix=".vcf.gz")
        os.close(fd)
    
    logger.info(f"Output will be written to {output_vcf}")
    
    # Open input VCF
    try:
        vcf_in = pysam.VariantFile(input_vcf)
        logger.info(f"Successfully opened input VCF with {sum(1 for _ in vcf_in)} variants")
        # Reset the file position
        vcf_in = pysam.VariantFile(input_vcf)
    except Exception as e:
        logger.error(f"Failed to open input VCF: {e}")
        raise
    
    # Open reference
    reference = pysam.FastaFile(reference_fasta)
    
    # Create output header
    header = vcf_in.header.copy()
    
    # Add preprocessing info to header
    header.add_line('##INFO=<ID=DECOMPOSED,Number=0,Type=Flag,'
                  'Description="Variant was decomposed from a multi-allelic variant">')
    header.add_line('##INFO=<ID=LEFTSHIFTED,Number=0,Type=Flag,'
                  'Description="Variant was left-shifted">')
    
    # Create output file
    vcf_out = pysam.VariantFile(output_vcf, mode="w", header=header)
    
    # Statistics
    stats = {
        "total": 0,
        "decomposed": 0,
        "left_shifted": 0,
        "skipped": 0
    }
    
    # Process variants
    for record in vcf_in:
        stats["total"] += 1
        
        # Skip non-PASS variants if requested
        if pass_only and not (len(record.filter) == 0 or "PASS" in record.filter):
            stats["skipped"] += 1
            continue
        
        # Process multi-allelic variants
        if decompose and len(record.alts) > 1:
            stats["decomposed"] += 1
            
            # Create a separate record for each alt allele
            for i, alt in enumerate(record.alts):
                new_record = vcf_out.new_record()
                new_record.contig = record.contig
                new_record.pos = record.pos
                new_record.id = record.id
                new_record.ref = record.ref
                new_record.alts = (alt,)
                
                # Copy filters
                for f in record.filter:
                    new_record.filter.add(f)
                
                # Copy sample information
                for sample in record.samples:
                    for field in record.format:
                        if field == "GT":
                            # Adjust genotype for bi-allelic variant
                            gt = record.samples[sample]["GT"]
                            new_gt = []
                            for g in gt:
                                if g is None:
                                    new_gt.append(None)
                                elif g == 0:
                                    new_gt.append(0)  # Reference
                                elif g == i + 1:
                                    new_gt.append(1)  # This alt
                                else:
                                    new_gt.append(0)  # Other alt
                            new_record.samples[sample]["GT"] = tuple(new_gt)
                        else:
                            # Copy other fields
                            try:
                                new_record.samples[sample][field] = record.samples[sample][field]
                            except:
                                pass
                
                # Left-shift if requested
                if left_shift:
                    # Basic left-shifting implementation
                    if len(new_record.ref) != len(alt):  # Only for indels
                        try:
                            # Simple check - try to shift one base left
                            if len(new_record.ref) > len(alt):  # Deletion
                                base_to_check = reference.fetch(
                                    new_record.contig, 
                                    new_record.pos - 2, 
                                    new_record.pos - 1
                                )
                                if base_to_check == new_record.ref[-1]:
                                    # Can shift left
                                    new_record.pos -= 1
                                    new_record.ref = base_to_check + new_record.ref[:-1]
                                    new_record.alts = (base_to_check + alt,)
                                    new_record.info["LEFTSHIFTED"] = True
                                    stats["left_shifted"] += 1
                            elif len(new_record.ref) < len(alt):  # Insertion
                                base_to_check = reference.fetch(
                                    new_record.contig, 
                                    new_record.pos - 2, 
                                    new_record.pos - 1
                                )
                                if base_to_check == alt[-1]:
                                    # Can shift left
                                    new_record.pos -= 1
                                    new_record.ref = base_to_check + new_record.ref
                                    new_record.alts = (base_to_check + alt[:-1],)
                                    new_record.info["LEFTSHIFTED"] = True
                                    stats["left_shifted"] += 1
                        except:
                            # Ignore errors in left-shifting
                            pass
                
                # Mark as decomposed
                new_record.info["DECOMPOSED"] = True
                
                # Write to output
                vcf_out.write(new_record)
        else:
            # Single-allelic variant
            if left_shift and len(record.ref) != len(record.alts[0]):
                # Basic left-shifting implementation for indels
                try:
                    # Simple check - try to shift one base left
                    if len(record.ref) > len(record.alts[0]):  # Deletion
                        base_to_check = reference.fetch(
                            record.contig, 
                            record.pos - 2, 
                            record.pos - 1
                        )
                        if base_to_check == record.ref[-1]:
                            # Can shift left
                            record.pos -= 1
                            record.ref = base_to_check + record.ref[:-1]
                            record.alts = (base_to_check + record.alts[0],)
                            record.info["LEFTSHIFTED"] = True
                            stats["left_shifted"] += 1
                    elif len(record.ref) < len(record.alts[0]):  # Insertion
                        base_to_check = reference.fetch(
                            record.contig, 
                            record.pos - 2, 
                            record.pos - 1
                        )
                        if base_to_check == record.alts[0][-1]:
                            # Can shift left
                            record.pos -= 1
                            record.ref = base_to_check + record.ref
                            record.alts = (base_to_check + record.alts[0][:-1],)
                            record.info["LEFTSHIFTED"] = True
                            stats["left_shifted"] += 1
                except:
                    # Ignore errors in left-shifting
                    pass
            
            # Write to output
            vcf_out.write(record)
    
    # Close files
    vcf_in.close()
    vcf_out.close()
    reference.close()
    
    # Index output if it's compressed
    if output_vcf.endswith(".vcf.gz"):
        try:
            pysam.tabix_index(output_vcf, preset="vcf")
        except:
            logger.warning("Failed to index output VCF")
    
    # Print statistics
    logger.info("Preprocessing statistics:")
    logger.info(f"  Total variants: {stats['total']}")
    logger.info(f"  Decomposed multi-allelic variants: {stats['decomposed']}")
    logger.info(f"  Left-shifted variants: {stats['left_shifted']}")
    logger.info(f"  Skipped variants: {stats['skipped']}")
    
    return output_vcf


def main():
    """Command line interface."""
    parser = argparse.ArgumentParser(description="Preprocess VCF files")
    
    parser.add_argument("input_vcf", help="Input VCF file")
    parser.add_argument("-r", "--reference", required=True, 
                       help="Reference FASTA file")
    parser.add_argument("-o", "--output", help="Output VCF file")
    parser.add_argument("--no-decompose", action="store_false", dest="decompose",
                       help="Don't decompose multi-allelic variants")
    parser.add_argument("--no-left-shift", action="store_false", dest="left_shift",
                       help="Don't left-shift variants")
    parser.add_argument("--pass-only", action="store_true",
                       help="Only include PASS variants")
    parser.add_argument("-v", "--verbose", action="store_true",
                       help="Enable verbose logging")
    
    args = parser.parse_args()
    
    # Set logging level
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    try:
        output_file = preprocess_vcf(
            args.input_vcf,
            args.reference,
            args.output,
            args.decompose,
            args.left_shift,
            args.pass_only
        )
        print(f"Preprocessing complete. Output written to {output_file}")
        return 0
    except Exception as e:
        logger.error(f"Error: {e}")
        return 1


if __name__ == "__main__":
    sys.exit(main())
