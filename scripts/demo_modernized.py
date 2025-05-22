#!/usr/bin/env python3
"""
Demonstration script for the modernized hap.py components.

This script shows how to use the Python implementations of 
blocksplit, vcfcheck, sequence utilities, and quantify.
"""

import argparse
import os
import sys
from pathlib import Path

# Add the src/python directory to the path so we can import from Haplo
sys.path.insert(0, str(Path(__file__).parent.parent / "src" / "python"))

from Haplo.python_blocksplit import BlockSplitter
from Haplo.python_preprocess import PreprocessEngine
from Haplo.python_quantify import QuantifyEngine
from Haplo.python_vcfcheck import VCFChecker
from Haplo.sequence_utils import FastaReader, SequenceUtils


def demo_blocksplit(vcf_path: str, output_path: str = None, block_size: int = 1000):
    """
    Demonstrate the BlockSplitter functionality.
    
    Args:
        vcf_path: Path to input VCF file
        output_path: Path to output BED file
        block_size: Block size to use
    """
    print("\n=== BlockSplit Demo ===\n")
    print(f"Input VCF: {vcf_path}")
    print(f"Block size: {block_size}")
    
    # Create the BlockSplitter
    splitter = BlockSplitter(block_size=block_size)
    
    # Process the VCF file
    blocks = splitter.process_file(vcf_path, output_path)
    
    # Print summary
    print(f"\nGenerated {len(blocks)} blocks")
    
    # Print a few blocks as example
    print("\nExample blocks:")
    for i, block in enumerate(blocks[:5]):
        print(f"  Block {i+1}: {block['chrom']}:{block['start']}-{block['end']} "
              f"({block['count']} variants, {block['span']} bp)")
    
    if output_path:
        print(f"\nBlocks written to: {output_path}")


def demo_vcfcheck(vcf_path: str, reference_path: str = None, output_path: str = None):
    """
    Demonstrate the VCFChecker functionality.
    
    Args:
        vcf_path: Path to input VCF file
        reference_path: Path to reference FASTA file
        output_path: Path to output issues file
    """
    print("\n=== VCFCheck Demo ===\n")
    print(f"Input VCF: {vcf_path}")
    if reference_path:
        print(f"Reference: {reference_path}")
    
    # Create the VCFChecker
    checker = VCFChecker(reference_path=reference_path)
    
    # Check the VCF file
    stats = checker.check_file(vcf_path, output_path)
    
    # Print summary
    checker.print_summary()
    
    if output_path:
        print(f"\nIssues written to: {output_path}")


def demo_quantify(truth_vcf: str, query_vcf: str, reference_path: str = None, output_prefix: str = None):
    """
    Demonstrate the QuantifyEngine functionality.
    
    Args:
        truth_vcf: Path to truth VCF file
        query_vcf: Path to query VCF file
        reference_path: Path to reference FASTA file
        output_prefix: Prefix for output files
    """
    print("\n=== Quantify Demo ===\n")
    print(f"Truth VCF: {truth_vcf}")
    print(f"Query VCF: {query_vcf}")
    if reference_path:
        print(f"Reference: {reference_path}")
    
    # Create the QuantifyEngine
    engine = QuantifyEngine(
        truth_vcf=truth_vcf,
        query_vcf=query_vcf,
        reference=reference_path,
        output_vtc=True
    )
    
    # Run quantification
    results = engine.quantify()
    
    # Print summary
    print("\nQuantification Results:")
    print(f"  Total TP: {results['metrics']['TP']}")
    print(f"  Total FP: {results['metrics']['FP']}")
    print(f"  Total FN: {results['metrics']['FN']}")
    print(f"  Precision: {results['metrics']['PRECISION']:.4f}")
    print(f"  Recall: {results['metrics']['RECALL']:.4f}")
    print(f"  F1 Score: {results['metrics']['F1']:.4f}")
    
    # Print stratifications
    print("\nStratifications:")
    for strat_type, strat_data in results['stratifications'].items():
        print(f"  {strat_type.upper()}:")
        for type_value, metrics in strat_data.items():
            print(f"    {type_value}: TP={metrics['TP']}, FP={metrics['FP']}, "
                  f"FN={metrics['FN']}, F1={metrics['F1']:.4f}")
    
    # Write results if output prefix is provided
    if output_prefix:
        engine.write_results(output_prefix)
        print(f"\nResults written to {output_prefix}.*.json/tsv")


def demo_sequence_utils(reference_path: str, chrom: str, start: int, end: int):
    """
    Demonstrate the sequence utilities functionality.
    
    Args:
        reference_path: Path to reference FASTA file
        chrom: Chromosome name
        start: Start position
        end: End position
    """
    print("\n=== Sequence Utilities Demo ===\n")
    print(f"Reference: {reference_path}")
    print(f"Region: {chrom}:{start}-{end}")
    
    # Create the FastaReader
    reader = FastaReader(reference_path)
    
    # Get the sequence
    sequence = reader.get_sequence(chrom, start, end)
    
    # Print the sequence
    print(f"\nSequence ({len(sequence)} bp):")
    print(sequence)
    
    # Demonstrate sequence operations
    print("\nComplement:")
    print(SequenceUtils.complement_sequence(sequence))
    
    print("\nReverse complement:")
    print(SequenceUtils.reverse_complement(sequence))
    
    print("\nNormalized sequence:")
    print(SequenceUtils.normalize_sequence(sequence))
    
    # Print chromosome information
    print("\nChromosome information:")
    print(f"  Available chromosomes: {', '.join(reader.get_chromosomes())}")
    print(f"  {chrom} length: {reader.get_chromosome_length(chrom)} bp")


def demo_preprocess(input_vcf: str, reference_path: str, output_path: str = None):
    """
    Demonstrate the PreprocessEngine functionality.
    
    Args:
        input_vcf: Path to input VCF file
        reference_path: Path to reference FASTA file
        output_path: Path to output VCF file
    """
    print("\n=== Preprocess Demo ===\n")
    print(f"Input VCF: {input_vcf}")
    print(f"Reference: {reference_path}")
    
    # Create the PreprocessEngine
    engine = PreprocessEngine(
        input_vcf=input_vcf,
        reference_fasta=reference_path,
        output_vcf=output_path,
        decompose_level=2,  # Aggressive decomposition
        left_shift=True
    )
    
    # Process the VCF file
    processed_vcf = engine.process()
    
    # Print summary
    print("\nPreprocessing Statistics:")
    print(f"  Total variants processed: {engine.stats['total_variants']}")
    print(f"  Decomposed variants: {engine.stats['decomposed_variants']}")
    print(f"  Left-shifted variants: {engine.stats['left_shifted_variants']}")
    print(f"  Normalized variants: {engine.stats['normalized_variants']}")
    
    if output_path:
        print(f"\nProcessed VCF written to: {processed_vcf}")
    
    return processed_vcf


def main():
    """Main function."""
    parser = argparse.ArgumentParser(
        description="Demonstrate the modernized hap.py components"
    )
    parser.add_argument(
        "--vcf", type=str,
        help="Path to input VCF file"
    )
    parser.add_argument(
        "--reference", type=str,
        help="Path to reference FASTA file"
    )
    parser.add_argument(
        "--truth-vcf", type=str,
        help="Path to truth VCF file for quantify"
    )
    parser.add_argument(
        "--query-vcf", type=str,
        help="Path to query VCF file for quantify"
    )
    parser.add_argument(
        "--blocksplit", action="store_true",
        help="Run BlockSplitter demo"
    )
    parser.add_argument(
        "--vcfcheck", action="store_true",
        help="Run VCFChecker demo"
    )
    parser.add_argument(
        "--sequence", action="store_true",
        help="Run sequence utilities demo"
    )
    parser.add_argument(
        "--quantify", action="store_true",
        help="Run quantify demo"
    )
    parser.add_argument(
        "--preprocess", action="store_true",
        help="Run preprocess demo"
    )
    parser.add_argument(
        "--preprocess", action="store_true",
        help="Run preprocess demo"
    )
    parser.add_argument(
        "--output-dir", type=str, default="demo_output",
        help="Directory for output files"
    )
    
    parser.add_argument(
        "--chrom", type=str, default="chr21",
        help="Chromosome for sequence demo"
    )
    parser.add_argument(
        "--start", type=int, default=1000,
        help="Start position for sequence demo"
    )
    parser.add_argument(
        "--end", type=int, default=1100,
        help="End position for sequence demo"
    )
    
    args = parser.parse_args()
    
    # Create output directory if it doesn't exist
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)
    
    # Run the requested demos
    if args.blocksplit and args.vcf:
        output_path = os.path.join(args.output_dir, "blocksplit_output.bed")
        demo_blocksplit(args.vcf, output_path)
    
    if args.vcfcheck and args.vcf:
        output_path = os.path.join(args.output_dir, "vcfcheck_output.txt")
        demo_vcfcheck(args.vcf, args.reference, output_path)
    
    if args.sequence and args.reference:
        demo_sequence_utils(args.reference, args.chrom, args.start, args.end)
    
    if args.quantify and args.truth_vcf and args.query_vcf:
        demo_quantify(args.truth_vcf, args.query_vcf, args.reference, args.output_dir)
    
    if args.preprocess and args.vcf and args.reference:
        output_path = os.path.join(args.output_dir, "preprocess_output.vcf")
        demo_preprocess(args.vcf, args.reference, output_path)
    
    # If no demos were requested, print help
    if not (args.blocksplit or args.vcfcheck or args.sequence or args.quantify or args.preprocess):
        parser.print_help()
        return 1
    
    return 0


if __name__ == "__main__":
    sys.exit(main())
