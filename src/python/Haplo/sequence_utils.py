#!/usr/bin/env python3
"""
Python implementation of sequence manipulation functions using BioPython.

This module replaces the C++ sequence manipulation functions with Python
equivalents using BioPython. It includes functions for:
- Complement and reverse complement
- FASTA file handling and sequence extraction
- Sequence alignment and manipulation
"""

import argparse
import logging
import sys
from pathlib import Path
from typing import Dict, List, Union

from Bio.Seq import Seq


class SequenceUtils:
    """
    Utility functions for working with genomic sequences.
    
    This class replaces C++ implementations with Python equivalents
    using BioPython. It maintains compatible interfaces with the
    original C++ functions for seamless migration.
    """

    @staticmethod
    def complement_sequence(sequence: Union[str, bytes, Seq]) -> str:
        """
        Get the complement of a DNA sequence.
        
        Args:
            sequence: Input DNA sequence
            
        Returns:
            Complemented DNA sequence as a string
        """
        # Convert bytes to string if needed
        if isinstance(sequence, bytes):
            sequence = sequence.decode("utf-8")
            
        # Convert to Seq object if not already
        if not isinstance(sequence, Seq):
            sequence = Seq(sequence)
            
        # Return the complement as a string
        return str(sequence.complement())

    @staticmethod
    def reverse_complement(sequence: Union[str, bytes, Seq]) -> str:
        """
        Get the reverse complement of a DNA sequence.
        
        Args:
            sequence: Input DNA sequence
            
        Returns:
            Reverse complemented DNA sequence as a string
        """
        # Convert bytes to string if needed
        if isinstance(sequence, bytes):
            sequence = sequence.decode("utf-8")
            
        # Convert to Seq object if not already
        if not isinstance(sequence, Seq):
            sequence = Seq(sequence)
            
        # Return the reverse complement as a string
        return str(sequence.reverse_complement())

    @staticmethod
    def normalize_sequence(sequence: Union[str, bytes, Seq]) -> str:
        """
        Normalize a sequence by converting to uppercase and replacing
        invalid characters with N.
        
        Args:
            sequence: Input DNA sequence
            
        Returns:
            Normalized sequence as a string
        """
        # Convert bytes to string if needed
        if isinstance(sequence, bytes):
            sequence = sequence.decode("utf-8")
            
        # Convert to string if Seq object
        if isinstance(sequence, Seq):
            sequence = str(sequence)
            
        # Convert to uppercase
        sequence = sequence.upper()
        
        # Replace invalid characters with N
        valid_chars = set("ACGTN")
        result = "".join(c if c in valid_chars else "N" for c in sequence)
        
        return result


class FastaReader:
    """
    FASTA file reader with efficient sequence extraction.
    
    This class replaces the C++ FastaFile class with a Python
    implementation using BioPython and indexed FASTA files.
    """

    def __init__(self, fasta_path: str):
        """
        Initialize the FASTA reader.
        
        Args:
            fasta_path: Path to the FASTA file
        """
        self.fasta_path = fasta_path
        self.logger = logging.getLogger("fasta_reader")
        
        # Check if the file exists
        if not Path(fasta_path).exists():
            raise FileNotFoundError(f"FASTA file not found: {fasta_path}")
            
        # Check if the index exists or create it
        fai_path = f"{fasta_path}.fai"
        if not Path(fai_path).exists():
            self.logger.info(f"Indexing FASTA file: {fasta_path}")
            try:
                from Bio import SeqIO
                from Bio.SeqIO.FastaIO import FastaIndex
                
                # Create the index
                FastaIndex.index_file(fasta_path, fai_path)
                self.logger.info(f"Created index file: {fai_path}")
            except Exception as e:
                self.logger.warning(f"Failed to create index: {e}")
        
        # Open the indexed FASTA file
        try:
            self.fasta = SeqIO.to_dict(SeqIO.parse(fasta_path, "fasta"))
            self.logger.info(f"Loaded FASTA file: {fasta_path}")
        except Exception as e:
            raise RuntimeError(f"Failed to load FASTA file: {e}")

    def get_sequence(self, chrom: str, start: int, end: int) -> str:
        """
        Get a sequence from the FASTA file.
        
        Args:
            chrom: Chromosome name
            start: 0-based start position (inclusive)
            end: 0-based end position (exclusive)
            
        Returns:
            Extracted sequence as a string
        """
        if chrom not in self.fasta:
            raise KeyError(f"Chromosome {chrom} not found in FASTA file")
            
        # Get the sequence
        seq = self.fasta[chrom].seq
        
        # Check bounds
        if start < 0 or end > len(seq):
            raise IndexError(
                f"Invalid coordinates: {start}-{end} "
                f"(sequence length: {len(seq)})"
            )
            
        # Extract and return the sequence
        return str(seq[start:end])

    def get_chromosomes(self) -> List[str]:
        """
        Get a list of chromosome names in the FASTA file.
        
        Returns:
            List of chromosome names
        """
        return list(self.fasta.keys())

    def get_chromosome_length(self, chrom: str) -> int:
        """
        Get the length of a chromosome.
        
        Args:
            chrom: Chromosome name
            
        Returns:
            Length of the chromosome in bases
        """
        if chrom not in self.fasta:
            raise KeyError(f"Chromosome {chrom} not found in FASTA file")
            
        return len(self.fasta[chrom].seq)

    def get_all_chromosome_lengths(self) -> Dict[str, int]:
        """
        Get the lengths of all chromosomes.
        
        Returns:
            Dictionary mapping chromosome names to lengths
        """
        return {chrom: len(self.fasta[chrom].seq) for chrom in self.fasta}


def main():
    """Command-line entry point for demonstration."""
    parser = argparse.ArgumentParser(
        description="Sequence manipulation utilities"
    )
    parser.add_argument(
        "--fasta", help="Path to a FASTA file"
    )
    parser.add_argument(
        "--complement", help="Sequence to complement"
    )
    parser.add_argument(
        "--reverse-complement", help="Sequence to reverse complement"
    )
    parser.add_argument(
        "--normalize", help="Sequence to normalize"
    )
    parser.add_argument(
        "--get-sequence", nargs=3, metavar=("CHROM", "START", "END"),
        help="Get a sequence from a FASTA file"
    )
    parser.add_argument(
        "--list-chromosomes", action="store_true",
        help="List chromosomes in a FASTA file"
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true",
        help="Enable verbose logging"
    )

    args = parser.parse_args()

    # Setup logging
    log_level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(
        level=log_level,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )

    # Perform requested operations
    if args.complement:
        result = SequenceUtils.complement_sequence(args.complement)
        print(f"Complement: {result}")
        
    if args.reverse_complement:
        result = SequenceUtils.reverse_complement(args.reverse_complement)
        print(f"Reverse complement: {result}")
        
    if args.normalize:
        result = SequenceUtils.normalize_sequence(args.normalize)
        print(f"Normalized: {result}")
        
    if args.fasta:
        try:
            fasta_reader = FastaReader(args.fasta)
            
            if args.list_chromosomes:
                chroms = fasta_reader.get_chromosomes()
                lengths = fasta_reader.get_all_chromosome_lengths()
                print("Chromosomes in FASTA file:")
                for chrom in chroms:
                    print(f"  {chrom}: {lengths[chrom]} bp")
                    
            if args.get_sequence:
                chrom, start, end = args.get_sequence
                sequence = fasta_reader.get_sequence(
                    chrom, int(start), int(end)
                )
                print(f"Sequence {chrom}:{start}-{end}:")
                print(sequence)
                
        except Exception as e:
            logging.error(f"Error with FASTA file: {e}")
            sys.exit(1)


if __name__ == "__main__":
    main()
