#!/usr/bin/env python3
"""
Python implementation of the blocksplit functionality.

This module replaces the C++ blocksplit tool, which divides a VCF file
into blocks of variants for more efficient processing. The Python version
uses pysam to achieve the same functionality while being more maintainable.
"""

import argparse
import logging
import sys
from typing import Dict, Iterator, List, Optional

from pysam import VariantFile, VariantRecord


class BlockSplitter:
    """Splits a VCF file into blocks of variants."""

    def __init__(
        self,
        block_size: int = 1000,
        min_distance: int = 1000,
        apply_filters: bool = False,
    ):
        """
        Initialize a BlockSplitter instance.

        Args:
            block_size: Maximum number of variants per block
            min_distance: Minimum distance between variants to start a new block
            apply_filters: Whether to apply filters (skip filtered variants)
        """
        self.block_size = block_size
        self.min_distance = min_distance
        self.apply_filters = apply_filters
        self.logger = logging.getLogger("blocksplit")

    def process_file(
        self, vcf_path: str, output_bed: Optional[str] = None
    ) -> List[Dict[str, any]]:
        """
        Process a VCF file and split it into blocks.

        Args:
            vcf_path: Path to the input VCF file
            output_bed: Optional path to write blocks to BED file

        Returns:
            List of block dictionaries with chromosome, start, end positions
        """
        self.logger.info(f"Processing {vcf_path}")
        blocks = list(self._generate_blocks(vcf_path))
        
        if output_bed:
            self._write_bed(blocks, output_bed)
            
        return blocks

    def _generate_blocks(self, vcf_path: str) -> Iterator[Dict[str, any]]:
        """
        Generate blocks from a VCF file.

        Args:
            vcf_path: Path to the input VCF file

        Yields:
            Block dictionaries with chromosome, start, end positions
        """
        try:
            vcf = VariantFile(vcf_path)
        except Exception as e:
            self.logger.error(f"Failed to open VCF file: {e}")
            return

        current_block: List[VariantRecord] = []
        previous_pos = -1
        current_chrom = ""

        for record in vcf:
            # Skip filtered records if requested
            if self.apply_filters and self._is_filtered(record):
                continue

            # Get chromosome and position
            chrom = record.chrom
            pos = record.pos

            # Handle chromosome change
            if chrom != current_chrom:
                if current_block:
                    yield self._create_block_dict(current_chrom, current_block)
                    current_block = []
                current_chrom = chrom
                previous_pos = -1

            # Check if we need to start a new block
            dist_to_previous = 0 if previous_pos < 0 else pos - previous_pos
            if (len(current_block) >= self.block_size or 
                    dist_to_previous > self.min_distance):
                if current_block:
                    yield self._create_block_dict(current_chrom, current_block)
                    current_block = []

            # Add record to current block
            current_block.append(record)
            previous_pos = pos

        # Handle the last block
        if current_block:
            yield self._create_block_dict(current_chrom, current_block)

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

    def _create_block_dict(
        self, chrom: str, records: List[VariantRecord]
    ) -> Dict[str, any]:
        """
        Create a block dictionary from a list of records.

        Args:
            chrom: Chromosome name
            records: List of variant records in the block

        Returns:
            Dictionary with block information
        """
        start = records[0].pos
        end = records[-1].pos
        
        return {
            "chrom": chrom,
            "start": start,
            "end": end,
            "count": len(records),
            "span": end - start + 1,
        }

    def _write_bed(self, blocks: List[Dict[str, any]], output_path: str) -> None:
        """
        Write blocks to a BED file.

        Args:
            blocks: List of block dictionaries
            output_path: Path to write the BED file
        """
        with open(output_path, "w") as f:
            for block in blocks:
                f.write(
                    f"{block['chrom']}\t{block['start']}\t{block['end']}\t"
                    f"variants={block['count']};span={block['span']}\n"
                )


def main():
    """Command-line entry point."""
    parser = argparse.ArgumentParser(
        description="Split a VCF file into blocks of variants"
    )
    parser.add_argument(
        "input_file", help="Input VCF/BCF file to process"
    )
    parser.add_argument(
        "-o", "--output-bed", help="Output BED file with block regions"
    )
    parser.add_argument(
        "-b", "--block-size", type=int, default=1000,
        help="Maximum number of variants per block"
    )
    parser.add_argument(
        "-d", "--min-distance", type=int, default=1000,
        help="Minimum distance between variants to start a new block"
    )
    parser.add_argument(
        "-f", "--apply-filters", action="store_true",
        help="Skip filtered variants (those with FILTER != PASS)"
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

    # Create and run the block splitter
    splitter = BlockSplitter(
        block_size=args.block_size,
        min_distance=args.min_distance,
        apply_filters=args.apply_filters,
    )
    
    try:
        blocks = splitter.process_file(args.input_file, args.output_bed)
        logging.info(f"Created {len(blocks)} blocks")
    except Exception as e:
        logging.error(f"Error processing file: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
