#
# Copyright (c) 2010-2015 Illumina, Inc.
# All rights reserved.
#
# This file is distributed under the simplified BSD license.
# The full text can be found here (and in LICENSE.txt in the root folder of
# this distribution):
#
# https://github.com/Illumina/licenses/blob/master/Simplified-BSD-License.txt
#
# 10/02/2017
#
# Run gvcf2bed - Python implementation to replace C++ binary

import logging
import tempfile
from typing import Optional

import pysam


def gvcf2bed(vcf: str, ref: str, regions: Optional[str] = None, scratch_prefix: Optional[str] = None) -> str:
    """Convert GVCF/VCF to BED regions representing confident call regions.
    
    This Python implementation replaces the original C++ gvcf2bed binary.
    It identifies confident regions in a VCF/GVCF file and outputs them as BED format.
    
    Args:
        vcf: Path to input VCF/GVCF file
        ref: Path to reference FASTA file (for compatibility, not used in this implementation)
        regions: Optional target regions (BED format file)
        scratch_prefix: Directory for temporary files
        
    Returns:
        Path to temporary BED file containing confident regions
    """
    logging.info(f"Converting VCF to BED regions: {vcf}")
    
    # Create temporary output file
    tf = tempfile.NamedTemporaryFile(dir=scratch_prefix, suffix=".bed", delete=False)
    output_path = tf.name
    tf.close()
    
    try:
        # Open VCF file with pysam
        vcf_file = pysam.VariantFile(vcf)
        
        # Parse target regions if provided
        target_intervals = None
        if regions:
            target_intervals = []
            try:
                import gzip
                # Handle both gzipped and regular BED files
                if regions.endswith('.gz'):
                    open_func = gzip.open
                    mode = 'rt'
                else:
                    open_func = open
                    mode = 'r'
                
                with open_func(regions, mode) as reg_file:
                    for line in reg_file:
                        line = line.strip()
                        if line and not line.startswith('#'):
                            parts = line.split('\t')
                            if len(parts) >= 3:
                                chrom, start, end = parts[0], int(parts[1]), int(parts[2])
                                target_intervals.append((chrom, start, end))
            except Exception as e:
                logging.warning(f"Could not parse regions file {regions}: {e}")
                target_intervals = None
        
        confident_regions = []
        
        # Process VCF records to identify confident regions
        for record in vcf_file:
            chrom = record.chrom
            start = record.start  # 0-based
            stop = record.stop   # 1-based, so this is the end position
            
            # Skip if target regions specified and this region is not in target
            if target_intervals:
                in_target = False
                for t_chrom, t_start, t_end in target_intervals:
                    if (chrom == t_chrom and 
                        not (stop <= t_start or start >= t_end)):  # Overlaps target
                        in_target = True
                        break
                if not in_target:
                    continue
            
            # Check if this is a confident call
            # For simplicity, we consider a region confident if:
            # 1. It has PASS filter or no filter
            # 2. It's not a missing/no-call genotype
            is_confident = True
            
            # Check filters
            if record.filter.keys():
                if 'PASS' not in record.filter.keys():
                    is_confident = False
            
            # Check genotypes for confident calls
            if is_confident and record.samples:
                for sample in record.samples.values():
                    gt = sample.get('GT', None)
                    if gt is None or None in gt:  # Missing or partial genotype
                        is_confident = False
                        break
            
            if is_confident:
                confident_regions.append((chrom, start, stop))
        
        vcf_file.close()
        
        # Merge overlapping regions and write to BED file
        if confident_regions:
            confident_regions.sort()
            merged_regions = []
            current_chrom, current_start, current_end = confident_regions[0]
            
            for chrom, start, end in confident_regions[1:]:
                if chrom == current_chrom and start <= current_end:
                    # Overlapping or adjacent - merge
                    current_end = max(current_end, end)
                else:
                    # Non-overlapping - save current and start new
                    merged_regions.append((current_chrom, current_start, current_end))
                    current_chrom, current_start, current_end = chrom, start, end
            
            # Add the last region
            merged_regions.append((current_chrom, current_start, current_end))
        else:
            merged_regions = []
        
        # Write BED file
        with open(output_path, 'w') as bed_file:
            for chrom, start, end in merged_regions:
                # BED format is 0-based start, 1-based end (exclusive)
                bed_file.write(f"{chrom}\t{start}\t{end}\n")
        
        logging.info(f"Generated {len(merged_regions)} confident regions in {output_path}")
        
    except Exception as e:
        logging.error(f"Error processing VCF file {vcf}: {e}")
        # Create an empty BED file as fallback
        with open(output_path, 'w') as bed_file:
            pass
    
    return output_path
