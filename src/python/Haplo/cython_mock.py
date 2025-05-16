"""
Mock implementations of Cython modules for testing.

This module provides pure Python implementations of the Cython modules
for testing and development purposes. They can be used when the actual
Cython modules are not available or during development.
"""
import numpy as np
from functools import cmp_to_key

# Mock implementation of sequence_utils.pyx functions
def complement_sequence(seq_input):
    """
    Return the complement of a DNA sequence.
    
    Args:
        seq_input: String or bytes sequence to complement
        
    Returns:
        Same type as input: Complemented sequence
    """
    # Check if input is bytes or string
    is_bytes = isinstance(seq_input, bytes)
    
    # Convert to string if bytes
    if is_bytes:
        seq = seq_input.decode('ascii')
    else:
        seq = seq_input
        
    # Create complement mapping
    comp_dict = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N', 
                 'a': 't', 't': 'a', 'g': 'c', 'c': 'g', 'n': 'n'}
                 
    # Create complement
    result = ''.join([comp_dict.get(base, base) for base in seq])
    
    # Convert back to bytes if input was bytes
    if is_bytes:
        return result.encode('ascii')
    return result

def reverse_complement(seq_input):
    """
    Return the reverse complement of a DNA sequence.
    
    Args:
        seq_input: String or bytes sequence to reverse complement
        
    Returns:
        Same type as input: Reverse complemented sequence
    """
    # First get complement
    comp_seq = complement_sequence(seq_input)
    
    # Then reverse it
    if isinstance(comp_seq, bytes):
        return comp_seq[::-1]
    else:
        return comp_seq[::-1]

# Mock implementation of happyroc.pyx functions
def compute_roc_points(tp, fp, fn, tn=None, resolution=100):
    """
    Compute ROC curve points.
    
    Args:
        tp: True positives array or list
        fp: False positives array or list
        fn: False negatives array or list
        tn: True negatives array or list (optional)
        resolution: Number of points in the ROC curve
    
    Returns:
        tuple: (recalls, precisions, f1_scores, specificities)
    """
    # Convert inputs to numpy arrays if they aren't already
    tp = np.asarray(tp)
    fp = np.asarray(fp)
    fn = np.asarray(fn)
    if tn is not None:
        tn = np.asarray(tn)
    
    # Calculate metrics
    with np.errstate(divide='ignore', invalid='ignore'):
        recalls = tp / (tp + fn)
        precisions = tp / (tp + fp)
        
    # Handle divide by zero
    recalls = np.nan_to_num(recalls)
    precisions = np.nan_to_num(precisions)
    
    # Calculate F1 scores
    f1_scores = 2 * precisions * recalls / (precisions + recalls)
    f1_scores = np.nan_to_num(f1_scores)
    
    # Calculate specificities if tn is provided
    if tn is not None:
        specificities = tn / (tn + fp)
        specificities = np.nan_to_num(specificities)
    else:
        specificities = np.zeros_like(recalls)
    
    return recalls, precisions, f1_scores, specificities

# Mock implementation of variant_processor.pyx classes/functions
class VariantProcessor:
    """Mock implementation of the VariantProcessor Cython class."""
    
    def __init__(self):
        """Initialize an empty variant processor."""
        self.variants = []
        
    def add_variant(self, variant):
        """
        Add a variant to the processor.
        
        Args:
            variant: Variant object with chrom, pos, ref, alt attributes
        """
        self.variants.append(variant)
        
    def get_variant_chrom(self, idx):
        """
        Get chromosome name for a variant.
        
        Args:
            idx: Index of the variant
            
        Returns:
            str: Chromosome name
        """
        if idx >= len(self.variants):
            raise IndexError(f"Index {idx} out of range")
        return self.variants[idx].chrom
        
    def process_variants(self, threads=1):
        """
        Process variants with mock implementation.
        
        Args:
            threads: Number of threads to use (ignored in mock)
            
        Returns:
            list: Processed variant data
        """
        results = []
        for variant in self.variants:
            # Create a simple result dictionary
            result = {
                'chrom': variant.chrom,
                'position': variant.pos,
                'ref': variant.ref,
                'alt': variant.alt,
                'processed': True
            }
            results.append(result)
        return results

# Mock chromosome sorting implementation compatible with Python 3
def cmp_chromosomes(a, b):
    """
    Compare two chromosomes for sorting.
    
    Args:
        a: First chromosome name
        b: Second chromosome name
        
    Returns:
        int: -1 if a < b, 0 if a == b, 1 if a > b
    """
    # Handle 'chr' prefix
    a_chr = a[3:] if a.startswith('chr') else a
    b_chr = b[3:] if b.startswith('chr') else b
    
    # Handle numeric chromosomes
    if a_chr.isdigit() and b_chr.isdigit():
        return int(a_chr) - int(b_chr)
    elif a_chr.isdigit():
        return -1
    elif b_chr.isdigit():
        return 1
    else:
        # Special chromosomes like X, Y, M
        special_order = {'X': 1, 'Y': 2, 'M': 3, 'MT': 3}
        a_val = special_order.get(a_chr, 99)
        b_val = special_order.get(b_chr, 99)
        
        if a_val != 99 and b_val != 99:
            return a_val - b_val
        elif a_val != 99:
            return -1
        elif b_val != 99:
            return 1
        else:
            return -1 if a < b else (1 if a > b else 0)

def sort_chromosomes(chrom_list):
    """
    Sort a list of chromosomes in a standard order.
    
    Args:
        chrom_list: List of chromosome names
        
    Returns:
        list: Sorted list of chromosome names
    """
    return sorted(chrom_list, key=cmp_to_key(cmp_chromosomes))