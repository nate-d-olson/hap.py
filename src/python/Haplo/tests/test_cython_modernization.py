#!/usr/bin/env python3
"""
Test script for Cython modernization and Python 3 compatibility.

This script verifies that the Cython modules work correctly with Python 3
and checks that the fallback implementations work as expected.
"""
import os
import sys
import argparse

# Add the parent directory to the path so we can import the Haplo module
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))

try:
    from Haplo import (
        complement_sequence,
        reverse_complement,
        VariantProcessor,
        compute_roc_points,
        cmp_chromosomes,
        sort_chromosomes,
        get_module_info
    )
except ImportError as e:
    print(f"Error importing Haplo module: {e}")
    sys.exit(1)

def test_sequence_functions():
    """Test the sequence-related functions."""
    print("\nTesting sequence functions...")
    
    # Test complement_sequence
    test_seq = "ACGTACGT"
    expected_comp = "TGCATGCA"
    actual_comp = complement_sequence(test_seq)
    
    if isinstance(actual_comp, bytes):
        actual_comp = actual_comp.decode('ascii')
        
    print(f"  complement_sequence('{test_seq}') = '{actual_comp}'")
    assert actual_comp == expected_comp, f"Expected '{expected_comp}', got '{actual_comp}'"
    
    # Test bytes input
    test_bytes = b"ACGTACGT"
    actual_comp_bytes = complement_sequence(test_bytes)
    
    if not isinstance(actual_comp_bytes, bytes):
        actual_comp_bytes = actual_comp_bytes.encode('ascii')
        
    print(f"  complement_sequence({test_bytes}) = {actual_comp_bytes}")
    assert actual_comp_bytes == b"TGCATGCA", f"Expected b'TGCATGCA', got {actual_comp_bytes}"
    
    # Test reverse_complement
    expected_revcomp = "ACGATGCA"
    actual_revcomp = reverse_complement("TGCATCGT")
    
    if isinstance(actual_revcomp, bytes):
        actual_revcomp = actual_revcomp.decode('ascii')
        
    print(f"  reverse_complement('TGCATCGT') = '{actual_revcomp}'")
    assert actual_revcomp == expected_revcomp, f"Expected '{expected_revcomp}', got '{actual_revcomp}'"
    
    print("Sequence functions test: PASSED")

def test_variant_processor():
    """Test the VariantProcessor class."""
    print("\nTesting VariantProcessor...")
    
    # Create a simple variant class for testing
    class SimpleVariant:
        def __init__(self, chrom, pos, ref, alt):
            self.chrom = chrom
            self.pos = pos
            self.ref = ref
            self.alt = alt
    
    # Create a processor
    processor = VariantProcessor()
    
    # Add some variants
    variants = [
        SimpleVariant("chr1", 100, "A", "T"),
        SimpleVariant("chr2", 200, "G", "C"),
        SimpleVariant("chrX", 300, "AT", "A")
    ]
    
    for v in variants:
        processor.add_variant(v)
        
    # Test get_variant_chrom
    chrom1 = processor.get_variant_chrom(0)
    print(f"  processor.get_variant_chrom(0) = '{chrom1}'")
    assert chrom1 == "chr1", f"Expected 'chr1', got '{chrom1}'"
    
    # Test processing
    results = processor.process_variants(threads=2)
    print(f"  processor.process_variants() returned {len(results)} results")
    assert len(results) == 3, f"Expected 3 results, got {len(results)}"
    
    # Check the first result
    if 'chrom' in results[0]:
        print(f"  results[0]['chrom'] = '{results[0]['chrom']}'")
        assert results[0]['chrom'] == "chr1", f"Expected 'chr1', got '{results[0]['chrom']}'"
        
    print("VariantProcessor test: PASSED")

def test_roc_functions():
    """Test the ROC-related functions."""
    print("\nTesting ROC functions...")
    
    # Create some test data
    tp = [10, 20, 30, 40, 50]
    fp = [5, 4, 3, 2, 1]
    fn = [1, 2, 3, 4, 5]
    
    # Compute ROC points
    recalls, precisions, f1_scores, specificities = compute_roc_points(tp, fp, fn)
    
    print(f"  recalls = {recalls}")
    print(f"  precisions = {precisions}")
    print(f"  f1_scores = {f1_scores}")
    
    # Check some basic expectations
    assert len(recalls) == 5, f"Expected 5 recall values, got {len(recalls)}"
    # In Python 3, division of integers returns a float, so using / is correct
    # But we can add explicit float() for clarity
    expected_recall = float(10) / 11
    expected_precision = float(10) / 15
    assert abs(recalls[0] - expected_recall) < 1e-10, f"Expected recall[0] to be {expected_recall}, got {recalls[0]}"
    assert abs(precisions[0] - expected_precision) < 1e-10, f"Expected precision[0] to be {expected_precision}, got {precisions[0]}"
    
    print("ROC functions test: PASSED")

def test_chromosome_sorting():
    """Test the chromosome sorting functions."""
    print("\nTesting chromosome sorting...")
    
    # Create a list of chromosomes to sort
    chroms = ["chrX", "chr11", "chr1", "chr20", "chrY", "chrM", "chr2"]
    
    # Expected order: numeric chromosomes in order, then alphabetic
    expected = ["chr1", "chr2", "chr11", "chr20", "chrX", "chrY", "chrM"]
    
    # Sort the chromosomes
    sorted_chroms = sort_chromosomes(chroms)
    
    print(f"  Original: {chroms}")
    print(f"  Sorted: {sorted_chroms}")
    print(f"  Expected: {expected}")
    
    assert sorted_chroms == expected, f"Expected {expected}, got {sorted_chroms}"
    
    print("Chromosome sorting test: PASSED")

def run_all_tests():
    """Run all the tests."""
    # Print module info
    info = get_module_info()
    print("\nModule information:")
    print(f"  Version: {info['version']}")
    print("  Implementation:")
    for name, details in info['modules'].items():
        impl_type = "Cython" if details['is_cython'] else "Pure Python"
        print(f"    {name}: {impl_type} from {details['module']}")
    
    # Run tests
    test_sequence_functions()
    test_variant_processor()
    test_roc_functions()
    test_chromosome_sorting()
    
    print("\nAll tests PASSED!")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Test Cython modernization and Python 3 compatibility")
    parser.add_argument("--verbose", action="store_true", help="Print verbose output")
    
    args = parser.parse_args()
    
    if args.verbose:
        print("Running in verbose mode")
        
    try:
        run_all_tests()
    except Exception as e:
        print(f"Error running tests: {e}")
        sys.exit(1)
        
    print("\nCython modernization tests completed successfully!")