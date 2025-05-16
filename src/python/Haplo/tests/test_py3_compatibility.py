# Testing utilities for Python 3 compatibility
import gc
import os
import tempfile


# Mock classes for testing Python/C++ integration without requiring C++ components
class MockVariant:
    """Mock implementation of C++ Variant for testing."""

    def __init__(self, chrom, pos, ref, alt, qual=None, info=None):
        self.chrom = chrom
        self.pos = pos
        self.ref = ref
        self.alt = alt
        self.qual = qual or []
        self.info = info or {}

    def __repr__(self):
        return f"MockVariant(chrom={self.chrom}, pos={self.pos}, ref={self.ref}, alt={self.alt})"

    def to_dict(self):
        return {
            "chrom": self.chrom,
            "pos": self.pos,
            "ref": self.ref,
            "alt": self.alt,
            "qual": self.qual,
            "info": self.info,
        }


class MockVariantProcessor:
    """Mock implementation of C++ VariantProcessor for testing."""

    def __init__(self):
        self.variants = []

    def add_variant(self, variant):
        """Add a variant to the processor."""
        self.variants.append(variant)

    def process_variants(self, threads=1):
        """Process variants (mock implementation)."""
        results = []
        for var in self.variants:
            result = var.to_dict()
            result["processed"] = True
            results.append(result)
        return results


# Test utilities for Python 3 string handling
def test_string_handling():
    """Test Python 3 string handling."""
    # Create test data
    test_str = "chromosome1"
    test_bytes = b"chromosome1"

    # Test string to bytes conversion
    from_str = test_str.encode("utf8")
    assert isinstance(from_str, bytes)
    assert from_str == test_bytes

    # Test bytes to string conversion
    from_bytes = test_bytes.decode("utf8")
    assert isinstance(from_bytes, str)
    assert from_bytes == test_str

    # Test file I/O with encoding
    with tempfile.NamedTemporaryFile(mode="w", encoding="utf8", delete=False) as f:
        f.write(test_str)
        filename = f.name

    try:
        # Read as text
        with open(filename, encoding="utf8") as f:
            content = f.read()
            assert content == test_str

        # Read as binary
        with open(filename, "rb") as f:
            binary_content = f.read()
            assert binary_content == test_bytes
    finally:
        os.unlink(filename)

    print("String handling tests passed")


# Test utilities for memory management
def test_memory_management():
    """Basic memory management test."""
    # Create a large number of objects
    large_list = []
    for i in range(10000):
        large_list.append(MockVariant(f"chr{i}", i, "A", "T"))

    # Process objects
    processor = MockVariantProcessor()
    for var in large_list:
        processor.add_variant(var)
    results = processor.process_variants()

    # Clean up
    del large_list
    del processor
    del results

    # Force garbage collection
    gc.collect()

    print("Memory management test completed")


# Test chromosome sorting for Python 3
def test_chromosome_sorting():
    """Test chromosome sorting with Python 3 compatibility."""
    from functools import cmp_to_key

    # Define chromosome comparison function
    def cmp_chromosomes(a, b):
        """Compare chromosomes with Python 3 compatibility."""
        # Handle numeric chromosomes
        a_is_num = a.replace("chr", "").isdigit()
        b_is_num = b.replace("chr", "").isdigit()

        if a_is_num and b_is_num:
            # Both numeric, compare as integers
            return int(a.replace("chr", "")) - int(b.replace("chr", ""))
        elif a_is_num:
            # a is numeric, b is not, a comes first
            return -1
        elif b_is_num:
            # b is numeric, a is not, b comes first
            return 1
        else:
            # Both non-numeric, compare as strings
            return -1 if a < b else (1 if a > b else 0)

    # Create test data
    chrom_list = ["chr10", "chr2", "chrX", "chr1", "chrM", "chrY"]

    # Sort with Python 3 compatible method
    sorted_chroms = sorted(chrom_list, key=cmp_to_key(cmp_chromosomes))

    # Expected order: numeric chromosomes in order, then alphabetic
    expected = ["chr1", "chr2", "chr10", "chrM", "chrX", "chrY"]

    assert sorted_chroms == expected, f"Expected {expected}, got {sorted_chroms}"
    print("Chromosome sorting tests passed")


# Run tests if file is executed directly
if __name__ == "__main__":
    print("Running Python 3 compatibility tests")
    test_string_handling()
    test_memory_management()
    test_chromosome_sorting()
    print("All tests passed!")
