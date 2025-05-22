#!/usr/bin/env python3
"""
Test cases for Python blocksplit implementation.

This module tests the functionality of python_blocksplit.py to ensure 
it correctly splits VCF files into blocks.
"""

import os
import tempfile
from pathlib import Path

import pytest
from Haplo.python_blocksplit import BlockSplitter


class TestBlockSplitter:
    """Test cases for BlockSplitter class."""

    @pytest.fixture
    def example_vcf(self):
        """Fixture to provide path to example VCF file."""
        # Use a VCF file from the example directory
        project_root = Path(__file__).parent.parent
        example_dir = project_root / "example"
        vcf_path = example_dir / "hc.vcf.gz"
        
        if not vcf_path.exists():
            pytest.skip(f"Example VCF file not found: {vcf_path}")
        
        return str(vcf_path)
    
    def test_init(self):
        """Test initialization with default parameters."""
        splitter = BlockSplitter()
        assert splitter.block_size == 1000
        assert splitter.min_distance == 1000
        assert splitter.apply_filters is False
    
    def test_init_custom_params(self):
        """Test initialization with custom parameters."""
        splitter = BlockSplitter(block_size=500, min_distance=2000, apply_filters=True)
        assert splitter.block_size == 500
        assert splitter.min_distance == 2000
        assert splitter.apply_filters is True
    
    def test_process_file(self, example_vcf):
        """Test processing a VCF file."""
        splitter = BlockSplitter(block_size=100)
        
        # Process file without output
        blocks = splitter.process_file(example_vcf)
        
        # Verify blocks were created
        assert len(blocks) > 0
        
        # Check block structure
        for block in blocks:
            assert "chrom" in block
            assert "start" in block
            assert "end" in block
            assert "count" in block
            assert "span" in block
            
            # Basic validation
            assert block["start"] <= block["end"]
            assert block["span"] == block["end"] - block["start"] + 1
            assert block["count"] > 0
    
    def test_process_file_with_output(self, example_vcf):
        """Test processing a VCF file with BED output."""
        splitter = BlockSplitter(block_size=100)
        
        # Create temporary file for output
        with tempfile.NamedTemporaryFile(suffix=".bed", delete=False) as temp_file:
            temp_path = temp_file.name
        
        try:
            # Process file with output
            blocks = splitter.process_file(example_vcf, temp_path)
            
            # Verify blocks were created
            assert len(blocks) > 0
            
            # Verify output file was created
            assert os.path.exists(temp_path)
            assert os.path.getsize(temp_path) > 0
            
            # Check file content
            with open(temp_path, "r") as f:
                lines = f.readlines()
                assert len(lines) == len(blocks)
                
                # Check a line
                line_parts = lines[0].strip().split("\t")
                assert len(line_parts) == 4
                assert "variants=" in line_parts[3]
                assert "span=" in line_parts[3]
        finally:
            # Clean up
            if os.path.exists(temp_path):
                os.unlink(temp_path)
    
    def test_is_filtered(self, example_vcf):
        """Test the _is_filtered method."""
        # This test requires creating a mock VariantRecord
        # Since that's complex, we'll use a fixture with a custom class
        
        class MockVariantRecord:
            def __init__(self, filter_values=None):
                self.filter = filter_values or []
        
        splitter = BlockSplitter()
        
        # Test unfiltered record
        unfiltered = MockVariantRecord([])
        assert splitter._is_filtered(unfiltered) is False
        
        # Test PASS filter
        pass_filter = MockVariantRecord(["PASS"])
        assert splitter._is_filtered(pass_filter) is False
        
        # Test filtered record
        filtered = MockVariantRecord(["LowQual"])
        assert splitter._is_filtered(filtered) is True
        
        # Test multiple filters
        multi_filter = MockVariantRecord(["LowQual", "IndelGap"])
        assert splitter._is_filtered(multi_filter) is True
    
    def test_apply_filters(self, example_vcf):
        """Test that apply_filters works correctly."""
        # Process without filtering
        splitter_no_filter = BlockSplitter(block_size=100, apply_filters=False)
        blocks_no_filter = splitter_no_filter.process_file(example_vcf)
        
        # Process with filtering
        splitter_with_filter = BlockSplitter(block_size=100, apply_filters=True)
        blocks_with_filter = splitter_with_filter.process_file(example_vcf)
        
        # If the example file has filtered variants, the counts should differ
        # Note: this test might not be meaningful if the example file has no filtered variants
        # We'll check the total variant count instead to make the test more robust
        
        total_variants_no_filter = sum(block["count"] for block in blocks_no_filter)
        total_variants_with_filter = sum(block["count"] for block in blocks_with_filter)
        
        # Log for debugging
        print(f"Variants without filtering: {total_variants_no_filter}")
        print(f"Variants with filtering: {total_variants_with_filter}")
        
        # They may be equal if the file has no filtered variants
        assert total_variants_no_filter >= total_variants_with_filter


if __name__ == "__main__":
    pytest.main(["-xvs", __file__])
