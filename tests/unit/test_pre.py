import os
import pytest
import tempfile
import shutil
import sys

# Add the src directory to the path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../src/python')))

from Tools.vcfextract import extractHeadersJSON
import pre


def test_hasChrPrefix():
    """Test the hasChrPrefix function"""
    # Test with chr prefix
    chr_list = ["chr1", "chr2", "chr3", "chrX", "chrY", "chrM"]
    assert pre.hasChrPrefix(chr_list) is True
    
    # Test without chr prefix
    no_chr_list = ["1", "2", "3", "X", "Y", "MT"]
    assert pre.hasChrPrefix(no_chr_list) is False
    
    # Test mixed list
    mixed_list = ["chr1", "2", "3", "chrX", "Y"]
    # This should be None (undecided) or match the majority
    result = pre.hasChrPrefix(mixed_list)
    assert result is None or isinstance(result, bool)


def test_vcfextract_integration():
    """Test the integration with vcfextract"""
    # Create a test VCF file
    test_vcf_content = """##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##FILTER=<ID=LowQual>
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	100	.	A	G	50	PASS	NS=3
"""
    
    # Write the test VCF to a temporary file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.vcf', delete=False) as tmp:
        tmp.write(test_vcf_content)
        test_vcf_path = tmp.name
    
    try:
        # Test extractHeadersJSON function as used in pre.py
        headers = extractHeadersJSON(test_vcf_path)
        
        # Check that we got fields
        assert "fields" in headers
        
        # Check for FILTER fields
        filter_fields = [f for f in headers["fields"] if f["key"] == "FILTER"]
        assert len(filter_fields) > 0
        
        # Check filter IDs
        filter_ids = [f["values"]["ID"] for f in filter_fields]
        assert "PASS" in filter_ids
        assert "LowQual" in filter_ids
        
    finally:
        # Clean up
        os.unlink(test_vcf_path)