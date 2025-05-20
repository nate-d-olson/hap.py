import os
import json
import tempfile
import pytest
from Tools.vcfextract import extractHeadersJSON, extract_header


def test_extractHeadersJSON_returns_data():
    """Test that extractHeadersJSON returns header data when called with only VCF path"""
    # Use a sample VCF file
    vcf_path = os.path.join(os.path.dirname(__file__), "test_filter.vcf")
    
    # Call function with only VCF path
    headers = extractHeadersJSON(vcf_path)
    
    # Check that fields list is included for compatibility
    assert "fields" in headers
    
    # Check that we got basic header structure
    assert isinstance(headers, dict)


def test_extractHeadersJSON_with_outfile():
    """Test that extractHeadersJSON writes to file when outfile is provided"""
    # Use a sample VCF file  
    vcf_path = os.path.join(os.path.dirname(__file__), "test_filter.vcf")
    
    # Create temp output file
    with tempfile.NamedTemporaryFile(suffix=".json", delete=False) as tmp:
        out_path = tmp.name
    
    try:
        # Call function with outfile path
        result = extractHeadersJSON(vcf_path, out_path)
        
        # Check that it returned data
        assert isinstance(result, dict)
        
        # Check that file was written
        assert os.path.exists(out_path)
        
        # Check file content
        with open(out_path, "r") as f:
            file_content = json.load(f)
        
        # Should match returned data
        assert file_content == result
        
    finally:
        # Clean up
        if os.path.exists(out_path):
            os.unlink(out_path)


def test_extract_header_filters():
    """Test that extract_header extracts FILTER fields when requested"""
    # Use a sample VCF file with filters
    vcf_path = os.path.join(os.path.dirname(__file__), "test_filter.vcf")
    
    # Extract with filters enabled
    headers = extract_header(vcf_path, extract_filters=True)
    
    # Check that fields includes filter data
    filter_fields = [f for f in headers["fields"] if f["key"] == "FILTER"]
    assert len(filter_fields) > 0
    
    # Check for specific filters
    filter_ids = [f["values"]["ID"] for f in filter_fields]
    assert "PASS" in filter_ids
    assert "LowQual" in filter_ids
    assert "q10" in filter_ids