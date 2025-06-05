"""Integration tests for GA4GH compliance functionality."""
import os
import pytest
import pandas as pd
from pathlib import Path

from hap_py.quantify.ga4gh import GA4GHFormatter, GA4GHMetrics, GA4GHStratification
from hap_py.quantify import QuantifyEngine


@pytest.fixture
def output_dir(tmp_path):
    """Create a temporary output directory for test results."""
    return tmp_path


@pytest.mark.integration
def test_ga4gh_output_generation(reference_file, rtg_executable, output_dir):
    """Test that GA4GH outputs are generated correctly."""
    # Setup test files
    truth_vcf = Path("example/integration/test.vcf.gz")
    query_vcf = Path("example/integration/test2.vcf.gz")
    
    # Skip if test files don't exist
    if not truth_vcf.exists() or not query_vcf.exists():
        pytest.skip("Test VCF files not available")
    
    # Run hap.py with GA4GH output
    output_prefix = output_dir / "test_ga4gh"
    cmd = [
        "hap.py",
        str(truth_vcf),
        str(query_vcf),
        "-r", str(reference_file),
        "-o", str(output_prefix),
        "--ga4gh",
        "--engine=vcfeval",
        f"--engine-vcfeval-path={rtg_executable}"
    ]
    
    # Execute command and check return code
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0, f"Command failed: {result.stderr}"
    
    # Verify output files exist
    assert (output_prefix.with_suffix(".ga4gh.json")).exists()
    assert (output_prefix.with_suffix(".ga4gh.tsv")).exists()
    
    # Validate JSON content
    with open(output_prefix.with_suffix(".ga4gh.json")) as f:
        data = json.load(f)
        assert "metrics" in data
        assert "snp" in data["metrics"]
        assert "precision" in data["metrics"]["snp"]
        assert "recall" in data["metrics"]["snp"]


@pytest.mark.integration
def test_ga4gh_stratification(reference_file, output_dir):
    """Test GA4GH stratification functionality."""
    # Create stratification object
    stratification = GA4GHStratification()
    
    # Add some test regions
    stratification.add_region("ALL", None)
    stratification.add_region("chr1", "1:1-10000")
    
    # Test stratification logic
    assert stratification.get_region_count() == 2
    assert "ALL" in stratification.get_region_names()
    assert "chr1" in stratification.get_region_names()


@pytest.mark.integration
def test_ga4gh_metrics_calculation():
    """Test calculation of GA4GH metrics."""
    # Create metrics object
    metrics = GA4GHMetrics()
    
    # Add test data
    metrics.add_variant_counts("ALL", "SNP", tp=90, fp=10, fn=10)
    
    # Calculate metrics
    results = metrics.calculate_metrics()
    
    # Verify metrics
    assert "ALL" in results
    assert "SNP" in results["ALL"]
    assert abs(results["ALL"]["SNP"]["precision"] - 0.9) < 0.001
    assert abs(results["ALL"]["SNP"]["recall"] - 0.9) < 0.001
    assert abs(results["ALL"]["SNP"]["f1"] - 0.9) < 0.001
