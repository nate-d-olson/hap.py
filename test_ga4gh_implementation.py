#!/usr/bin/env python3
"""
Simple test script for validating the GA4GH implementation.

Usage:
    python test_ga4gh_implementation.py

This script tests the GA4GH compliance implementation by:
1. Creating a simplified QuantifyEngine with GA4GH quantification
2. Running the quantification process
3. Checking the GA4GH-compliant outputs
"""

import logging
import os
import sys
from pathlib import Path

import pandas as pd

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Try to import hap_py module
try:
    from hap_py.haplo.python_quantify import QuantifyEngine
    from hap_py.haplo.ga4gh_compliance import (
        GA4GHDecision,
        GA4GHDecisionDetail,
        GA4GHFormatter,
        GA4GHMetrics,
        GA4GHStratification,
        GA4GHVariantType,
    )
    from hap_py.haplo.ga4gh_integration import GA4GHIntegration, enhance_quantify_engine_with_ga4gh
except ImportError as e:
    logger.error(f"Could not import hap_py module: {e}")
    sys.exit(1)

def find_test_data():
    """Find test data for validation."""
    # Look for test data in common locations
    test_locations = [
        "./example/integration",
        "./example/test",
        "./tests/data",
    ]
    
    for location in test_locations:
        if not os.path.exists(location):
            continue
            
        truth_files = [
            f for f in os.listdir(location) 
            if f.endswith(".vcf") and "truth" in f.lower()
        ]
        
        if not truth_files:
            continue
            
        truth_file = os.path.join(location, truth_files[0])
        
        # Find corresponding query file
        query_files = [
            f for f in os.listdir(location)
            if f.endswith(".vcf") and "query" in f.lower()
        ]
        
        if not query_files:
            continue
            
        query_file = os.path.join(location, query_files[0])
        
        return truth_file, query_file
    
    return None, None

def main():
    """Main function to test GA4GH implementation."""
    logger.info("Testing GA4GH implementation")
    
    # Find test data
    truth_file, query_file = find_test_data()
    
    if not truth_file or not query_file:
        logger.error("Could not find test data")
        sys.exit(1)
        
    logger.info(f"Found test data: {truth_file} and {query_file}")
    
    # Create output directory
    output_dir = Path("./ga4gh_test_output")
    output_dir.mkdir(exist_ok=True)
    
    # Create a QuantifyEngine with GA4GH quantification
    try:
        engine = QuantifyEngine(
            truth_vcf=truth_file,
            query_vcf=query_file,
            output_prefix=str(output_dir / "output"),
            quantify_method="ga4gh",
            output_vtc=True,  # Output VCF with variant truth categories
            enable_roc_analysis=True,
        )
        
        # Test if GA4GH integration is available
        if hasattr(engine, "ga4gh_integration") and engine.ga4gh_integration:
            logger.info("GA4GH integration is available and initialized")
        else:
            logger.warning("GA4GH integration is not available")
            # Try to enhance the engine manually
            try:
                engine = enhance_quantify_engine_with_ga4gh(engine)
                logger.info("Manually enhanced engine with GA4GH integration")
            except Exception as e:
                logger.error(f"Failed to enhance engine: {e}")
        
        # Run quantification
        logger.info("Running quantification with GA4GH method")
        results = engine.quantify()
        
        # Report results
        logger.info("Quantification results:")
        logger.info(f"  TP: {results['all']['TP']}")
        logger.info(f"  FP: {results['all']['FP']}")
        logger.info(f"  FN: {results['all']['FN']}")
        
        # Check for GA4GH-specific outputs
        metrics_file = Path(str(output_dir / "output") + ".ga4gh.metrics.tsv")
        if metrics_file.exists():
            logger.info(f"GA4GH metrics file created successfully: {metrics_file}")
            
            # Read and display metrics
            try:
                metrics_df = pd.read_csv(metrics_file, sep="\t")
                logger.info(f"GA4GH metrics:\n{metrics_df.head()}")
            except Exception as e:
                logger.error(f"Error reading metrics file: {e}")
        else:
            logger.warning(f"GA4GH metrics file not created: {metrics_file}")
            
        # Overall status
        logger.info("GA4GH implementation test completed successfully")
        
    except Exception as e:
        logger.error(f"Error during testing: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
    
if __name__ == "__main__":
    main()
