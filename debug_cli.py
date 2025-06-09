#!/usr/bin/env python3
"""
Debug script to find the argparse formatting issue
"""
import sys
sys.path.insert(0, 'src')

try:
    from hap_py.hap import _setup_imports
    # Import all necessary modules
    pre, qfy, get_rtg_path, gvcf2bed, vcfeval, bcftools, vcfextract, bedOverlapCheck, fastaContigLengths, getPool, sessionInfo, version = _setup_imports()
    print("Imports successful")
except Exception as e:
    print(f"Import error: {e}")
    import traceback
    traceback.print_exc()

import argparse

# Test creating the argument parser step by step
try:
    parser = argparse.ArgumentParser("Haplotype Comparison")
    print("Basic parser created")
    
    # Add version argument
    parser.add_argument(
        "-v",
        "--version",
        dest="version",
        action="store_true",
        help="Show version number and exit.",
    )
    print("Version argument added")
    
    # Test parsing version
    try:
        args = parser.parse_args(['--version'])
        print("Version parsing successful")
    except Exception as e:
        print(f"Version parsing error: {e}")
    
    # Test help
    try:
        args = parser.parse_args(['--help'])
    except SystemExit:
        print("Help displayed successfully")
    except Exception as e:
        print(f"Help error: {e}")
        
except Exception as e:
    print(f"Parser creation error: {e}")
    import traceback
    traceback.print_exc()
