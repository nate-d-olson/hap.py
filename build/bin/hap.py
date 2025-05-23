#!/usr/bin/env python3
"""
Wrapper script for hap.py that calls the module directly.
"""
import os
import sys

# Add the project root to the Python path
project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, project_root)

# Import and run the main module
from src.hap_py.hap import main

if __name__ == "__main__":
    main()