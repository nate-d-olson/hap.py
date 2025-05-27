#!/usr/bin/env python3
"""
Wrapper script for hap.py main module.
This script provides compatibility with the build/bin/hap.py path expected by tests.
"""

import sys
from pathlib import Path

# Add the src directory to Python path
project_root = Path(__file__).parent.parent.parent
src_dir = project_root / "src"
sys.path.insert(0, str(src_dir))

# Import and run the main hap.py module
from hap_py.hap import main

if __name__ == "__main__":
    sys.exit(main())
