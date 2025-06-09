#!/usr/bin/env python3
"""
Entry point script for hap.py CLI tool.
This provides a convenient way to run hap.py from the root directory.
"""

import sys
from pathlib import Path

# Add the source directory to Python path
src_dir = Path(__file__).parent / "src"
sys.path.insert(0, str(src_dir))

# Import and run the main CLI
from hap_py.hap import main

if __name__ == "__main__":
    sys.exit(main())
