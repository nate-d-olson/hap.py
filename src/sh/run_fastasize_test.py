#!/usr/bin/env python3
# Test module for fastasize.py functionality with Python 3

import logging
import os
import sys
from pathlib import Path

logging.getLogger().setLevel(logging.INFO)

# Use pathlib for more readable path operations
script_dir = Path(__file__).resolve().parent
tools_dir = script_dir.parent / "python" / "Tools"
sys.path.append(str(tools_dir))

from fastasize import calculateLength


def main():
    locations = "chrMT chrY:1-10"
    fastacontiglengths = "{'chrY': 59373566, 'chrM': 16571}"
    total_length = calculateLength(fastacontiglengths, locations)
    if total_length == 10:
        logging.info("fastasize test SUCCEEDED!")
    else:
        logging.error("fastasize test FAILED!")


if __name__ == "__main__":
    main()
