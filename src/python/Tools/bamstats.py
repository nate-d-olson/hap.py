# Tools/bamstats.py
"""Mock module to enable testing without external dependencies."""

import logging

import pandas


def bamStats(filename):
    """Mock function to simulate reading stats from a BAM file."""
    logging.warning(f"Mock bamStats called with {filename}")
    return pandas.DataFrame(
        {
            "CHROM": ["1", "2", "X"],
            "COVERAGE": [30.0, 28.5, 15.0],
            "MAPQ": [60.0, 60.0, 60.0],
        }
    )
