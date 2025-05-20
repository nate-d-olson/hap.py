"""Mock Somatic module for testing."""

import logging
import pandas as pd
from typing import Dict, Optional, Union


class FeatureSet:
    """Mock feature set class for somatic variant analysis."""

    sets = {
        "default": "Default feature set",
        "svs": "Structural variants feature set",
        "indels": "Indel variants feature set",
        "snps": "SNP variants feature set",
        "generic": "Generic feature set",  # Added generic set referenced in ftx.py
    }

    @staticmethod
    def make(feature_set_name):
        """Create a feature set."""
        if feature_set_name not in FeatureSet.sets:
            logging.warning(f"Unknown feature set: {feature_set_name}, using default")
            feature_set_name = "default"
        return FeatureSet()

    def __init__(self):
        self.chr_depths = {}

    def setChrDepths(self, depths: Optional[Dict[str, float]]) -> None:
        """Set chromosome depths for normalization."""
        self.chr_depths = depths or {}

    def collect(self, vcf_path: str, label: str = "QUERY") -> pd.DataFrame:
        """Mock collection of features from a VCF file."""
        logging.info(f"Collecting features from {vcf_path} with label {label}")

        # Return sample data that matches the structure expected by ftx.py
        return pd.DataFrame(
            {
                "CHROM": ["1", "2", "X"],
                "POS": [1000, 2000, 3000],
                "FEATURE": ["VAF", "VAF", "VAF"],
                "VALUE": [0.5, 0.3, 0.4],
                "LABEL": [label, label, label],
            }
        )
