#!/usr/bin/env python3
"""
This module provides a mock implementation of the Haplo.cython._internal module.
It can be used for testing Python code that depends on the C++ functionality
without requiring the actual C++ components to be built.

This is especially useful during the Python 2 to 3 migration process.
"""

from typing import Any, Dict, Union


class MockVariantRecord:
    """Mock implementation of a variant record"""

    def __init__(self, chrom: str, pos: int, ref: str, alt: str, qual: float = 0):
        self.chrom = chrom
        self.pos = pos
        self.ref = ref
        self.alt = alt
        self.qual = qual
        self.info = {}
        self.format = {}

    def __str__(self) -> str:
        return f"{self.chrom}:{self.pos} {self.ref}>{self.alt}"


class MockHaploCompare:
    """Mock implementation of the haplotype comparison functionality"""

    def __init__(self, reference_file: str = None):
        self.reference_file = reference_file
        self.truth_variants = []
        self.query_variants = []
        self.results = None

    def add_truth_variant(self, variant: Union[MockVariantRecord, Dict]) -> None:
        """Add a truth variant"""
        if isinstance(variant, dict):
            variant = MockVariantRecord(**variant)
        self.truth_variants.append(variant)

    def add_query_variant(self, variant: Union[MockVariantRecord, Dict]) -> None:
        """Add a query variant"""
        if isinstance(variant, dict):
            variant = MockVariantRecord(**variant)
        self.query_variants.append(variant)

    def compare(self) -> Dict[str, Any]:
        """Perform the comparison"""
        # Simple mock implementation
        self.results = {
            "total_truth": len(self.truth_variants),
            "total_query": len(self.query_variants),
            "tp": 0,
            "fp": 0,
            "fn": 0,
        }

        # Count exact matches as TP
        for qvar in self.query_variants:
            matched = False
            for tvar in self.truth_variants:
                if (
                    qvar.chrom == tvar.chrom
                    and qvar.pos == tvar.pos
                    and qvar.ref == tvar.ref
                    and qvar.alt == tvar.alt
                ):
                    matched = True
                    self.results["tp"] += 1
                    break
            if not matched:
                self.results["fp"] += 1

        # Count unmatched truth variants as FN
        for tvar in self.truth_variants:
            matched = False
            for qvar in self.query_variants:
                if (
                    qvar.chrom == tvar.chrom
                    and qvar.pos == tvar.pos
                    and qvar.ref == tvar.ref
                    and qvar.alt == tvar.alt
                ):
                    matched = True
                    break
            if not matched:
                self.results["fn"] += 1

        # Calculate metrics
        if self.results["tp"] + self.results["fn"] > 0:
            self.results["recall"] = self.results["tp"] / (
                self.results["tp"] + self.results["fn"]
            )
        else:
            self.results["recall"] = 0

        if self.results["tp"] + self.results["fp"] > 0:
            self.results["precision"] = self.results["tp"] / (
                self.results["tp"] + self.results["fp"]
            )
        else:
            self.results["precision"] = 0

        if self.results["precision"] + self.results["recall"] > 0:
            self.results["f1"] = (
                2
                * (self.results["precision"] * self.results["recall"])
                / (self.results["precision"] + self.results["recall"])
            )
        else:
            self.results["f1"] = 0

        return self.results


# Mock version functions
def get_version() -> str:
    """Get the hap.py version string (mock implementation)"""
    return "0.3.15-mock"


def get_build_time() -> str:
    """Get the hap.py build timestamp (mock implementation)"""
    return "2023-01-01 00:00:00"


def test_module() -> Dict[str, str]:
    """Test if the module is working properly (mock implementation)"""
    return {"version": get_version(), "build_time": get_build_time()}