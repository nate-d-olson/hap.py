"""Pure Python fallback for the VariantProcessor Cython module."""

from dataclasses import dataclass
from typing import List

@dataclass
class Variant:
    chrom: str
    pos: int
    ref: str
    alt: str

class VariantProcessor:
    """Simplified Python implementation of VariantProcessor."""

    def __init__(self) -> None:
        self.variants: List[Variant] = []

    def add_variant(self, variant: Variant) -> None:
        self.variants.append(variant)

    def get_variant_chrom(self, idx: int) -> str:
        if idx >= len(self.variants):
            raise IndexError(f"Index {idx} out of range")
        return self.variants[idx].chrom

    def process_variants(self, threads: int = 1):
        results = []
        for var in self.variants:
            results.append(
                {
                    "chrom": var.chrom,
                    "position": var.pos,
                    "ref": var.ref,
                    "alt": var.alt,
                    "processed": True,
                }
            )
        return results
