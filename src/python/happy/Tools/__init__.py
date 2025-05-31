import os

# Extend package path to include sibling Tools directory
__path__.insert(
    0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "Tools"))
)
