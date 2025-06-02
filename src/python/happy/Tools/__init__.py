import os

# Extend package path to include external Tools directory under src/python/Tools
__path__.insert(
    0,
    os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "Tools")),
)
