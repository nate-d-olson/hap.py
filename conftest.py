import os
import sys

# Add src/python directory to sys.path for imports during testing
ROOT = os.path.dirname(__file__)
SRC_PYTHON = os.path.join(ROOT, "src", "python")
if os.path.isdir(SRC_PYTHON):
    sys.path.insert(0, SRC_PYTHON)
