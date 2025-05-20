"""
Unit tests for hap.py Tools module.
"""

import os
import sys

# Add src/python to path for imports during tests
sys.path.insert(
    0,
    os.path.join(
        os.path.dirname(os.path.dirname(os.path.dirname(__file__))),
        "src",
        "python",
    ),
)

import Tools


def test_which():
    """Test the Tools.which function."""
    # Should be able to find 'python' in PATH
    python_path = Tools.which("python")
    assert python_path is not None
    assert os.path.exists(python_path)

    # Should not find a non-existent command
    assert Tools.which("this_command_does_not_exist_12345") is None


def test_defaultReference():
    """Test the defaultReference function."""
    # This is just testing the function returns something or None
    ref = Tools.defaultReference()
    assert isinstance(ref, str) or ref is None
