"""
Integration tests for the fastasize functionality.
Migrated from src/sh/run_fastasize_test.py
"""

import sys
import pytest
from pathlib import Path


# Add the paths needed for importing the modules
@pytest.fixture(scope="module", autouse=True)
def setup_path():
    script_dir = Path(__file__).resolve().parent.parent.parent  # project root
    tools_dir = script_dir / "src" / "python" / "Tools"
    sys.path.append(str(tools_dir))


@pytest.mark.integration
def test_fastasize_calculation():
    """Test fastasize's calculateLength function"""
    # Import the function after the path has been set up
    from Tools.fastasize import calculateLength

    # Test with the same parameters as in the original script
    locations = "chrMT chrY:1-10"
    fastacontiglengths = "{'chrY': 59373566, 'chrM': 16571}"
    total_length = calculateLength(fastacontiglengths, locations)

    # Verify the result is correct
    assert total_length == 10, f"Expected length 10, got {total_length}"
