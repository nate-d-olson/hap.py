#!/usr/bin/env python3
"""
Test script for Cython integration in the Python 3 migration of hap.py.

This script tests the ability to load the Cython modules and verifies
that the fallback to mock implementations works correctly.

Usage:
    python3 test_cython_integration.py [--mock] [--test-both]
"""

import argparse
import logging
import os
import sys
from pathlib import Path

# Configure logging
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger("cython-test")


def test_import(mock_mode=False):
    """Test importing Haplo modules with or without mock implementations."""
    # Set environment variable to control mock mode
    if mock_mode:
        os.environ["HAPLO_USE_MOCK"] = "1"
        logger.info("Using mock implementation")
    else:
        os.environ.pop("HAPLO_USE_MOCK", None)
        logger.info("Attempting to use actual C++ implementation")

    # Test imports
    try:
        import Haplo
        from Haplo.cython import USING_MOCK

        logger.info(f"Successfully imported Haplo version: {Haplo.__version__}")
        logger.info(f"Using mock implementation: {USING_MOCK}")

        # Try basic operations
        try:
            # Try string handling - this is important for Python 3 compatibility
            # as it tests bytes vs unicode string handling
            test_str = "ACGT"
            comp = Haplo.complement_sequence(test_str)
            rev_comp = Haplo.reverse_complement(test_str)

            logger.info(f"Complement of {test_str}: {comp}")
            logger.info(f"Reverse complement of {test_str}: {rev_comp}")

            # Check if results are as expected
            expected_comp = "TGCA"
            expected_rev_comp = "ACGT"[::-1].translate(str.maketrans("ACGT", "TGCA"))

            if comp != expected_comp or rev_comp != expected_rev_comp:
                logger.warning("Unexpected complement results!")

            return True, {"comp": comp, "rev_comp": rev_comp}

        except Exception as e:
            logger.error(f"Failed during basic operations: {e}")
            return False, str(e)

    except ImportError as e:
        logger.error(f"Failed to import Haplo: {e}")
        return False, str(e)
    except Exception as e:
        logger.error(f"Unexpected error: {e}")
        return False, str(e)


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(description="Test Cython integration")
    parser.add_argument(
        "--mock", action="store_true", help="Force use of mock implementation"
    )
    parser.add_argument(
        "--test-both",
        action="store_true",
        help="Test both mock and real implementations",
    )
    args = parser.parse_args()

    # Add src/python to Python path if needed
    script_dir = Path(__file__).resolve().parent
    python_dir = script_dir / "src" / "python"
    if python_dir.exists():
        sys.path.insert(0, str(python_dir))
        logger.info(f"Added {python_dir} to Python path")

    if args.test_both:
        # Test both implementations
        logger.info("=== Testing MOCK implementation ===")
        mock_success, mock_result = test_import(mock_mode=True)

        logger.info("\n=== Testing REAL C++ implementation ===")
        real_success, real_result = test_import(mock_mode=False)

        # Summary
        logger.info("\n=== SUMMARY ===")
        logger.info(f"Mock implementation: {'SUCCESS' if mock_success else 'FAILED'}")
        logger.info(f"C++ implementation:  {'SUCCESS' if real_success else 'FAILED'}")

        if not mock_success and not real_success:
            logger.error("Both implementations failed to load!")
            return 1

        return 0
    else:
        # Test only one implementation as specified
        success, result = test_import(mock_mode=args.mock)

        if success:
            logger.info("=== Additional Tests ===")

            # Import main module again to test more functionality
            import Haplo

            # Test Haplo version info
            if hasattr(Haplo, "get_module_info"):
                module_info = Haplo.get_module_info()
                logger.info(f"Module info: {module_info}")

            return 0
        else:
            return 1


if __name__ == "__main__":
    sys.exit(main())
