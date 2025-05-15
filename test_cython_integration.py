#!/usr/bin/env python3
"""
Test script for Cython integration in hap.py Python 3 migration

This script helps diagnose and verify the integration between Python and C++
components in the Python 3 version of hap.py.

Usage:
    python3 test_cython_integration.py [--mock]
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
    """Test importing Haplo modules"""
    if mock_mode:
        os.environ["HAPLO_USE_MOCK"] = "1"
        logger.info("Using mock implementation")
    else:
        os.environ.pop("HAPLO_USE_MOCK", None)
        logger.info("Attempting to use actual C++ implementation")

    # Test the import
    try:
        import Haplo

        logger.info(f"Successfully imported Haplo version: {Haplo.__version__}")

        # Try to access the Cython module
        try:
            from Haplo.cython import test_module

            result = test_module()
            logger.info(f"Cython module test result: {result}")
            return True, result
        except (ImportError, AttributeError) as e:
            logger.error(f"Failed to import or use Cython module: {e}")
            return False, str(e)

    except ImportError as e:
        logger.error(f"Failed to import Haplo: {e}")
        return False, str(e)


def main():
    """Main entry point"""
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
        # Test only one implementation
        success, result = test_import(mock_mode=args.mock)

        if success:
            logger.info("=== Additional Tests ===")

            # Import main module again
            import Haplo

            # Test additional components if needed
            logger.info("Testing specific functionality...")
            if hasattr(Haplo, "quantify"):
                logger.info("Found quantify module")

            return 0
        else:
            return 1


if __name__ == "__main__":
    sys.exit(main())
